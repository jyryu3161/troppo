"""
TINITReconstructor - One-stop pipeline for tINIT-based tissue-specific model reconstruction.

Usage:
    from troppo import TINITReconstructor

    recon = TINITReconstructor(
        model_path="Human-GEM.xml",
        expression_path="expression_data.csv",
        solver="CPLEX"
    )
    results = recon.run()
    results.summary()
    results.plot_overview(save_dir="figures/")
"""

import warnings
from typing import Optional, List, Callable, Union

import numpy as np
import pandas as pd


class TINITReconstructor:
    """tINIT-based tissue-specific metabolic model reconstruction pipeline.

    Parameters
    ----------
    model_path : str
        Path to SBML or JSON metabolic model file.
    expression_path : str, optional
        Path to CSV/TSV expression data file.
    expression_data : pd.DataFrame, optional
        Expression data as DataFrame (genes as index, samples as columns).
        Provide either expression_path or expression_data.
    gene_id_type : str
        Gene ID type in expression data. 'auto' for auto-detection.
        Options: 'auto', 'symbol', 'entrez_id', 'ensembl_gene_id', 'hgnc_id'
    target_id_type : str
        Gene ID type in the model. 'auto' for auto-detection.
    solver : str
        Solver to use: 'CPLEX' or 'GUROBI'.
    n_jobs : int
        Number of parallel jobs. Currently only n_jobs=1 is supported
        (tINIT solver calls are not thread-safe).
    verbose : bool
        Print progress information.
    """

    def __init__(
        self,
        model_path: str,
        expression_path: str = None,
        expression_data: pd.DataFrame = None,
        gene_id_type: str = "auto",
        target_id_type: str = "auto",
        solver: str = "CPLEX",
        n_jobs: int = 1,
        verbose: bool = True,
    ):
        self.model_path = model_path
        self.expression_path = expression_path
        self.solver = solver.upper()
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.gene_id_type = gene_id_type
        self.target_id_type = target_id_type

        # Configuration defaults
        self._production_weight = 0.5
        self._allow_excretion = False
        self._no_reverse_loops = True
        self._essential_reactions = None
        self._protected_subsystems = None
        self._biomass_reaction = "auto"
        self._and_func = min
        self._or_func = max
        self._score_integration = "adjusted"
        self._custom_integration_func = None

        # Load data
        self._model = None
        self._expression_df = None
        self._model_wrapper = None

        self._load_model()
        self._load_expression_data(expression_data)

    def _load_model(self):
        """Load COBRA model."""
        from troppo.io import load_model
        if self.verbose:
            print(f"Loading model from {self.model_path}...")
        self._model = load_model(self.model_path, biomass_reaction=self._biomass_reaction)
        if self.verbose:
            print(f"  Model: {len(self._model.reactions)} reactions, "
                  f"{len(self._model.metabolites)} metabolites, "
                  f"{len(self._model.genes)} genes")

    def _load_expression_data(self, expression_data: pd.DataFrame = None):
        """Load expression data from file or DataFrame."""
        if expression_data is not None:
            self._expression_df = expression_data.copy()
        elif self.expression_path is not None:
            from troppo.io import load_expression_csv
            if self.verbose:
                print(f"Loading expression data from {self.expression_path}...")
            self._expression_df = load_expression_csv(
                self.expression_path,
                gene_id_type=self.gene_id_type,
            )
        else:
            self._expression_df = None
            return

        if self.verbose and self._expression_df is not None:
            print(f"  Expression data: {self._expression_df.shape[0]} genes, "
                  f"{self._expression_df.shape[1]} samples")

    def configure(
        self,
        production_weight: float = None,
        allow_excretion: bool = None,
        no_reverse_loops: bool = None,
        essential_reactions: list = None,
        protected_subsystems: list = None,
        biomass_reaction: str = None,
        and_func: Callable = None,
        or_func: Callable = None,
        score_integration: str = None,
        custom_integration_func: Callable = None,
    ) -> 'TINITReconstructor':
        """Configure reconstruction parameters.

        Parameters
        ----------
        production_weight : float
            Weight for production reactions (default: 0.5).
        allow_excretion : bool
            Allow excretion in the model (default: False).
        no_reverse_loops : bool
            Prevent reverse loops (default: True).
        essential_reactions : list
            Reaction IDs to always include.
        protected_subsystems : list
            Subsystem names whose reactions are protected.
        biomass_reaction : str
            Biomass reaction ID ('auto' for auto-detection).
        and_func : callable
            Function for AND in GPR rules (default: min).
        or_func : callable
            Function for OR in GPR rules (default: max).
        score_integration : str
            Integration strategy: 'adjusted', 'continuous', or 'custom'.
        custom_integration_func : callable
            Custom integration function (when score_integration='custom').

        Returns
        -------
        TINITReconstructor
            Self, for method chaining.
        """
        if production_weight is not None:
            self._production_weight = production_weight
        if allow_excretion is not None:
            self._allow_excretion = allow_excretion
        if no_reverse_loops is not None:
            self._no_reverse_loops = no_reverse_loops
        if essential_reactions is not None:
            self._essential_reactions = essential_reactions
        if protected_subsystems is not None:
            self._protected_subsystems = protected_subsystems
        if biomass_reaction is not None:
            self._biomass_reaction = biomass_reaction
            self._load_model()
        if and_func is not None:
            self._and_func = and_func
        if or_func is not None:
            self._or_func = or_func
        if score_integration is not None:
            self._score_integration = score_integration
        if custom_integration_func is not None:
            self._custom_integration_func = custom_integration_func
        return self

    def preview_data(self) -> dict:
        """Preview data summary before running reconstruction.

        Returns
        -------
        dict
            Summary with gene count, sample count, ID match rate, etc.
        """
        info = {
            'model_reactions': len(self._model.reactions),
            'model_metabolites': len(self._model.metabolites),
            'model_genes': len(self._model.genes),
        }

        if self._expression_df is not None:
            info['expression_genes'] = self._expression_df.shape[0]
            info['expression_samples'] = self._expression_df.shape[1]
            info['sample_names'] = list(self._expression_df.columns)

            # Check gene ID overlap
            model_gene_ids = {g.id for g in self._model.genes}
            expr_gene_ids = set(self._expression_df.index.astype(str))
            overlap = model_gene_ids & expr_gene_ids
            info['gene_overlap'] = len(overlap)
            info['gene_match_rate'] = round(len(overlap) / max(len(expr_gene_ids), 1) * 100, 1)
        else:
            info['expression_data'] = 'Not loaded'

        return info

    def get_model_info(self) -> dict:
        """Get model information.

        Returns
        -------
        dict
            Model details (reactions, genes, subsystems).
        """
        subsystems = set()
        for rxn in self._model.reactions:
            sub = getattr(rxn, 'subsystem', '') or 'Unknown'
            subsystems.add(sub)

        return {
            'reactions': len(self._model.reactions),
            'metabolites': len(self._model.metabolites),
            'genes': len(self._model.genes),
            'subsystems': sorted(subsystems),
            'n_subsystems': len(subsystems),
            'objective': str(self._model.objective),
        }

    def _get_essential_reaction_indices(self) -> list:
        """Get indices of essential reactions (from IDs + protected subsystems)."""
        r_ids = list(self._model_wrapper.model_reader.r_ids)
        essential_idx = []

        if self._essential_reactions:
            for rid in self._essential_reactions:
                if rid in r_ids:
                    essential_idx.append(r_ids.index(rid))

        if self._protected_subsystems:
            for rxn in self._model.reactions:
                sub = getattr(rxn, 'subsystem', '') or ''
                if sub in self._protected_subsystems and rxn.id in r_ids:
                    idx = r_ids.index(rxn.id)
                    if idx not in essential_idx:
                        essential_idx.append(idx)

        return essential_idx

    def _prepare_integration_strategy(self):
        """Prepare the integration strategy object."""
        from troppo.omics.integration import (
            ContinuousScoreIntegrationStrategy,
            AdjustedScoreIntegrationStrategy,
            CustomSelectionIntegrationStrategy,
        )

        if self._score_integration == 'continuous':
            return ContinuousScoreIntegrationStrategy(self._custom_integration_func)
        elif self._score_integration == 'adjusted':
            # Collect protected reaction IDs
            protected = list(self._essential_reactions or [])
            if self._protected_subsystems:
                for rxn in self._model.reactions:
                    sub = getattr(rxn, 'subsystem', '') or ''
                    if sub in self._protected_subsystems:
                        protected.append(rxn.id)
            return AdjustedScoreIntegrationStrategy(protected)
        elif self._score_integration == 'custom' and self._custom_integration_func:
            return CustomSelectionIntegrationStrategy({'custom': self._custom_integration_func})
        else:
            return ContinuousScoreIntegrationStrategy()

    def _convert_gene_ids_if_needed(self, expression_dict: dict) -> dict:
        """Convert gene IDs in expression data to match model if needed."""
        from troppo.io import detect_gene_id_type, convert_gene_ids

        model_gene_ids = {g.id for g in self._model.genes}
        expr_gene_ids = set(str(k) for k in expression_dict.keys())

        overlap = model_gene_ids & expr_gene_ids
        match_rate = len(overlap) / max(len(expr_gene_ids), 1)

        if match_rate > 0.3:
            if self.verbose:
                print(f"  Gene ID match rate: {match_rate*100:.1f}% ({len(overlap)}/{len(expr_gene_ids)}) - using as-is")
            return expression_dict

        # Try conversion
        if self.verbose:
            print(f"  Low gene ID match rate ({match_rate*100:.1f}%), attempting ID conversion...")

        from_type = self.gene_id_type
        if from_type == "auto":
            from_type = detect_gene_id_type(list(expr_gene_ids)[:50])
            if self.verbose:
                print(f"  Detected expression gene ID type: {from_type}")

        to_type = self.target_id_type
        if to_type == "auto":
            to_type = detect_gene_id_type(list(model_gene_ids)[:50])
            if self.verbose:
                print(f"  Detected model gene ID type: {to_type}")

        if from_type == to_type:
            return expression_dict

        try:
            id_map = convert_gene_ids(
                list(expression_dict.keys()),
                from_type=from_type,
                to_type=to_type,
                report=self.verbose,
            )
            converted = {}
            for old_id, value in expression_dict.items():
                new_id = id_map.get(old_id) or id_map.get(str(old_id))
                if new_id:
                    converted[new_id] = value
            if self.verbose:
                print(f"  Converted {len(converted)}/{len(expression_dict)} gene IDs")
            return converted if converted else expression_dict
        except Exception as e:
            warnings.warn(f"Gene ID conversion failed: {e}")
            return expression_dict

    def _load_tasks(self, task_file: str) -> list:
        """task 파일을 로드한다. .xlsx/.xls → ExcelTaskIO, .json → JSONTaskIO.

        Parameters
        ----------
        task_file : str
            Task 파일 경로.

        Returns
        -------
        list of Task
        """
        from troppo.tasks.task_io import ExcelTaskIO, JSONTaskIO

        ext = task_file.lower().rsplit('.', 1)[-1]
        if ext in ('xlsx', 'xls'):
            tasks = ExcelTaskIO().read_task(task_file)
        elif ext == 'json':
            tasks = JSONTaskIO().read_task(task_file)
        else:
            raise ValueError(
                f"Unsupported task file format: .{ext} (supported: .xlsx, .xls, .json)"
            )
        if self.verbose:
            print(f"  Loaded {len(tasks)} tasks from {task_file}")
        return tasks

    def _get_medium_reaction_ids(self) -> set:
        """모델에서 열려있는 exchange 반응 ID를 반환한다.

        Exchange 반응 판별 기준:
        - 반응 ID가 'EX_'로 시작하거나,
        - 단일 metabolite 반응이면서 lower_bound < 0 (uptake 가능)

        Returns
        -------
        set of str
            Medium에 포함된 반응 ID 집합.
        """
        medium_rxn_ids = set()
        for rxn in self._model.reactions:
            is_exchange = rxn.id.startswith('EX_') or (
                len(rxn.metabolites) == 1 and rxn.lower_bound < 0
            )
            if is_exchange and rxn.lower_bound < 0:
                medium_rxn_ids.add(rxn.id)
        if self.verbose:
            print(f"  Identified {len(medium_rxn_ids)} open exchange reactions (medium)")
        return medium_rxn_ids

    def _identify_medium_related_tasks(self, tasks: list, medium_rxn_ids: set) -> set:
        """배지 조건 반응과 관련된 task 이름을 반환한다.

        판별 기준:
        - task.inflow_dict 의 metabolite ID가 medium exchange 반응의 metabolite와 매칭되거나
        - task.flux_constraints 에 medium exchange 반응 ID가 포함된 경우

        Parameters
        ----------
        tasks : list of Task
        medium_rxn_ids : set of str

        Returns
        -------
        set of str
            배지-연관 task 이름 집합.
        """
        medium_metabolites = set()
        for rxn in self._model.reactions:
            if rxn.id in medium_rxn_ids:
                for met in rxn.metabolites:
                    medium_metabolites.add(met.id)
                    # compartment 제거 버전도 추가 (e.g. 'glc__D_e' → 'glc__D')
                    if '_' in met.id:
                        medium_metabolites.add(met.id.rsplit('_', 1)[0])

        medium_related = set()
        for task in tasks:
            inflow_mets = set(task.inflow_dict.keys())
            if inflow_mets & medium_metabolites:
                medium_related.add(task.name)
                continue
            if set(task.flux_constraints.keys()) & medium_rxn_ids:
                medium_related.add(task.name)

        if self.verbose:
            print(f"  Medium-related tasks: {len(medium_related)}/{len(tasks)}")
        return medium_related

    def _build_bound_changes(self, result_dict: dict) -> dict:
        """비활성 반응을 (0, 0)으로 설정한 bound_changes 딕셔너리를 생성한다.

        TaskEvaluator.batch_evaluate()의 bound_changes 리스트의 각 원소 형태.
        활성 반응은 포함하지 않아 원래 bounds를 유지한다.

        Parameters
        ----------
        result_dict : dict
            {reaction_id: bool} - tINIT reconstruction 결과.

        Returns
        -------
        dict
            {reaction_id: (0.0, 0.0)} for inactive reactions.
        """
        return {
            rxn.id: (0.0, 0.0)
            for rxn in self._model.reactions
            if not result_dict.get(rxn.id, False)
        }

    def _validate_tasks_for_sample(
        self,
        sample_name: str,
        result_dict: dict,
        tasks: list,
        medium_related_names: set,
    ):
        """단일 샘플에 대한 metabolic task 검증을 수행한다.

        Parameters
        ----------
        sample_name : str
        result_dict : dict
            {reaction_id: bool} - tINIT reconstruction 결과.
        tasks : list of Task
        medium_related_names : set of str
            배지-연관 task 이름 집합.

        Returns
        -------
        TaskValidationResult
        """
        from troppo.tasks.core import TaskEvaluator, TaskValidationResult

        if self.verbose:
            print(f"    Validating {len(tasks)} tasks for {sample_name}...")

        tev = TaskEvaluator(model=self._model, tasks=tasks)
        bound_changes = [self._build_bound_changes(result_dict)]

        raw_results = tev.batch_evaluate(bound_changes=bound_changes, threads=1)

        results = {
            task_name: truth
            for (_, task_name), (truth, _, _) in raw_results[0].items()
        }

        return TaskValidationResult(
            sample_name=sample_name,
            results=results,
            medium_related=medium_related_names,
        )

    def _print_task_report(self, tvr) -> None:
        """콘솔에 task 검증 요약 리포트를 출력한다.

        Parameters
        ----------
        tvr : TaskValidationResult
        """
        sep = "─" * 52
        print(f"\n    {sep}")
        print(f"    Task Validation: {tvr.sample_name}")
        print(f"    {sep}")
        print(f"    Tasks evaluated : {tvr.total}")
        if tvr.total > 0:
            print(f"    Passed          : {tvr.passed}/{tvr.total} "
                  f"({tvr.passed / tvr.total * 100:.1f}%)")
        if tvr.medium_total > 0:
            print(f"    Medium-related  : {tvr.medium_passed}/{tvr.medium_total} "
                  f"({tvr.medium_passed / tvr.medium_total * 100:.1f}%)")
        if tvr.failed_tasks:
            print(f"    Failed tasks:")
            for task_name in sorted(tvr.failed_tasks):
                tag = " [MEDIUM-RELATED]" if task_name in tvr.medium_related else ""
                print(f"      - {task_name}{tag}")
        print(f"    {sep}")

    def _run_single_sample(self, sample_name: str, expression_dict: dict) -> dict:
        """Run tINIT reconstruction for a single sample.

        Returns
        -------
        dict
            {reaction_id: True/False} for all reactions in the model.
        """
        from troppo.methods_wrappers import ReconstructionWrapper
        from troppo.omics.core import OmicsContainer

        if self.verbose:
            print(f"  Reconstructing: {sample_name}...")

        # Convert IDs if needed
        converted_data = self._convert_gene_ids_if_needed(expression_dict)

        # Create OmicsContainer
        oc = OmicsContainer(
            omicstype='transcriptomics',
            condition=sample_name,
            data=converted_data,
            nomenclature=self.gene_id_type if self.gene_id_type != 'auto' else None,
        )

        # Build model wrapper
        wrapper = ReconstructionWrapper(self._model)
        self._model_wrapper = wrapper

        # Get essential reaction indices
        essential_idx = self._get_essential_reaction_indices()

        # Build integration strategy
        strat = self._prepare_integration_strategy()

        # Run reconstruction via wrapper
        kwargs = {
            'solver': self.solver,
            'production_weight': self._production_weight,
            'allow_excretion': self._allow_excretion,
            'no_reverse_loops': self._no_reverse_loops,
            'essential_reactions': essential_idx,
        }

        try:
            result = wrapper.run_from_omics(
                omics_data=oc,
                integration_strategy=strat,
                and_or_funcs=(self._and_func, self._or_func),
                raise_errors=True,
                **kwargs,
            )
            return result
        except Exception as e:
            warnings.warn(f"Reconstruction failed for sample {sample_name}: {e}")
            return {r.id: False for r in self._model.reactions}

    def run(
        self,
        samples: list = None,
        gap_fill: bool = True,
        validate_tasks: bool = False,
        task_file: str = None,
    ) -> 'ReconstructionResults':
        """Run tINIT reconstruction for all (or selected) samples.

        Parameters
        ----------
        samples : list, optional
            Sample names to reconstruct. None = all samples.
        gap_fill : bool
            Whether to apply gap-filling after reconstruction.
        validate_tasks : bool
            Whether to validate metabolic tasks.
        task_file : str, optional
            Path to metabolic task Excel file.

        Returns
        -------
        ReconstructionResults
            Results container with reaction matrix and analysis methods.
        """
        from troppo.results import ReconstructionResults

        if self._expression_df is None:
            raise ValueError("No expression data loaded. Provide expression_path or expression_data.")

        # Select samples
        all_samples = list(self._expression_df.columns)
        if samples is not None:
            selected = [s for s in samples if s in all_samples]
            if len(selected) < len(samples):
                missing = set(samples) - set(selected)
                warnings.warn(f"Samples not found in expression data: {missing}")
        else:
            selected = all_samples

        if validate_tasks and not task_file:
            raise ValueError(
                "task_file must be provided when validate_tasks=True. "
                "Example: run(validate_tasks=True, task_file='metabolicTasks_Essential.xlsx')"
            )

        if self.verbose:
            print(f"\nRunning tINIT reconstruction for {len(selected)} samples "
                  f"(solver={self.solver})...")

        # Prepare tasks if validation is requested
        tasks = None
        medium_related_names = set()
        if validate_tasks:
            tasks = self._load_tasks(task_file)
            medium_rxn_ids = self._get_medium_reaction_ids()
            medium_related_names = self._identify_medium_related_tasks(tasks, medium_rxn_ids)
            if not medium_related_names and self.verbose:
                warnings.warn(
                    "No medium-related tasks identified. "
                    "Check that task inflow metabolites match model exchange reaction metabolites."
                )

        # Run per sample
        results_dict = {}
        task_results = {}
        for i, sample in enumerate(selected):
            expression_dict = self._expression_df[sample].dropna().to_dict()
            # Ensure keys are strings
            expression_dict = {str(k): float(v) for k, v in expression_dict.items()}
            result = self._run_single_sample(sample, expression_dict)
            results_dict[sample] = result

            if self.verbose:
                n_active = sum(1 for v in result.values() if v)
                print(f"    [{i+1}/{len(selected)}] {sample}: {n_active} active reactions")

            # Task validation immediately after reconstruction
            if validate_tasks and tasks:
                tvr = self._validate_tasks_for_sample(
                    sample, result, tasks, medium_related_names
                )
                task_results[sample] = tvr
                if self.verbose:
                    self._print_task_report(tvr)

        # Build metadata
        metadata = {
            'solver': self.solver,
            'production_weight': self._production_weight,
            'allow_excretion': self._allow_excretion,
            'no_reverse_loops': self._no_reverse_loops,
            'score_integration': self._score_integration,
            'model_path': self.model_path,
            'expression_path': self.expression_path,
        }

        results = ReconstructionResults(
            results_dict=results_dict,
            model=self._model,
            metadata=metadata,
            task_results=task_results if task_results else None,
        )

        if self.verbose:
            summary = results.summary()
            print(f"\nReconstruction complete:")
            print(f"  Samples: {results.n_samples}")
            print(f"  Mean active reactions: {summary['n_active_reactions'].mean():.0f}")
            print(f"  Common reactions: {len(results.get_common_reactions())}")

        return results

    def __repr__(self):
        expr_info = ""
        if self._expression_df is not None:
            expr_info = f", expression=({self._expression_df.shape[0]} genes, {self._expression_df.shape[1]} samples)"
        return (f"TINITReconstructor(model={self.model_path}, solver={self.solver}"
                f"{expr_info})")
