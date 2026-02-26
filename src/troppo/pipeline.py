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
        # Use model-side coverage: are enough model genes covered by expression data?
        model_match_rate = len(overlap) / max(len(model_gene_ids), 1)
        expr_match_rate = len(overlap) / max(len(expr_gene_ids), 1)
        match_rate = max(model_match_rate, expr_match_rate)

        if match_rate > 0.3:
            if self.verbose:
                print(f"  Gene ID match rate: {model_match_rate*100:.1f}% of model genes "
                      f"({len(overlap)}/{len(model_gene_ids)}) - using as-is")
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
        """Load a task file. .xlsx/.xls → ExcelTaskIO, .json → JSONTaskIO.

        Excel task files may use 'Name[compartment]' metabolite IDs,
        which are automatically converted to the model's internal IDs.

        Parameters
        ----------
        task_file : str
            Path to the task file.

        Returns
        -------
        list of Task
        """
        from troppo.tasks.task_io import ExcelTaskIO, JSONTaskIO

        ext = task_file.lower().rsplit('.', 1)[-1]
        if ext in ('xlsx', 'xls'):
            tasks = ExcelTaskIO().read_task(task_file)
            tasks = self._translate_task_metabolite_ids(tasks)
        elif ext == 'json':
            tasks = JSONTaskIO().read_task(task_file)
        else:
            raise ValueError(
                f"Unsupported task file format: .{ext} (supported: .xlsx, .xls, .json)"
            )
        if self.verbose:
            print(f"  Loaded {len(tasks)} tasks from {task_file}")
        return tasks

    def _translate_task_metabolite_ids(self, tasks: list) -> list:
        """Translate 'Name[compartment]' metabolite IDs in Excel task files to model IDs.

        Models such as Human-GEM use internal HMDB IDs (e.g. m02040s), while
        Excel task files use human-readable names like 'H2O[s]'. This method
        builds a mapping from metabolite name + compartment and translates IDs
        automatically. Untranslatable IDs are left unchanged.

        Parameters
        ----------
        tasks : list of Task

        Returns
        -------
        list of Task
            Tasks with metabolite IDs translated to model format.
        """
        # Build mapping: "name[compartment]" (lowercase) -> model_metabolite_id
        name_map = {}
        for m in self._model.metabolites:
            key = f"{m.name.lower()}[{m.compartment}]"
            name_map[key] = m.id
            # Also add without compartment brackets for fallback
            name_map[m.name.lower()] = m.id

        def translate_met_id(met_id: str) -> str:
            return name_map.get(met_id.lower(), name_map.get(met_id, met_id))

        translated = 0
        for task in tasks:
            # Translate inflow_dict
            new_inflow = {}
            for met_id, bounds in task.inflow_dict.items():
                new_id = translate_met_id(met_id)
                if new_id != met_id:
                    translated += 1
                new_inflow[new_id] = bounds
            task.inflow_dict = new_inflow

            # Translate outflow_dict
            new_outflow = {}
            for met_id, bounds in task.outflow_dict.items():
                new_id = translate_met_id(met_id)
                if new_id != met_id:
                    translated += 1
                new_outflow[new_id] = bounds
            task.outflow_dict = new_outflow

            # Translate reaction_dict stoichiometry (EQU column)
            # reaction_dict: {rxn_name: (stoich_dict, bounds)}
            # stoich_dict: {met_id: stoich_value}
            new_rxn_dict = {}
            for rxn_name, rxn_val in task.reaction_dict.items():
                if isinstance(rxn_val, (list, tuple)) and len(rxn_val) >= 1:
                    stoich = rxn_val[0]
                    bounds = rxn_val[1] if len(rxn_val) > 1 else (0, 1000)
                    new_stoich = {}
                    for met_id, coeff in stoich.items():
                        new_id = translate_met_id(met_id)
                        if new_id != met_id:
                            translated += 1
                        new_stoich[new_id] = coeff
                    new_rxn_dict[rxn_name] = (new_stoich, bounds)
                else:
                    new_rxn_dict[rxn_name] = rxn_val
            task.reaction_dict = new_rxn_dict

        if self.verbose and translated > 0:
            print(f"  Translated {translated} metabolite IDs to model format")
        return tasks

    def _get_medium_reaction_ids(self) -> set:
        """Return IDs of open exchange reactions (medium) in the model.

        A reaction is considered an exchange (medium) reaction if:
        - its ID starts with 'EX_', or
        - it involves a single metabolite and has lower_bound < 0 (uptake possible).

        Returns
        -------
        set of str
            Set of reaction IDs included in the medium.
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
        """Return names of tasks related to medium (exchange) conditions.

        A task is considered medium-related if:
        - any metabolite in task.inflow_dict matches a metabolite of a medium reaction, or
        - task.flux_constraints contains a medium exchange reaction ID.

        Parameters
        ----------
        tasks : list of Task
        medium_rxn_ids : set of str

        Returns
        -------
        set of str
            Set of medium-related task names.
        """
        medium_metabolites = set()
        for rxn in self._model.reactions:
            if rxn.id in medium_rxn_ids:
                for met in rxn.metabolites:
                    medium_metabolites.add(met.id)
                    # Also add compartment-stripped version (e.g. 'glc__D_e' → 'glc__D')
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
        """Build a bound_changes dict that sets inactive reactions to (0, 0).

        Each entry is an element of the bound_changes list passed to
        TaskEvaluator.batch_evaluate(). Active reactions are not included,
        preserving their original bounds.

        Parameters
        ----------
        result_dict : dict
            {reaction_id: bool} tINIT reconstruction result.

        Returns
        -------
        dict
            {reaction_id: (0.0, 0.0)} for all inactive reactions.
        """
        return {
            rxn.id: (0.0, 0.0)
            for rxn in self._model.reactions
            if not result_dict.get(rxn.id, False)
        }

    def _filter_compatible_tasks(self, tasks: list) -> list:
        """Remove tasks that are incompatible with the model.

        A task is incompatible if:
        - any metabolite ID in inflow_dict or outflow_dict is absent from the model, or
        - any reaction ID in mandatory_activity is absent from the model.

        Parameters
        ----------
        tasks : list of Task

        Returns
        -------
        list of Task
        """
        _ALLMETSIN_PLACEHOLDERS = {'ALLMETSIN[s]', 'ALLMETSIN[e]', 'ALLMETSIN'}
        model_met_ids = {m.id for m in self._model.metabolites}
        model_rxn_ids = {r.id for r in self._model.reactions}
        compatible = []
        skipped = 0
        for task in tasks:
            ok = True
            # Check inflow/outflow metabolites (skip ALLMETSIN placeholders)
            for met_id in list(task.inflow_dict.keys()) + list(task.outflow_dict.keys()):
                if met_id in _ALLMETSIN_PLACEHOLDERS:
                    continue
                if met_id not in model_met_ids:
                    ok = False
                    break
            # Check reaction_dict stoichiometry metabolites
            if ok:
                for rxn_val in task.reaction_dict.values():
                    if isinstance(rxn_val, (list, tuple)) and len(rxn_val) >= 1:
                        stoich = rxn_val[0]
                        for met_id in stoich.keys():
                            if met_id not in model_met_ids:
                                ok = False
                                break
                    if not ok:
                        break
            # Check mandatory_activity reactions
            if ok:
                for rxn_id in task.mandatory_activity:
                    if rxn_id not in model_rxn_ids:
                        ok = False
                        break
            if ok:
                compatible.append(task)
            else:
                skipped += 1
        if skipped > 0 and self.verbose:
            print(f"  Skipped {skipped} tasks incompatible with model "
                  f"({len(compatible)} compatible)")
        return compatible

    def _validate_tasks_for_sample(
        self,
        sample_name: str,
        result_dict: dict,
        tasks: list,
        medium_related_names: set,
    ):
        """Run metabolic task validation for a single sample.

        Parameters
        ----------
        sample_name : str
        result_dict : dict
            {reaction_id: bool} tINIT reconstruction result.
        tasks : list of Task
        medium_related_names : set of str
            Names of medium-related tasks.

        Returns
        -------
        TaskValidationResult
        """
        from troppo.tasks.core import TaskEvaluator, TaskValidationResult

        if self.verbose:
            print(f"    Validating {len(tasks)} tasks for {sample_name}...")

        tev = TaskEvaluator(model=self._model, tasks=tasks, solver=self.solver)
        bound_changes = [self._build_bound_changes(result_dict)]

        raw_results = tev.batch_evaluate(bound_changes=bound_changes, threads=1)

        # raw_results[0] = {task_name: (truth, expected, sol)}
        results = {
            task_name: truth
            for task_name, (truth, _, _) in raw_results[0].items()
        }

        return TaskValidationResult(
            sample_name=sample_name,
            results=results,
            medium_related=medium_related_names,
        )

    def _derive_task_essential_reactions(self, tasks: list) -> list:
        """Extract reaction IDs that must be included in reconstruction from task definitions.

        Reactions are selected by:
        - ``flux_constraints`` entries where lb > 0 (must carry flux), and
        - reactions listed in ``mandatory_activity``.

        Tasks with should_fail=True are excluded.

        Parameters
        ----------
        tasks : list of Task

        Returns
        -------
        list of str
            Reaction IDs to add to tINIT essential_reactions.
        """
        model_rxn_ids = {r.id for r in self._model.reactions}
        essential = set()

        for task in tasks:
            if task.should_fail:
                continue
            for rxn_id, bounds in task.flux_constraints.items():
                lb = bounds[0] if isinstance(bounds, (list, tuple)) else bounds
                if lb > 0 and rxn_id in model_rxn_ids:
                    essential.add(rxn_id)
            for rxn_id in task.mandatory_activity:
                if rxn_id in model_rxn_ids:
                    essential.add(rxn_id)

        if self.verbose:
            print(f"  Task-derived essential reactions: {len(essential)}")
        return list(essential)

    def derive_biomass_essential_reactions(
        self,
        biomass_rxn_id: str = None,
        min_fraction: float = 0.01,
        tol: float = 1e-6,
        processes: int = 1,
    ) -> tuple:
        """Identify biomass-essential reactions via FVA.

        Runs FVA on the reference model while constraining
        ``biomass_flux >= min_fraction * optimal``, and returns reactions
        that must carry non-zero flux to maintain that level.

        Pass the returned reaction IDs to ``configure(essential_reactions=...)``
        to force them into the tINIT reconstruction.

        Parameters
        ----------
        biomass_rxn_id : str, optional
            Biomass reaction ID. Auto-detected from model objective if None.
        min_fraction : float
            Minimum fraction of optimal biomass to maintain (default 0.01 = 1%).
            Lower values yield fewer essential reactions; higher values are stricter.
        tol : float
            Tolerance for considering flux non-zero (default 1e-6).
        processes : int
            Number of parallel FVA processes (default 1). Gurobi uses internal
            multi-threading, so 1 is usually sufficient.

        Returns
        -------
        tuple[list[str], pd.DataFrame]
            - ``essential_rxn_ids``: Reaction IDs to add to tINIT essential_reactions.
            - ``fva_df``: Full FVA result DataFrame.
              Columns: ``rxn_id``, ``min_flux``, ``max_flux``, ``is_essential``,
              ``rxn_name``, ``rxn_equation``.

        Examples
        --------
        >>> recon = TINITReconstructor(model_path="Human-GEM.xml", ...)
        >>> essential, fva_df = recon.derive_biomass_essential_reactions(min_fraction=0.01)
        >>> print(f"FVA essential reactions: {len(essential)}")
        >>> recon.configure(essential_reactions=essential)
        >>> results = recon.run(...)
        """
        import cobra.flux_analysis

        if biomass_rxn_id is None:
            biomass_rxn_id = self._detect_biomass_reaction()
        if biomass_rxn_id is None:
            raise ValueError(
                "Biomass reaction could not be detected. "
                "Pass biomass_rxn_id explicitly."
            )

        if self.verbose:
            print(f"  Detecting biomass-essential reactions via FVA ...")
            print(f"    Biomass reaction : {biomass_rxn_id}")

        with self._model as m:
            # 1. Compute optimal biomass
            m.objective = biomass_rxn_id
            opt_sol = m.optimize()

            if opt_sol.status != 'optimal' or opt_sol.objective_value < tol:
                raise RuntimeError(
                    f"Reference model cannot produce biomass "
                    f"(status: {opt_sol.status}, value: {opt_sol.objective_value})"
                )

            opt_biomass = opt_sol.objective_value
            min_biomass = opt_biomass * min_fraction

            if self.verbose:
                print(f"    Optimal biomass  : {opt_biomass:.6f}")
                print(f"    Min constraint   : {min_biomass:.8f} "
                      f"({min_fraction * 100:.1f}% of optimal)")

            # 2. Enforce minimum biomass flux
            m.reactions.get_by_id(biomass_rxn_id).lower_bound = min_biomass

            # 3. Run FVA
            #    fraction_of_optimum=0.0: the biomass lb already constrains feasibility,
            #    so we explore the full flux range without an additional objective constraint.
            if self.verbose:
                print(f"    Running FVA on {len(m.reactions)} reactions "
                      f"(processes={processes}) ...")

            fva = cobra.flux_analysis.flux_variability_analysis(
                m,
                reaction_list=None,
                fraction_of_optimum=0.0,
                processes=processes,
            )

        # 4. Extract essential reactions
        #    min_flux > tol  : must carry forward flux
        #    max_flux < -tol : must carry reverse flux
        essential = []
        rows = []
        for rxn_id, row in fva.iterrows():
            min_f = float(row['minimum'])
            max_f = float(row['maximum'])
            is_ess = (min_f > tol) or (max_f < -tol)
            if is_ess:
                essential.append(rxn_id)
            try:
                rxn      = self._model.reactions.get_by_id(rxn_id)
                rxn_name = rxn.name
                rxn_eq   = rxn.reaction
            except Exception:
                rxn_name = ''
                rxn_eq   = ''
            rows.append({
                'rxn_id':       rxn_id,
                'min_flux':     round(min_f, 8),
                'max_flux':     round(max_f, 8),
                'is_essential': is_ess,
                'rxn_name':     rxn_name,
                'rxn_equation': rxn_eq,
            })

        fva_df = pd.DataFrame(rows)

        if self.verbose:
            print(f"    FVA essential reactions: {len(essential)} "
                  f"/ {len(fva_df)} total")

        return essential, fva_df

    def _detect_biomass_reaction(self) -> str:
        """Detect the biomass reaction ID from the model objective or reaction names."""
        try:
            from cobra.util.solver import linear_reaction_coefficients
            coeffs = linear_reaction_coefficients(self._model)
            if coeffs:
                return max(coeffs, key=lambda r: abs(coeffs[r])).id
        except Exception:
            pass
        for rxn in self._model.reactions:
            if 'biomass' in rxn.id.lower():
                return rxn.id
        return None

    def _gap_fill_for_failed_tasks(
        self,
        result_dict: dict,
        tasks: list,
        failed_task_names: set,
    ) -> tuple:
        """Apply targeted gap-filling for failed metabolic tasks.

        Strategy (applied in order):
        1. Activate reactions from ``flux_constraints`` where lb > 0 but currently inactive.
        2. Activate reactions listed in ``mandatory_activity`` that are inactive.
        3. Recursively activate producer reactions for ``outflow_dict`` metabolites.

        Parameters
        ----------
        result_dict : dict
            {reaction_id: bool} current reconstruction result.
        tasks : list of Task
        failed_task_names : set of str
            Names of tasks to gap-fill.

        Returns
        -------
        tuple[dict, dict]
            (updated_result_dict, {task_name: set_of_added_rxn_ids})
        """
        model_rxn_ids = {r.id for r in self._model.reactions}
        met_map = {m.id: m for m in self._model.metabolites}
        met_base_map = {}  # base_id (compartment stripped) → cobra Metabolite
        for m in self._model.metabolites:
            base = m.id.rsplit('_', 1)[0] if '_' in m.id else m.id
            if base not in met_base_map:
                met_base_map[base] = m

        updated = dict(result_dict)
        gap_fill_log = {}  # task_name → set of added rxn IDs

        def _activate_producer(met_obj, depth, max_depth, visited_mets, added_local):
            """Recursively activate producer reactions for a metabolite.

            Returns immediately if an active producer already exists.
            Otherwise selects the first candidate and recursively activates its substrates.
            """
            if depth > max_depth or met_obj.id in visited_mets:
                return
            visited_mets.add(met_obj.id)

            # Return early if an active producer already exists
            for rxn in met_obj.reactions:
                coeff = rxn.metabolites[met_obj]
                can_produce = (coeff > 0 and rxn.upper_bound > 0) or (
                    coeff < 0 and rxn.lower_bound < 0
                )
                if can_produce and updated.get(rxn.id, False):
                    return

            # Select the first feasible producer and activate it
            for rxn in met_obj.reactions:
                coeff = rxn.metabolites[met_obj]
                can_produce = (coeff > 0 and rxn.upper_bound > 0) or (
                    coeff < 0 and rxn.lower_bound < 0
                )
                if not can_produce:
                    continue
                # Activate this reaction
                if not updated.get(rxn.id, False):
                    updated[rxn.id] = True
                    added_local.add(rxn.id)
                # Recursively activate its substrates
                for sub_met, sub_coeff in rxn.metabolites.items():
                    produces_sub = (sub_coeff > 0 and rxn.lower_bound < 0) or (
                        sub_coeff < 0 and rxn.upper_bound > 0
                    )
                    is_substrate = (coeff > 0 and sub_coeff < 0) or (
                        coeff < 0 and sub_coeff > 0
                    )
                    if is_substrate:
                        _activate_producer(
                            sub_met, depth + 1, max_depth, visited_mets, added_local
                        )
                break  # process only one producer reaction

        for task in tasks:
            if task.name not in failed_task_names or task.should_fail:
                continue

            added = set()

            # 1. flux_constraints lb > 0 → must carry flux, activate
            for rxn_id, bounds in task.flux_constraints.items():
                lb = bounds[0] if isinstance(bounds, (list, tuple)) else bounds
                if lb > 0 and rxn_id in model_rxn_ids and not updated.get(rxn_id, False):
                    updated[rxn_id] = True
                    added.add(rxn_id)

            # 2. Activate mandatory_activity reactions
            for rxn_id in task.mandatory_activity:
                if rxn_id in model_rxn_ids and not updated.get(rxn_id, False):
                    updated[rxn_id] = True
                    added.add(rxn_id)

            # 3. Recursively activate producers for outflow_dict metabolites (depth limit: 10)
            for met_id in task.outflow_dict:
                met_obj = met_map.get(met_id) or met_base_map.get(
                    met_id.rsplit('_', 1)[0] if '_' in met_id else met_id
                )
                if met_obj is None:
                    continue
                _activate_producer(met_obj, depth=0, max_depth=10,
                                   visited_mets=set(), added_local=added)

            if added:
                gap_fill_log[task.name] = added

        return updated, gap_fill_log

    def _validate_tasks_streaming(
        self,
        sample_name: str,
        result_dict: dict,
        tasks: list,
        medium_related_names: set,
        gap_filled_tasks: set = None,
    ):
        """Evaluate tasks and print each result in real time.

        Parameters
        ----------
        sample_name : str
        result_dict : dict
        tasks : list of Task
        medium_related_names : set of str
        gap_filled_tasks : set of str, optional
            Task names recovered by gap-filling (shown as a tag in output).

        Returns
        -------
        TaskValidationResult
        """
        from troppo.tasks.core import TaskEvaluator, TaskValidationResult

        gap_filled_tasks = gap_filled_tasks or set()
        bound_changes = [self._build_bound_changes(result_dict)]
        tev = TaskEvaluator(model=self._model, tasks=tasks, solver=self.solver)

        # Full batch evaluation (preserves batch optimization)
        raw_results = tev.batch_evaluate(bound_changes=bound_changes, threads=1)

        results = {}
        sep = "─" * 56
        print(f"\n    {sep}")
        print(f"    Task Validation (real-time): {sample_name}")
        print(f"    {sep}")

        passed_count = 0
        for task in tasks:
            truth = raw_results[0].get(task.name, (False, {}, None))[0]
            results[task.name] = truth
            if truth:
                passed_count += 1

            status = "✅" if truth else "❌"
            med_tag = " [MED]" if task.name in medium_related_names else ""
            gf_tag = " [GAP-FILLED]" if task.name in gap_filled_tasks else ""
            print(f"    {status} {task.name}{med_tag}{gf_tag}")

        total = len(tasks)
        print(f"    {sep}")
        print(f"    Passed: {passed_count}/{total} ({passed_count/max(total,1)*100:.1f}%)")
        med_tasks = [t.name for t in tasks if t.name in medium_related_names]
        if med_tasks:
            med_passed = sum(1 for t in med_tasks if results.get(t, False))
            print(f"    Medium-related: {med_passed}/{len(med_tasks)}")
        print(f"    {sep}")

        return TaskValidationResult(
            sample_name=sample_name,
            results=results,
            medium_related=medium_related_names,
        )

    def _print_task_report(self, tvr, gap_filled_tasks: set = None) -> None:
        """Print a task validation summary report to the console.

        Parameters
        ----------
        tvr : TaskValidationResult
        gap_filled_tasks : set of str, optional
            Task names recovered by gap-filling.
        """
        gap_filled_tasks = gap_filled_tasks or set()
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
        if gap_filled_tasks:
            print(f"    Gap-filled tasks: {len(gap_filled_tasks)}")
        if tvr.failed_tasks:
            print(f"    Failed tasks:")
            for task_name in sorted(tvr.failed_tasks):
                tags = []
                if task_name in tvr.medium_related:
                    tags.append("MEDIUM")
                tag_str = f" [{', '.join(tags)}]" if tags else ""
                print(f"      - {task_name}{tag_str}")
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
        task_guided: bool = False,
        stream_tasks: bool = True,
    ) -> 'ReconstructionResults':
        """Run tINIT reconstruction for all (or selected) samples.

        Parameters
        ----------
        samples : list, optional
            Sample names to reconstruct. None = all samples.
        gap_fill : bool
            Whether to apply gap-filling after reconstruction.
        validate_tasks : bool
            Whether to validate metabolic tasks after reconstruction.
        task_file : str, optional
            Path to metabolic task Excel/JSON file.
        task_guided : bool
            When True, task constraints inform the tINIT optimization:

            1. **Pre-tINIT**: reactions from ``flux_constraints`` (lb > 0) and
               ``mandatory_activity`` are added to ``essential_reactions``,
               forcing them into the reconstructed model.
            2. **Post-tINIT gap-fill**: for each failed task, the minimum set
               of reactions needed (flux_constraints, mandatory_activity,
               outflow producers) is added and tasks are re-validated.
            3. **Real-time display**: each task result is printed immediately.

            Implies ``validate_tasks=True``; ``task_file`` is required.
        stream_tasks : bool
            When True (default), print each task result line-by-line in
            real time. When False, print only the summary at the end.

        Returns
        -------
        ReconstructionResults
            Results container with reaction matrix and analysis methods.
        """
        from troppo.results import ReconstructionResults

        if self._expression_df is None:
            raise ValueError("No expression data loaded. Provide expression_path or expression_data.")

        # task_guided implies validate_tasks
        if task_guided:
            validate_tasks = True

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
                "task_file must be provided when validate_tasks=True or task_guided=True. "
                "Example: run(validate_tasks=True, task_file='metabolicTasks_Essential.xlsx')"
            )

        if self.verbose:
            mode_tag = " [task-guided]" if task_guided else ""
            print(f"\nRunning tINIT reconstruction for {len(selected)} samples "
                  f"(solver={self.solver}){mode_tag}...")

        # Prepare tasks if validation (or guidance) is requested
        tasks = None
        medium_related_names = set()
        task_essential_rxn_ids = []
        if validate_tasks:
            tasks = self._load_tasks(task_file)
            tasks = self._filter_compatible_tasks(tasks)
            medium_rxn_ids = self._get_medium_reaction_ids()
            medium_related_names = self._identify_medium_related_tasks(tasks, medium_rxn_ids)
            if not medium_related_names and self.verbose:
                warnings.warn(
                    "No medium-related tasks identified. "
                    "Check that task inflow metabolites match model exchange reaction metabolites."
                )

            # Pre-tINIT: derive essential reactions from task definitions
            if task_guided:
                task_essential_rxn_ids = self._derive_task_essential_reactions(tasks)
                # Temporarily extend essential_reactions for this run
                _orig_essential = self._essential_reactions
                self._essential_reactions = list(self._essential_reactions or []) + [
                    r for r in task_essential_rxn_ids
                    if not self._essential_reactions or r not in self._essential_reactions
                ]
                if self.verbose and task_essential_rxn_ids:
                    print(f"  Added {len(task_essential_rxn_ids)} task-essential reactions to reconstruction")

        # Run per sample
        results_dict = {}
        task_results = {}
        try:
            for i, sample in enumerate(selected):
                expression_dict = self._expression_df[sample].dropna().to_dict()
                # Ensure keys are strings
                expression_dict = {str(k): float(v) for k, v in expression_dict.items()}
                result = self._run_single_sample(sample, expression_dict)

                if self.verbose:
                    n_active = sum(1 for v in result.values() if v)
                    print(f"    [{i+1}/{len(selected)}] {sample}: {n_active} active reactions")

                # --- Task validation & gap-fill ---
                gap_filled_tasks: set = set()
                if validate_tasks and tasks:
                    # Initial validation
                    if stream_tasks and task_guided:
                        tvr = self._validate_tasks_streaming(
                            sample, result, tasks, medium_related_names
                        )
                    else:
                        if self.verbose:
                            print(f"    Validating {len(tasks)} tasks for {sample}...")
                        tvr = self._validate_tasks_for_sample(
                            sample, result, tasks, medium_related_names
                        )

                    # task_guided: gap-fill for failed tasks then re-validate
                    if task_guided and tvr.failed_tasks:
                        if self.verbose:
                            print(f"\n    Gap-filling for {len(tvr.failed_tasks)} failed tasks...")
                        result, gap_fill_log = self._gap_fill_for_failed_tasks(
                            result, tasks, set(tvr.failed_tasks)
                        )
                        gap_filled_tasks = set(gap_fill_log.keys())
                        total_added = sum(len(v) for v in gap_fill_log.values())
                        if self.verbose and gap_fill_log:
                            for task_name, rxns in gap_fill_log.items():
                                rxn_preview = ', '.join(sorted(rxns)[:4])
                                ellipsis = '...' if len(rxns) > 4 else ''
                                print(f"      [{task_name}] +{len(rxns)} rxns: {rxn_preview}{ellipsis}")

                        # Re-validate after gap-fill
                        if self.verbose:
                            print(f"    Re-validating after gap-fill (+{total_added} reactions)...")
                        if stream_tasks:
                            tvr = self._validate_tasks_streaming(
                                sample, result, tasks, medium_related_names,
                                gap_filled_tasks=gap_filled_tasks,
                            )
                        else:
                            tvr = self._validate_tasks_for_sample(
                                sample, result, tasks, medium_related_names
                            )
                            if self.verbose:
                                self._print_task_report(tvr, gap_filled_tasks)
                    elif self.verbose and not stream_tasks:
                        self._print_task_report(tvr)
                    elif self.verbose and not task_guided:
                        self._print_task_report(tvr)

                    task_results[sample] = tvr

                results_dict[sample] = result
        finally:
            # Restore essential_reactions if task_guided modified them
            if task_guided and validate_tasks:
                self._essential_reactions = _orig_essential

        # Build metadata
        metadata = {
            'solver': self.solver,
            'production_weight': self._production_weight,
            'allow_excretion': self._allow_excretion,
            'no_reverse_loops': self._no_reverse_loops,
            'score_integration': self._score_integration,
            'model_path': self.model_path,
            'expression_path': self.expression_path,
            'task_guided': task_guided,
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
            print(f"  Samples        : {results.n_samples}")
            print(f"  Active rxns    : {summary['n_active_reactions'].mean():.0f} (mean)")
            print(f"  Common rxns    : {len(results.get_common_reactions())}")
            if task_results:
                ts = results.task_summary()
                print(f"  Task pass rate : {ts['pass_rate'].mean():.1f}% (mean)")
                if 'medium_pass_rate' in ts.columns:
                    print(f"  Medium tasks   : {ts['medium_pass_rate'].mean():.1f}% (mean)")

        return results

    def __repr__(self):
        expr_info = ""
        if self._expression_df is not None:
            expr_info = f", expression=({self._expression_df.shape[0]} genes, {self._expression_df.shape[1]} samples)"
        return (f"TINITReconstructor(model={self.model_path}, solver={self.solver}"
                f"{expr_info})")
