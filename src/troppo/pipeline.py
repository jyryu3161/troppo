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

        if self.verbose:
            print(f"\nRunning tINIT reconstruction for {len(selected)} samples "
                  f"(solver={self.solver})...")

        # Run per sample
        results_dict = {}
        for i, sample in enumerate(selected):
            expression_dict = self._expression_df[sample].dropna().to_dict()
            # Ensure keys are strings
            expression_dict = {str(k): float(v) for k, v in expression_dict.items()}
            result = self._run_single_sample(sample, expression_dict)
            results_dict[sample] = result

            if self.verbose:
                n_active = sum(1 for v in result.values() if v)
                print(f"    [{i+1}/{len(selected)}] {sample}: {n_active} active reactions")

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
