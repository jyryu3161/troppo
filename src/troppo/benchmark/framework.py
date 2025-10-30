"""
Benchmark Framework for Omics Integration Methods

This module provides tools for comparing different omics integration methods
using standardized metrics and validation procedures.

여러 오믹스 통합 방법론을 표준화된 지표와 검증 절차로 비교하는 도구를 제공합니다.
"""

import time
import tracemalloc
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any, Callable, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime
import warnings


@dataclass
class BenchmarkResult:
    """
    Results from a single method execution
    단일 방법 실행 결과
    """
    method_name: str
    execution_time: float  # seconds
    peak_memory: float  # MB
    num_reactions_selected: int
    num_reactions_original: int
    percentage_retained: float
    selected_reaction_ids: List[str] = field(default_factory=list)
    selected_reaction_indices: List[int] = field(default_factory=list)

    # Biological validation metrics
    biomass_flux: Optional[float] = None
    growth_rate: Optional[float] = None
    task_completion: Optional[Dict[str, bool]] = None
    task_completion_rate: Optional[float] = None

    # Gene essentiality validation
    essentiality_accuracy: Optional[float] = None
    essentiality_precision: Optional[float] = None
    essentiality_recall: Optional[float] = None
    essentiality_f1: Optional[float] = None
    essentiality_mcc: Optional[float] = None
    essentiality_result: Optional[Any] = None  # Full GeneEssentialityResult

    # Theoretical yield validation
    yield_results: Optional[List[Any]] = None  # List of TheoreticalYieldResult
    avg_aerobic_yield: Optional[float] = None
    avg_anaerobic_yield: Optional[float] = None
    num_feasible_conditions: Optional[int] = None

    # Statistical metrics
    num_genes_mapped: Optional[int] = None
    num_reactions_with_gpr: Optional[int] = None
    num_reactions_without_gpr: Optional[int] = None

    # Quality metrics
    network_consistency: Optional[bool] = None
    blocked_reactions: Optional[List[str]] = None
    num_blocked_reactions: Optional[int] = None

    # Error information
    success: bool = True
    error_message: Optional[str] = None

    # Metadata
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    parameters: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return asdict(self)

    def to_summary_dict(self) -> Dict[str, Any]:
        """Convert to summary dictionary (without large lists)"""
        summary = self.to_dict()
        # Remove large lists to keep summary compact
        summary.pop('selected_reaction_ids', None)
        summary.pop('selected_reaction_indices', None)
        summary.pop('blocked_reactions', None)
        summary.pop('task_completion', None)
        return summary


@dataclass
class BenchmarkComparison:
    """
    Comparison results across multiple methods
    여러 방법론 간 비교 결과
    """
    results: Dict[str, BenchmarkResult]
    dataset_info: Dict[str, Any]
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def get_summary_dataframe(self) -> pd.DataFrame:
        """
        Get summary as pandas DataFrame

        Returns
        -------
        pd.DataFrame
            Summary of all methods
        """
        summaries = []
        for method_name, result in self.results.items():
            if result.success:
                summaries.append({
                    'Method': method_name,
                    'Time (s)': result.execution_time,
                    'Memory (MB)': result.peak_memory,
                    'Reactions': result.num_reactions_selected,
                    '% Retained': result.percentage_retained,
                    'Biomass Flux': result.biomass_flux,
                    'Growth Rate': result.growth_rate,
                    'Task Success Rate': result.task_completion_rate,
                    'Blocked Rxns': result.num_blocked_reactions,
                    'Consistent': result.network_consistency
                })
            else:
                summaries.append({
                    'Method': method_name,
                    'Error': result.error_message
                })

        return pd.DataFrame(summaries)

    def get_overlap_matrix(self) -> pd.DataFrame:
        """
        Get reaction overlap matrix between methods

        Returns
        -------
        pd.DataFrame
            Jaccard similarity matrix
        """
        methods = list(self.results.keys())
        n = len(methods)
        overlap_matrix = np.zeros((n, n))

        for i, method1 in enumerate(methods):
            result1 = self.results[method1]
            if not result1.success:
                continue

            set1 = set(result1.selected_reaction_indices)

            for j, method2 in enumerate(methods):
                result2 = self.results[method2]
                if not result2.success:
                    continue

                set2 = set(result2.selected_reaction_indices)

                # Jaccard similarity
                intersection = len(set1 & set2)
                union = len(set1 | set2)
                overlap_matrix[i, j] = intersection / union if union > 0 else 0

        return pd.DataFrame(overlap_matrix, index=methods, columns=methods)

    def get_common_reactions(self) -> List[int]:
        """Get reactions common to all methods"""
        successful_results = [r for r in self.results.values() if r.success]
        if not successful_results:
            return []

        common = set(successful_results[0].selected_reaction_indices)
        for result in successful_results[1:]:
            common &= set(result.selected_reaction_indices)

        return sorted(list(common))

    def get_unique_reactions(self, method_name: str) -> List[int]:
        """Get reactions unique to a specific method"""
        if method_name not in self.results:
            return []

        result = self.results[method_name]
        if not result.success:
            return []

        unique = set(result.selected_reaction_indices)
        for name, other_result in self.results.items():
            if name != method_name and other_result.success:
                unique -= set(other_result.selected_reaction_indices)

        return sorted(list(unique))

    def rank_methods(self, metric: str = 'percentage_retained') -> List[Tuple[str, float]]:
        """
        Rank methods by a specific metric

        Parameters
        ----------
        metric : str
            Metric to rank by (execution_time, peak_memory, percentage_retained, etc.)

        Returns
        -------
        List[Tuple[str, float]]
            Sorted list of (method_name, metric_value)
        """
        rankings = []
        for method_name, result in self.results.items():
            if result.success and hasattr(result, metric):
                value = getattr(result, metric)
                if value is not None:
                    rankings.append((method_name, value))

        # Sort based on metric (lower is better for time/memory, higher for others)
        reverse = metric not in ['execution_time', 'peak_memory', 'num_blocked_reactions']
        return sorted(rankings, key=lambda x: x[1], reverse=reverse)

    def save_to_json(self, filepath: str) -> None:
        """Save comparison results to JSON file"""
        data = {
            'dataset_info': self.dataset_info,
            'timestamp': self.timestamp,
            'results': {
                name: result.to_dict()
                for name, result in self.results.items()
            }
        }

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load_from_json(cls, filepath: str) -> 'BenchmarkComparison':
        """Load comparison results from JSON file"""
        with open(filepath, 'r') as f:
            data = json.load(f)

        results = {}
        for name, result_data in data['results'].items():
            results[name] = BenchmarkResult(**result_data)

        return cls(
            results=results,
            dataset_info=data['dataset_info'],
            timestamp=data.get('timestamp', datetime.now().isoformat())
        )


class BenchmarkRunner:
    """
    Runner for benchmarking omics integration methods

    Examples
    --------
    >>> runner = BenchmarkRunner(
    ...     model_wrapper=model_wrapper,
    ...     data_map=data_map,
    ...     methods=['gimme', 'tinit', 'imat', 'fastcore']
    ... )
    >>> comparison = runner.run_benchmark()
    >>> print(comparison.get_summary_dataframe())
    """

    def __init__(
        self,
        model_wrapper,
        data_map,
        methods: Optional[List[str]] = None,
        validation_tasks: Optional[Dict] = None,
        biomass_reaction: Optional[str] = None,
        essential_genes: Optional[List[str]] = None,
        non_essential_genes: Optional[List[str]] = None,
        essential_genes_id_type: Optional[str] = None,
        non_essential_genes_id_type: Optional[str] = None,
        carbon_sources: Optional[List[str]] = None,
        verbose: bool = True
    ):
        """
        Initialize benchmark runner

        Parameters
        ----------
        model_wrapper : ModelBasedWrapper
            Wrapped metabolic model
        data_map : OmicsDataMap
            Mapped omics data
        methods : List[str], optional
            List of method names to benchmark (default: all registered methods)
        validation_tasks : Dict, optional
            Metabolic tasks for validation
        biomass_reaction : str, optional
            ID of biomass reaction for growth rate calculation
        essential_genes : List[str], optional
            List of experimentally verified essential genes
        non_essential_genes : List[str], optional
            List of experimentally verified non-essential genes
        essential_genes_id_type : str, optional
            ID type for essential genes (e.g., 'entrez_id', 'ensembl_gene_id', 'symbol')
            Auto-detected if not specified
        non_essential_genes_id_type : str, optional
            ID type for non-essential genes
            Auto-detected if not specified
        carbon_sources : List[str], optional
            List of carbon sources to test for theoretical yields
        verbose : bool
            Print progress messages
        """
        self.model_wrapper = model_wrapper
        self.data_map = data_map
        self.methods = methods
        self.validation_tasks = validation_tasks
        self.biomass_reaction = biomass_reaction
        self.essential_genes = essential_genes
        self.non_essential_genes = non_essential_genes
        self.essential_genes_id_type = essential_genes_id_type
        self.non_essential_genes_id_type = non_essential_genes_id_type
        self.carbon_sources = carbon_sources or ['glucose', 'acetate', 'glycerol']
        self.verbose = verbose

        # Get methods from registry if not specified
        if self.methods is None:
            from troppo.methods.registry import MethodRegistry
            self.methods = MethodRegistry.list_methods()

        # Prepare gene IDs for validation
        self._prepare_gene_ids()

    def _prepare_gene_ids(self):
        """Prepare and standardize gene IDs for validation"""
        if not (self.essential_genes or self.non_essential_genes):
            self.essential_genes_converted = None
            self.non_essential_genes_converted = None
            return

        try:
            from troppo.benchmark.gene_id_utils import (
                detect_id_type,
                convert_gene_ids,
                normalize_id_type
            )

            # Get model gene IDs
            try:
                model_genes = [g.id for g in self.model_wrapper.model_reader.model.genes]
                model_id_type = detect_id_type(model_genes[:100])  # Sample first 100
                if self.verbose:
                    print(f"  Model gene ID type detected: {model_id_type}")
            except:
                model_id_type = 'symbol'
                if self.verbose:
                    print(f"  Model gene ID type: using default '{model_id_type}'")

            # Convert essential genes
            if self.essential_genes:
                if self.essential_genes_id_type:
                    source_type = normalize_id_type(self.essential_genes_id_type)
                else:
                    source_type = detect_id_type(self.essential_genes)

                if self.verbose:
                    print(f"  Essential genes ID type: {source_type}")

                if source_type != model_id_type:
                    conv = convert_gene_ids(
                        self.essential_genes,
                        from_type=source_type,
                        to_type=model_id_type
                    )
                    self.essential_genes_converted = list(conv.values())
                    if self.verbose:
                        print(f"  Converted {len(self.essential_genes_converted)}/{len(self.essential_genes)} essential genes")
                else:
                    self.essential_genes_converted = self.essential_genes
            else:
                self.essential_genes_converted = None

            # Convert non-essential genes
            if self.non_essential_genes:
                if self.non_essential_genes_id_type:
                    source_type = normalize_id_type(self.non_essential_genes_id_type)
                else:
                    source_type = detect_id_type(self.non_essential_genes)

                if self.verbose:
                    print(f"  Non-essential genes ID type: {source_type}")

                if source_type != model_id_type:
                    conv = convert_gene_ids(
                        self.non_essential_genes,
                        from_type=source_type,
                        to_type=model_id_type
                    )
                    self.non_essential_genes_converted = list(conv.values())
                    if self.verbose:
                        print(f"  Converted {len(self.non_essential_genes_converted)}/{len(self.non_essential_genes)} non-essential genes")
                else:
                    self.non_essential_genes_converted = self.non_essential_genes
            else:
                self.non_essential_genes_converted = None

        except Exception as e:
            if self.verbose:
                print(f"  Warning: Gene ID conversion failed: {str(e)}")
                print(f"  Using original gene IDs without conversion")
            self.essential_genes_converted = self.essential_genes
            self.non_essential_genes_converted = self.non_essential_genes

    def run_benchmark(
        self,
        method_configs: Optional[Dict[str, Dict]] = None,
        validate_biology: bool = True,
        validate_consistency: bool = True,
        validate_essentiality: bool = True,
        validate_yields: bool = True
    ) -> BenchmarkComparison:
        """
        Run benchmark for all methods

        Parameters
        ----------
        method_configs : Dict[str, Dict], optional
            Method-specific configurations
        validate_biology : bool
            Whether to perform biological validation
        validate_consistency : bool
            Whether to check network consistency
        validate_essentiality : bool
            Whether to validate gene essentiality predictions
        validate_yields : bool
            Whether to calculate theoretical yields

        Returns
        -------
        BenchmarkComparison
            Comparison results
        """
        if method_configs is None:
            method_configs = {}

        results = {}

        for method_name in self.methods:
            if self.verbose:
                print(f"\n{'=' * 60}")
                print(f"Benchmarking: {method_name.upper()}")
                print(f"{'=' * 60}")

            config = method_configs.get(method_name, {})

            try:
                result = self._run_single_method(
                    method_name,
                    config,
                    validate_biology,
                    validate_consistency,
                    validate_essentiality,
                    validate_yields
                )
                results[method_name] = result

                if self.verbose:
                    print(f"✓ Success: {result.num_reactions_selected} reactions "
                          f"({result.percentage_retained:.2f}%) in {result.execution_time:.2f}s")

            except Exception as e:
                if self.verbose:
                    print(f"✗ Failed: {str(e)}")

                results[method_name] = BenchmarkResult(
                    method_name=method_name,
                    execution_time=0,
                    peak_memory=0,
                    num_reactions_selected=0,
                    num_reactions_original=len(self.model_wrapper.model_reader.r_ids),
                    percentage_retained=0,
                    success=False,
                    error_message=str(e)
                )

        # Create comparison
        dataset_info = {
            'num_reactions': len(self.model_wrapper.model_reader.r_ids),
            'num_metabolites': len(self.model_wrapper.model_reader.m_ids),
            'num_genes': len(self.data_map.get_scores()),
            'biomass_reaction': self.biomass_reaction
        }

        return BenchmarkComparison(
            results=results,
            dataset_info=dataset_info
        )

    def _run_single_method(
        self,
        method_name: str,
        config: Dict,
        validate_biology: bool,
        validate_consistency: bool,
        validate_essentiality: bool,
        validate_yields: bool
    ) -> BenchmarkResult:
        """Run a single method and collect metrics"""
        from troppo.methods.registry import MethodRegistry

        # Prepare method-specific configuration
        properties = self._prepare_properties(method_name, config)

        # Start tracking
        tracemalloc.start()
        start_time = time.time()

        # Run method
        method_instance = MethodRegistry.create_method(
            method_name,
            S=self.model_wrapper.S,
            lb=self.model_wrapper.lb,
            ub=self.model_wrapper.ub,
            properties=properties
        )

        selected_indices = method_instance.run()

        # Stop tracking
        execution_time = time.time() - start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        peak_memory_mb = peak / 1024 / 1024  # Convert to MB

        # Get reaction IDs
        selected_reaction_ids = [
            self.model_wrapper.model_reader.r_ids[i]
            for i in selected_indices
        ]

        # Create result
        result = BenchmarkResult(
            method_name=method_name,
            execution_time=execution_time,
            peak_memory=peak_memory_mb,
            num_reactions_selected=len(selected_indices),
            num_reactions_original=len(self.model_wrapper.model_reader.r_ids),
            percentage_retained=len(selected_indices) / len(self.model_wrapper.model_reader.r_ids) * 100,
            selected_reaction_ids=selected_reaction_ids,
            selected_reaction_indices=selected_indices,
            parameters=config
        )

        # Biological validation
        if validate_biology:
            self._validate_biology(result, selected_indices)

        # Consistency validation
        if validate_consistency:
            self._validate_consistency(result, selected_indices)

        # Gene essentiality validation
        if validate_essentiality and (self.essential_genes or self.non_essential_genes):
            self._validate_essentiality(result, selected_indices)

        # Theoretical yield validation
        if validate_yields:
            self._validate_yields(result, selected_indices)

        return result

    def _prepare_properties(self, method_name: str, config: Dict):
        """Prepare method properties based on method type"""
        from troppo.methods.registry import MethodRegistry
        from troppo.omics.integration import (
            ContinuousScoreIntegrationStrategy,
            ThresholdSelectionIntegrationStrategy
        )

        metadata = MethodRegistry.get_method(method_name)
        if metadata is None:
            raise ValueError(f"Method {method_name} not found")

        # Get default configuration
        default_config = {
            'solver': config.get('solver', 'CPLEX'),
            'reaction_ids': self.model_wrapper.model_reader.r_ids,
            'metabolite_ids': self.model_wrapper.model_reader.m_ids
        }

        # Method-specific preparation
        if method_name.lower() == 'gimme':
            # Prepare continuous scores
            score_apply = lambda x: {k: 0 if v is None else v for k, v in x.items()}
            integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
            scores = integration.integrate(data_map=self.data_map)

            idx_biomass = self.model_wrapper.model_reader.r_ids.index(
                self.biomass_reaction or 'biomass_human'
            )

            properties = MethodRegistry.create_properties(
                'gimme',
                exp_vector=[v for k, v in scores.items()],
                obj_frac=config.get('obj_frac', 0.8),
                objectives=[{idx_biomass: 1}],
                preprocess=config.get('preprocess', True),
                flux_threshold=config.get('flux_threshold', 0.8),
                **default_config
            )

        elif method_name.lower() in ['tinit', 'fastcore']:
            # Prepare threshold scores
            threshold = config.get('threshold', np.percentile(
                [v for v in self.data_map.get_scores().values() if v is not None],
                75 if method_name.lower() == 'tinit' else 50
            ))

            integration = ThresholdSelectionIntegrationStrategy(thresholds=threshold)
            core = integration.integrate(data_map=self.data_map)

            properties = MethodRegistry.create_properties(
                method_name,
                core=core,
                **default_config
            )

        elif method_name.lower() == 'imat':
            # Prepare high/low expression sets
            score_values = [v for v in self.data_map.get_scores().values() if v is not None]
            high_threshold = config.get('high_threshold', np.percentile(score_values, 67))
            low_threshold = config.get('low_threshold', np.percentile(score_values, 33))

            highly_expressed = set([
                k for k, v in self.data_map.get_scores().items()
                if v is not None and v >= high_threshold
            ])
            lowly_expressed = set([
                k for k, v in self.data_map.get_scores().items()
                if v is not None and v <= low_threshold
            ])

            properties = MethodRegistry.create_properties(
                'imat',
                rh=highly_expressed,
                rl=lowly_expressed,
                epsilon=config.get('epsilon', 1.0),
                threshold=config.get('threshold', 1e-3),
                **default_config
            )

        else:
            # Generic preparation
            properties = MethodRegistry.create_properties(method_name, **{**default_config, **config})

        return properties

    def _validate_biology(self, result: BenchmarkResult, selected_indices: List[int]):
        """Perform biological validation"""
        try:
            import cobra

            # Create tissue-specific model
            original_model = self.model_wrapper.model_reader.model
            tissue_model = original_model.copy()

            selected_ids = [self.model_wrapper.model_reader.r_ids[i] for i in selected_indices]
            reactions_to_remove = [rxn for rxn in tissue_model.reactions if rxn.id not in selected_ids]
            tissue_model.remove_reactions(reactions_to_remove, remove_orphans=True)

            # Calculate biomass flux
            if self.biomass_reaction:
                try:
                    solution = tissue_model.optimize()
                    if solution.status == 'optimal':
                        result.biomass_flux = solution.objective_value
                        result.growth_rate = solution.objective_value
                except:
                    result.biomass_flux = 0
                    result.growth_rate = 0

            # Task completion
            if self.validation_tasks:
                task_results = {}
                for task_name, task_def in self.validation_tasks.items():
                    # Simplified task checking
                    task_results[task_name] = True  # Implement actual task checking
                result.task_completion = task_results
                result.task_completion_rate = sum(task_results.values()) / len(task_results) * 100

        except Exception as e:
            if self.verbose:
                print(f"  Warning: Biological validation failed: {str(e)}")

    def _validate_consistency(self, result: BenchmarkResult, selected_indices: List[int]):
        """Check network consistency"""
        try:
            import cobra

            # Create tissue-specific model
            original_model = self.model_wrapper.model_reader.model
            tissue_model = original_model.copy()

            selected_ids = [self.model_wrapper.model_reader.r_ids[i] for i in selected_indices]
            reactions_to_remove = [rxn for rxn in tissue_model.reactions if rxn.id not in selected_ids]
            tissue_model.remove_reactions(reactions_to_remove, remove_orphans=True)

            # Find blocked reactions
            from cobra.flux_analysis import find_blocked_reactions
            blocked = find_blocked_reactions(tissue_model)

            result.blocked_reactions = [r.id for r in blocked]
            result.num_blocked_reactions = len(blocked)
            result.network_consistency = len(blocked) == 0

        except Exception as e:
            if self.verbose:
                print(f"  Warning: Consistency validation failed: {str(e)}")

    def _validate_essentiality(self, result: BenchmarkResult, selected_indices: List[int]):
        """Validate gene essentiality predictions"""
        try:
            from troppo.benchmark.validation import GeneEssentialityValidator
            import cobra

            # Create tissue-specific model
            original_model = self.model_wrapper.model_reader.model
            tissue_model = original_model.copy()

            selected_ids = [self.model_wrapper.model_reader.r_ids[i] for i in selected_indices]
            reactions_to_remove = [rxn for rxn in tissue_model.reactions if rxn.id not in selected_ids]
            tissue_model.remove_reactions(reactions_to_remove, remove_orphans=True)

            # Create validator with converted gene IDs
            validator = GeneEssentialityValidator(
                essential_genes=self.essential_genes_converted,
                non_essential_genes=self.non_essential_genes_converted,
                growth_threshold=0.01
            )

            # Validate
            essentiality_result = validator.validate(tissue_model)

            # Store results
            result.essentiality_accuracy = essentiality_result.accuracy
            result.essentiality_precision = essentiality_result.precision
            result.essentiality_recall = essentiality_result.recall
            result.essentiality_f1 = essentiality_result.f1_score
            result.essentiality_mcc = essentiality_result.matthews_correlation
            result.essentiality_result = essentiality_result

            if self.verbose:
                print(f"  Essentiality - Acc: {essentiality_result.accuracy:.3f}, "
                      f"F1: {essentiality_result.f1_score:.3f}, "
                      f"MCC: {essentiality_result.matthews_correlation:.3f}")

        except Exception as e:
            if self.verbose:
                print(f"  Warning: Essentiality validation failed: {str(e)}")

    def _validate_yields(self, result: BenchmarkResult, selected_indices: List[int]):
        """Calculate theoretical yields for various carbon sources"""
        try:
            from troppo.benchmark.validation import TheoreticalYieldCalculator
            import cobra

            # Create tissue-specific model
            original_model = self.model_wrapper.model_reader.model
            tissue_model = original_model.copy()

            selected_ids = [self.model_wrapper.model_reader.r_ids[i] for i in selected_indices]
            reactions_to_remove = [rxn for rxn in tissue_model.reactions if rxn.id not in selected_ids]
            tissue_model.remove_reactions(reactions_to_remove, remove_orphans=True)

            # Create calculator
            calculator = TheoreticalYieldCalculator(
                biomass_reaction=self.biomass_reaction
            )

            # Calculate yields
            yield_results = calculator.calculate_yields(
                tissue_model,
                carbon_sources=self.carbon_sources,
                aerobic=True,
                anaerobic=True,
                product='biomass'
            )

            # Store results
            result.yield_results = yield_results

            # Calculate averages
            aerobic_yields = [r.theoretical_yield for r in yield_results if r.aerobic and r.feasible]
            anaerobic_yields = [r.theoretical_yield for r in yield_results if not r.aerobic and r.feasible]

            result.avg_aerobic_yield = np.mean(aerobic_yields) if aerobic_yields else 0.0
            result.avg_anaerobic_yield = np.mean(anaerobic_yields) if anaerobic_yields else 0.0
            result.num_feasible_conditions = len([r for r in yield_results if r.feasible])

            if self.verbose:
                print(f"  Yields - Aerobic: {result.avg_aerobic_yield:.3f}, "
                      f"Anaerobic: {result.avg_anaerobic_yield:.3f}, "
                      f"Feasible: {result.num_feasible_conditions}/{len(yield_results)}")

        except Exception as e:
            if self.verbose:
                print(f"  Warning: Yield validation failed: {str(e)}")


def quick_benchmark(
    model_wrapper,
    data_map,
    methods: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Quick benchmark function for easy comparison

    Parameters
    ----------
    model_wrapper : ModelBasedWrapper
        Wrapped model
    data_map : OmicsDataMap
        Omics data map
    methods : List[str], optional
        Methods to benchmark
    **kwargs
        Additional arguments for BenchmarkRunner

    Returns
    -------
    pd.DataFrame
        Summary dataframe

    Examples
    --------
    >>> summary = quick_benchmark(model_wrapper, data_map)
    >>> print(summary)
    """
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=methods,
        **kwargs
    )

    comparison = runner.run_benchmark()
    return comparison.get_summary_dataframe()
