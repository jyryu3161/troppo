"""
Biological Validation Modules for Benchmark Framework

This module provides tools for validating reconstructed models using:
1. Gene essentiality data
2. Theoretical yield calculations for various carbon sources

생물학적 검증 모듈: 유전자 필수성 및 이론적 수율 계산
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Optional, Tuple, Any
from dataclasses import dataclass, field
import warnings


@dataclass
class GeneEssentialityResult:
    """
    Results from gene essentiality analysis
    유전자 필수성 분석 결과
    """
    # Input data
    experimental_essential: Set[str] = field(default_factory=set)
    experimental_non_essential: Set[str] = field(default_factory=set)

    # Predictions
    predicted_essential: Set[str] = field(default_factory=set)
    predicted_non_essential: Set[str] = field(default_factory=set)

    # Confusion matrix components
    true_positives: Set[str] = field(default_factory=set)
    false_positives: Set[str] = field(default_factory=set)
    true_negatives: Set[str] = field(default_factory=set)
    false_negatives: Set[str] = field(default_factory=set)

    # Metrics
    accuracy: float = 0.0
    precision: float = 0.0
    recall: float = 0.0
    f1_score: float = 0.0
    specificity: float = 0.0
    matthews_correlation: float = 0.0

    # Additional info
    num_tested: int = 0
    num_experimental_essential: int = 0
    num_predicted_essential: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary (excluding gene sets for compactness)"""
        return {
            'accuracy': self.accuracy,
            'precision': self.precision,
            'recall': self.recall,
            'f1_score': self.f1_score,
            'specificity': self.specificity,
            'matthews_correlation': self.matthews_correlation,
            'num_tested': self.num_tested,
            'num_experimental_essential': self.num_experimental_essential,
            'num_predicted_essential': self.num_predicted_essential,
            'tp': len(self.true_positives),
            'fp': len(self.false_positives),
            'tn': len(self.true_negatives),
            'fn': len(self.false_negatives)
        }


@dataclass
class TheoreticalYieldResult:
    """
    Results from theoretical yield calculations
    이론적 수율 계산 결과
    """
    carbon_source: str
    aerobic: bool

    # Yield values
    theoretical_yield: float = 0.0  # mol product / mol substrate
    biomass_yield: float = 0.0      # gDW / mol substrate

    # Flux values
    substrate_uptake: float = 0.0
    product_formation: float = 0.0
    oxygen_uptake: float = 0.0
    co2_production: float = 0.0

    # Solution status
    feasible: bool = False
    optimal: bool = False

    # Metabolite IDs used
    substrate_id: str = ""
    product_id: str = ""
    oxygen_id: str = ""
    co2_id: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'carbon_source': self.carbon_source,
            'aerobic': self.aerobic,
            'theoretical_yield': self.theoretical_yield,
            'biomass_yield': self.biomass_yield,
            'substrate_uptake': self.substrate_uptake,
            'product_formation': self.product_formation,
            'oxygen_uptake': self.oxygen_uptake if self.aerobic else 0.0,
            'co2_production': self.co2_production,
            'feasible': self.feasible,
            'optimal': self.optimal
        }


class GeneEssentialityValidator:
    """
    Validate model predictions against experimental gene essentiality data

    Examples
    --------
    >>> validator = GeneEssentialityValidator(
    ...     essential_genes=['gene1', 'gene2'],
    ...     non_essential_genes=['gene3', 'gene4']
    ... )
    >>> result = validator.validate(model, minimal_medium)
    >>> print(f"Accuracy: {result.accuracy:.3f}")
    """

    def __init__(
        self,
        essential_genes: Optional[List[str]] = None,
        non_essential_genes: Optional[List[str]] = None,
        essentiality_data: Optional[pd.DataFrame] = None,
        growth_threshold: float = 0.01
    ):
        """
        Initialize validator

        Parameters
        ----------
        essential_genes : List[str], optional
            List of experimentally verified essential genes
        non_essential_genes : List[str], optional
            List of experimentally verified non-essential genes
        essentiality_data : pd.DataFrame, optional
            DataFrame with columns ['gene_id', 'essential'] where essential is boolean
        growth_threshold : float
            Growth rate threshold below which gene is considered essential (default: 0.01)
        """
        self.growth_threshold = growth_threshold

        if essentiality_data is not None:
            self.experimental_essential = set(
                essentiality_data[essentiality_data['essential']]['gene_id'].tolist()
            )
            self.experimental_non_essential = set(
                essentiality_data[~essentiality_data['essential']]['gene_id'].tolist()
            )
        else:
            self.experimental_essential = set(essential_genes or [])
            self.experimental_non_essential = set(non_essential_genes or [])

    def validate(
        self,
        model,
        condition: Optional[Dict[str, float]] = None,
        method: str = 'single_deletion'
    ) -> GeneEssentialityResult:
        """
        Validate gene essentiality predictions

        Parameters
        ----------
        model : cobra.Model
            Metabolic model to test
        condition : Dict[str, float], optional
            Medium conditions (exchange reaction bounds)
        method : str
            Method to use ('single_deletion', 'moma', 'room')

        Returns
        -------
        GeneEssentialityResult
            Validation results
        """
        try:
            import cobra
            from cobra.flux_analysis import single_gene_deletion
        except ImportError:
            raise ImportError("COBRApy required for gene essentiality validation")

        # Apply condition if provided
        if condition:
            with model:
                for rxn_id, bound in condition.items():
                    if rxn_id in model.reactions:
                        rxn = model.reactions.get_by_id(rxn_id)
                        rxn.lower_bound = min(bound, 0)
                        rxn.upper_bound = max(bound, 0)

        # Get wild-type growth rate
        wt_solution = model.optimize()
        if wt_solution.status != 'optimal':
            warnings.warn("Wild-type model has no feasible solution")
            return GeneEssentialityResult()

        wt_growth = wt_solution.objective_value

        # Perform single gene deletion
        deletion_results = single_gene_deletion(model, method=method)

        # Classify genes as essential or non-essential based on growth
        predicted_essential = set()
        predicted_non_essential = set()

        for gene_id in model.genes:
            gene_id_str = str(gene_id.id)

            # Find deletion result for this gene
            gene_result = deletion_results[deletion_results['ids'] == gene_id_str]

            if len(gene_result) > 0:
                growth = gene_result['growth'].values[0]

                # Check if growth is significantly reduced
                if growth < self.growth_threshold * wt_growth:
                    predicted_essential.add(gene_id_str)
                else:
                    predicted_non_essential.add(gene_id_str)

        # Calculate metrics
        result = self._calculate_metrics(
            predicted_essential,
            predicted_non_essential
        )

        return result

    def _calculate_metrics(
        self,
        predicted_essential: Set[str],
        predicted_non_essential: Set[str]
    ) -> GeneEssentialityResult:
        """Calculate performance metrics"""

        # Find genes in common between experimental and predicted
        tested_genes = (
            (self.experimental_essential | self.experimental_non_essential) &
            (predicted_essential | predicted_non_essential)
        )

        # Calculate confusion matrix
        tp = self.experimental_essential & predicted_essential & tested_genes
        fp = self.experimental_non_essential & predicted_essential & tested_genes
        tn = self.experimental_non_essential & predicted_non_essential & tested_genes
        fn = self.experimental_essential & predicted_non_essential & tested_genes

        # Calculate metrics
        tp_count = len(tp)
        fp_count = len(fp)
        tn_count = len(tn)
        fn_count = len(fn)

        total = tp_count + fp_count + tn_count + fn_count

        if total == 0:
            warnings.warn("No genes in common between experimental and predicted")
            return GeneEssentialityResult()

        accuracy = (tp_count + tn_count) / total if total > 0 else 0.0
        precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
        recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
        specificity = tn_count / (tn_count + fp_count) if (tn_count + fp_count) > 0 else 0.0

        f1_score = (
            2 * precision * recall / (precision + recall)
            if (precision + recall) > 0 else 0.0
        )

        # Matthews Correlation Coefficient
        numerator = (tp_count * tn_count) - (fp_count * fn_count)
        denominator = np.sqrt(
            (tp_count + fp_count) *
            (tp_count + fn_count) *
            (tn_count + fp_count) *
            (tn_count + fn_count)
        )
        mcc = numerator / denominator if denominator > 0 else 0.0

        return GeneEssentialityResult(
            experimental_essential=self.experimental_essential,
            experimental_non_essential=self.experimental_non_essential,
            predicted_essential=predicted_essential,
            predicted_non_essential=predicted_non_essential,
            true_positives=tp,
            false_positives=fp,
            true_negatives=tn,
            false_negatives=fn,
            accuracy=accuracy,
            precision=precision,
            recall=recall,
            f1_score=f1_score,
            specificity=specificity,
            matthews_correlation=mcc,
            num_tested=total,
            num_experimental_essential=len(self.experimental_essential & tested_genes),
            num_predicted_essential=len(predicted_essential & tested_genes)
        )


class TheoreticalYieldCalculator:
    """
    Calculate theoretical yields for various carbon sources under different conditions

    Examples
    --------
    >>> calculator = TheoreticalYieldCalculator()
    >>> results = calculator.calculate_yields(
    ...     model,
    ...     carbon_sources=['glucose', 'acetate'],
    ...     conditions=['aerobic', 'anaerobic']
    ... )
    """

    # Default metabolite IDs (can be overridden)
    DEFAULT_METABOLITES = {
        'glucose': {'id': 'glc__D_e', 'name': 'D-Glucose'},
        'acetate': {'id': 'ac_e', 'name': 'Acetate'},
        'glycerol': {'id': 'glyc_e', 'name': 'Glycerol'},
        'lactate': {'id': 'lac__L_e', 'name': 'L-Lactate'},
        'pyruvate': {'id': 'pyr_e', 'name': 'Pyruvate'},
        'succinate': {'id': 'succ_e', 'name': 'Succinate'},
        'xylose': {'id': 'xyl__D_e', 'name': 'D-Xylose'},
        'fructose': {'id': 'fru_e', 'name': 'Fructose'},
        'galactose': {'id': 'gal_e', 'name': 'Galactose'},
        'mannose': {'id': 'man_e', 'name': 'Mannose'},
        'oxygen': {'id': 'o2_e', 'name': 'Oxygen'},
        'co2': {'id': 'co2_e', 'name': 'CO2'}
    }

    def __init__(
        self,
        metabolite_mapping: Optional[Dict[str, Dict[str, str]]] = None,
        uptake_rate: float = 10.0,
        biomass_reaction: Optional[str] = None
    ):
        """
        Initialize calculator

        Parameters
        ----------
        metabolite_mapping : Dict, optional
            Custom metabolite ID mapping
        uptake_rate : float
            Default uptake rate for carbon sources (mmol/gDW/h)
        biomass_reaction : str, optional
            ID of biomass reaction
        """
        self.metabolite_mapping = metabolite_mapping or self.DEFAULT_METABOLITES
        self.uptake_rate = uptake_rate
        self.biomass_reaction = biomass_reaction

    def calculate_yields(
        self,
        model,
        carbon_sources: Optional[List[str]] = None,
        aerobic: bool = True,
        anaerobic: bool = True,
        product: str = 'biomass'
    ) -> List[TheoreticalYieldResult]:
        """
        Calculate theoretical yields for multiple carbon sources

        Parameters
        ----------
        model : cobra.Model
            Metabolic model
        carbon_sources : List[str], optional
            List of carbon sources to test (default: glucose, acetate, glycerol)
        aerobic : bool
            Include aerobic conditions
        anaerobic : bool
            Include anaerobic conditions
        product : str
            Product to optimize ('biomass', 'atp', or metabolite ID)

        Returns
        -------
        List[TheoreticalYieldResult]
            Results for each condition
        """
        if carbon_sources is None:
            carbon_sources = ['glucose', 'acetate', 'glycerol']

        results = []

        for carbon_source in carbon_sources:
            if aerobic:
                result = self._calculate_single_yield(
                    model, carbon_source, aerobic=True, product=product
                )
                results.append(result)

            if anaerobic:
                result = self._calculate_single_yield(
                    model, carbon_source, aerobic=False, product=product
                )
                results.append(result)

        return results

    def _calculate_single_yield(
        self,
        model,
        carbon_source: str,
        aerobic: bool,
        product: str
    ) -> TheoreticalYieldResult:
        """Calculate yield for a single condition"""

        try:
            import cobra
        except ImportError:
            raise ImportError("COBRApy required for yield calculations")

        # Get metabolite IDs
        if carbon_source not in self.metabolite_mapping:
            warnings.warn(f"Unknown carbon source: {carbon_source}")
            return TheoreticalYieldResult(carbon_source=carbon_source, aerobic=aerobic)

        substrate_id = self.metabolite_mapping[carbon_source]['id']
        oxygen_id = self.metabolite_mapping['oxygen']['id']
        co2_id = self.metabolite_mapping['co2']['id']

        # Find exchange reactions
        substrate_exchange = self._find_exchange_reaction(model, substrate_id)
        oxygen_exchange = self._find_exchange_reaction(model, oxygen_id)
        co2_exchange = self._find_exchange_reaction(model, co2_id)

        if substrate_exchange is None:
            warnings.warn(f"No exchange reaction found for {carbon_source}")
            return TheoreticalYieldResult(carbon_source=carbon_source, aerobic=aerobic)

        # Set up model conditions
        with model:
            # Close all exchange reactions
            for rxn in model.exchanges:
                rxn.lower_bound = 0
                rxn.upper_bound = 1000

            # Set substrate uptake
            substrate_exchange.lower_bound = -self.uptake_rate
            substrate_exchange.upper_bound = 0

            # Set oxygen availability
            if aerobic and oxygen_exchange:
                oxygen_exchange.lower_bound = -1000
            elif oxygen_exchange:
                oxygen_exchange.lower_bound = 0

            # Allow CO2 secretion
            if co2_exchange:
                co2_exchange.lower_bound = -1000
                co2_exchange.upper_bound = 1000

            # Optimize for product
            if product == 'biomass':
                # Use existing objective (usually biomass)
                solution = model.optimize()
            else:
                # Set custom objective
                product_exchange = self._find_exchange_reaction(model, product)
                if product_exchange:
                    model.objective = product_exchange.id
                    solution = model.optimize()
                else:
                    warnings.warn(f"Product exchange reaction not found: {product}")
                    return TheoreticalYieldResult(carbon_source=carbon_source, aerobic=aerobic)

            # Extract results
            if solution.status == 'optimal':
                substrate_flux = substrate_exchange.flux
                biomass_flux = solution.objective_value

                # Calculate yield
                if abs(substrate_flux) > 1e-6:
                    theoretical_yield = biomass_flux / abs(substrate_flux)
                else:
                    theoretical_yield = 0.0

                # Get other fluxes
                oxygen_flux = oxygen_exchange.flux if oxygen_exchange else 0.0
                co2_flux = co2_exchange.flux if co2_exchange else 0.0

                return TheoreticalYieldResult(
                    carbon_source=carbon_source,
                    aerobic=aerobic,
                    theoretical_yield=theoretical_yield,
                    biomass_yield=theoretical_yield,
                    substrate_uptake=abs(substrate_flux),
                    product_formation=biomass_flux,
                    oxygen_uptake=abs(oxygen_flux) if oxygen_flux < 0 else 0.0,
                    co2_production=co2_flux if co2_flux > 0 else 0.0,
                    feasible=True,
                    optimal=True,
                    substrate_id=substrate_id,
                    product_id=product,
                    oxygen_id=oxygen_id if aerobic else "",
                    co2_id=co2_id
                )
            else:
                return TheoreticalYieldResult(
                    carbon_source=carbon_source,
                    aerobic=aerobic,
                    feasible=False,
                    optimal=False,
                    substrate_id=substrate_id
                )

    def _find_exchange_reaction(self, model, metabolite_id: str):
        """Find exchange reaction for a metabolite"""
        try:
            # Try to find metabolite
            if metabolite_id in model.metabolites:
                met = model.metabolites.get_by_id(metabolite_id)

                # Find exchange reactions
                for rxn in met.reactions:
                    if rxn.id.startswith('EX_') or '_EX_' in rxn.id or rxn in model.exchanges:
                        return rxn

            # Alternative: search by pattern
            for rxn in model.exchanges:
                if metabolite_id.replace('_e', '') in rxn.id:
                    return rxn

            return None
        except:
            return None

    def results_to_dataframe(self, results: List[TheoreticalYieldResult]) -> pd.DataFrame:
        """Convert results to pandas DataFrame"""
        data = []
        for result in results:
            data.append(result.to_dict())
        return pd.DataFrame(data)


def compare_yields_across_methods(
    results_dict: Dict[str, List[TheoreticalYieldResult]],
    metric: str = 'theoretical_yield'
) -> pd.DataFrame:
    """
    Compare yields across different reconstruction methods

    Parameters
    ----------
    results_dict : Dict[str, List[TheoreticalYieldResult]]
        Dictionary mapping method names to their yield results
    metric : str
        Metric to compare ('theoretical_yield', 'biomass_yield', etc.)

    Returns
    -------
    pd.DataFrame
        Comparison table
    """
    comparison_data = []

    for method_name, results in results_dict.items():
        for result in results:
            comparison_data.append({
                'method': method_name,
                'carbon_source': result.carbon_source,
                'condition': 'Aerobic' if result.aerobic else 'Anaerobic',
                metric: getattr(result, metric),
                'feasible': result.feasible
            })

    df = pd.DataFrame(comparison_data)

    # Pivot for easier comparison
    pivot = df.pivot_table(
        index=['carbon_source', 'condition'],
        columns='method',
        values=metric,
        aggfunc='first'
    )

    return pivot


def create_essentiality_confusion_matrix(result: GeneEssentialityResult) -> pd.DataFrame:
    """
    Create confusion matrix as DataFrame

    Parameters
    ----------
    result : GeneEssentialityResult
        Essentiality validation result

    Returns
    -------
    pd.DataFrame
        Confusion matrix
    """
    return pd.DataFrame(
        [
            [len(result.true_positives), len(result.false_negatives)],
            [len(result.false_positives), len(result.true_negatives)]
        ],
        index=['Predicted Essential', 'Predicted Non-Essential'],
        columns=['Actually Essential', 'Actually Non-Essential']
    )
