"""
Example: How to Create and Register a Custom Omics Integration Method

This example demonstrates how to create a new custom method and integrate it
with the Troppo framework using the method registry system.

새로운 커스텀 방법론을 만들고 Troppo 프레임워크에 통합하는 방법을 보여줍니다.
"""

import numpy as np
from troppo.methods.base import ContextSpecificModelReconstructionAlgorithm
from troppo.methods.registry import MethodRegistry, register_method
from cobamp.utilities.property_management import PropertyDictionary


# ============================================================================
# Step 1: Define Properties Class
# ============================================================================

class CustomMethodProperties(PropertyDictionary):
    """
    Properties for custom reconstruction method

    This class defines the parameters needed by your method.
    """

    def __init__(
        self,
        expression_scores: dict,
        threshold: float = 0.5,
        weight_factor: float = 1.0,
        solver: str = 'CPLEX',
        reaction_ids: list = None,
        metabolite_ids: list = None,
        **kwargs
    ):
        """
        Initialize properties

        Parameters
        ----------
        expression_scores : dict
            Dictionary of {reaction_id: expression_score}
        threshold : float
            Expression threshold for reaction selection
        weight_factor : float
            Weight factor for scoring function
        solver : str
            LP/MILP solver to use
        reaction_ids : list
            List of reaction IDs
        metabolite_ids : list
            List of metabolite IDs
        """
        # Define mandatory and optional properties
        mandatory = {
            'expression_scores': dict,
            'solver': str
        }

        optional = {
            'threshold': float,
            'weight_factor': float,
            'reaction_ids': list,
            'metabolite_ids': list
        }

        super().__init__(mandatory, optional)

        # Set properties
        self.expression_scores = expression_scores
        self.threshold = threshold
        self.weight_factor = weight_factor
        self.solver = solver
        self.reaction_ids = reaction_ids or []
        self.metabolite_ids = metabolite_ids or []

        # Store additional kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)


# ============================================================================
# Step 2: Implement the Method Class
# ============================================================================

class CustomMethod(ContextSpecificModelReconstructionAlgorithm):
    """
    Custom omics integration method

    This is a simple example method that selects reactions based on
    a weighted scoring function combining expression data and network properties.

    간단한 예제 방법론으로, 발현 데이터와 네트워크 속성을 결합한 가중 점수 함수를
    기반으로 반응을 선택합니다.

    Algorithm:
    1. Calculate weighted score for each reaction
    2. Select reactions above threshold
    3. Ensure network consistency using FVA
    """

    def __init__(self, S, lb, ub, properties: CustomMethodProperties):
        """
        Initialize custom method

        Parameters
        ----------
        S : array-like
            Stoichiometric matrix
        lb : array-like
            Lower bounds
        ub : array-like
            Upper bounds
        properties : CustomMethodProperties
            Method properties
        """
        self.S = S
        self.lb = lb
        self.ub = ub
        self.properties = properties

        # Store solver
        self.solver = properties.solver

    def run(self):
        """
        Execute the custom method

        Returns
        -------
        list
            List of selected reaction indices
        """
        print(f"Running CustomMethod with threshold={self.properties.threshold}")

        # Step 1: Calculate scores
        scores = self._calculate_scores()

        # Step 2: Select reactions above threshold
        selected = self._select_reactions(scores)

        # Step 3: Ensure consistency (optional)
        selected = self._ensure_consistency(selected)

        print(f"Selected {len(selected)} reactions")

        return selected

    def _calculate_scores(self):
        """
        Calculate weighted scores for reactions

        This is where you implement your scoring logic.
        """
        scores = []
        expression_scores = self.properties.expression_scores
        weight = self.properties.weight_factor

        for idx, rxn_id in enumerate(self.properties.reaction_ids):
            # Get expression score (default to 0 if not found)
            expr_score = expression_scores.get(rxn_id, 0)

            # Calculate network-based score (example: based on bounds)
            network_score = self._calculate_network_score(idx)

            # Combine scores
            combined_score = weight * expr_score + (1 - weight) * network_score
            scores.append(combined_score)

        return np.array(scores)

    def _calculate_network_score(self, reaction_idx):
        """
        Calculate network-based score for a reaction

        This is an example - you can implement any scoring logic here.
        """
        # Example: reactions with larger flux ranges get higher scores
        lb = self.lb[reaction_idx]
        ub = self.ub[reaction_idx]

        flux_range = abs(ub - lb)
        # Normalize to 0-1
        max_range = 1000  # Example maximum
        return min(flux_range / max_range, 1.0)

    def _select_reactions(self, scores):
        """Select reactions above threshold"""
        threshold = self.properties.threshold
        selected_indices = np.where(scores >= threshold)[0].tolist()
        return selected_indices

    def _ensure_consistency(self, selected_indices):
        """
        Ensure selected reactions form a consistent network

        This is a simplified example. In practice, you might want to:
        - Check for dead-end metabolites
        - Ensure biomass production is possible
        - Remove blocked reactions
        """
        # For this example, just return the selected reactions
        # In a real implementation, you would perform consistency checks
        return selected_indices

    @property
    def properties_class(self):
        """Return the properties class"""
        return CustomMethodProperties


# ============================================================================
# Step 3: Register the Method
# ============================================================================

# Option 1: Manual registration
MethodRegistry.register(
    name='custom_method',
    method_class=CustomMethod,
    properties_class=CustomMethodProperties,
    description='Custom method combining expression and network scores',
    category='reconstruction',
    requires_continuous_scores=True,
    speed='fast',
    memory='low',
    complexity='low',
    reference='Your Name et al. (2024)'
)

print("✓ CustomMethod registered successfully!")


# ============================================================================
# Step 4: Use the Method
# ============================================================================

def example_usage():
    """
    Example of using the custom method
    """
    import pandas as pd
    import cobra
    from troppo.omics.readers.generic import TabularReader
    from troppo.methods_wrappers import ModelBasedWrapper
    from troppo.omics.integration import ContinuousScoreIntegrationStrategy

    # Load model and data (example)
    # model = cobra.io.read_sbml_model('path/to/model.xml')
    # omics_data = pd.read_csv('path/to/data.csv', index_col=0)

    # Create containers (example)
    # reader = TabularReader(path_or_df=omics_data, ...)
    # omics_container = reader.to_containers()[0]
    # model_wrapper = ModelBasedWrapper(model=model, ttg_ratio=9999)

    # Get integrated scores
    # data_map = omics_container.get_integrated_data_map(...)

    # Prepare scores for custom method
    # integration = ContinuousScoreIntegrationStrategy()
    # scores = integration.integrate(data_map=data_map)

    # Create properties
    properties = CustomMethodProperties(
        expression_scores={},  # Your scores here
        threshold=0.5,
        weight_factor=0.7,
        solver='CPLEX',
        reaction_ids=[],  # Your reaction IDs
        metabolite_ids=[]  # Your metabolite IDs
    )

    # Create and run method using registry
    # custom = MethodRegistry.create_method(
    #     'custom_method',
    #     S=model_wrapper.S,
    #     lb=model_wrapper.lb,
    #     ub=model_wrapper.ub,
    #     properties=properties
    # )
    #
    # result = custom.run()

    print("Example usage completed!")


# ============================================================================
# Step 5: Use in Benchmark
# ============================================================================

def example_benchmark():
    """
    Example of including custom method in benchmark
    """
    from troppo.benchmark import BenchmarkRunner

    # Your model_wrapper and data_map here

    # Run benchmark including your custom method
    # runner = BenchmarkRunner(
    #     model_wrapper=model_wrapper,
    #     data_map=data_map,
    #     methods=['gimme', 'tinit', 'imat', 'custom_method']  # Include your method!
    # )
    #
    # comparison = runner.run_benchmark()
    # print(comparison.get_summary_dataframe())

    print("Benchmark example completed!")


# ============================================================================
# Alternative: Using Decorator
# ============================================================================

@register_method(
    name='custom_method_v2',
    description='Custom method v2 using decorator',
    requires_continuous_scores=True,
    speed='fast'
)
class CustomMethodV2(ContextSpecificModelReconstructionAlgorithm):
    """
    Alternative way to register a method using decorator
    """

    properties_class = CustomMethodProperties  # Required for decorator

    def __init__(self, S, lb, ub, properties):
        self.S = S
        self.lb = lb
        self.ub = ub
        self.properties = properties

    def run(self):
        # Your implementation here
        return list(range(len(self.properties.reaction_ids) // 2))


print("✓ CustomMethodV2 registered via decorator!")


# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    print("=" * 70)
    print("Custom Method Registration Example")
    print("=" * 70)

    # Show all registered methods
    print("\nRegistered Methods:")
    MethodRegistry.print_registry()

    print("\n" + "=" * 70)
    print("Example completed! Your custom method is ready to use.")
    print("=" * 70)
