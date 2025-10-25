"""
Example: Benchmark with Gene Essentiality and Theoretical Yield Validation

This example demonstrates how to use the enhanced benchmark framework with
biological validation metrics.

Gene essentiality와 이론적 수율 검증을 포함한 벤치마크 예제
"""

import pandas as pd
import cobra
import re
import numpy as np

# Troppo imports
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper
from troppo.benchmark import (
    BenchmarkRunner,
    GeneEssentialityValidator,
    TheoreticalYieldCalculator
)
from troppo.benchmark.viz_extended import (
    plot_essentiality_comparison,
    plot_confusion_matrices,
    plot_yield_comparison,
    plot_yield_heatmap,
    plot_comprehensive_validation
)


def main():
    """
    Complete example of benchmark with validation
    """

    print("=" * 80)
    print("Troppo Benchmark with Biological Validation Example")
    print("=" * 80)

    # ========================================================================
    # 1. Load Model and Data
    # ========================================================================
    print("\n[1/6] Loading model and omics data...")

    # GPR parsing
    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    # Load model
    model = cobra.io.read_sbml_model('data/HumanGEM_Consistent_COVID19_HAM.xml')
    print(f"  Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    # Load omics data
    omics_data = pd.read_csv('data/Desai-GTEx_ensembl.csv', index_col=0)
    print(f"  Omics: {omics_data.shape[0]} samples, {omics_data.shape[1]} genes")

    # ========================================================================
    # 2. Prepare Data
    # ========================================================================
    print("\n[2/6] Preparing data...")

    # Create omics container
    reader = TabularReader(
        path_or_df=omics_data,
        nomenclature='ensemble_gene_id',
        omics_type='transcriptomics'
    )
    omics_container = reader.to_containers()[0]

    # Create model wrapper
    model_wrapper = ModelBasedWrapper(
        model=model,
        ttg_ratio=9999,
        gpr_gene_parse_function=replace_alt_transcripts
    )

    # Map genes to reactions
    data_map = omics_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=min,
        or_func=sum
    )

    print(f"  Gene-reaction mapping complete: {len(data_map.get_scores())} scores")

    # ========================================================================
    # 3. Load Validation Data
    # ========================================================================
    print("\n[3/6] Loading validation data...")

    # Example: Load gene essentiality data
    # In real use, load from experimental data
    essential_genes_example = [
        'ENSG00000111640',  # GAPDH
        'ENSG00000100292',  # HMGCR
        'ENSG00000108479',  # GALE
        # Add more essential genes...
    ]

    non_essential_genes_example = [
        'ENSG00000139618',  # BRCA2
        'ENSG00000141510',  # TP53
        # Add more non-essential genes...
    ]

    # Example: Carbon sources to test
    carbon_sources = ['glucose', 'acetate', 'glycerol']

    print(f"  Essential genes: {len(essential_genes_example)}")
    print(f"  Non-essential genes: {len(non_essential_genes_example)}")
    print(f"  Carbon sources: {carbon_sources}")

    # ========================================================================
    # 4. Run Comprehensive Benchmark
    # ========================================================================
    print("\n[4/6] Running comprehensive benchmark...")
    print("  This includes:")
    print("    - Performance metrics (time, memory)")
    print("    - Model quality (network consistency)")
    print("    - Gene essentiality validation")
    print("    - Theoretical yields (aerobic & anaerobic)")

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat', 'fastcore'],
        biomass_reaction='biomass_human',
        essential_genes=essential_genes_example,
        non_essential_genes=non_essential_genes_example,
        carbon_sources=carbon_sources,
        verbose=True
    )

    comparison = runner.run_benchmark(
        validate_biology=True,
        validate_consistency=True,
        validate_essentiality=True,
        validate_yields=True
    )

    # ========================================================================
    # 5. Analyze Results
    # ========================================================================
    print("\n[5/6] Analyzing results...")

    # Overall summary
    print("\n" + "=" * 80)
    print("OVERALL SUMMARY")
    print("=" * 80)
    print(comparison.get_summary_dataframe())

    # Gene essentiality summary
    print("\n" + "=" * 80)
    print("GENE ESSENTIALITY VALIDATION")
    print("=" * 80)
    for method_name, result in comparison.results.items():
        if result.success and result.essentiality_accuracy is not None:
            print(f"\n{method_name.upper()}:")
            print(f"  Accuracy:  {result.essentiality_accuracy:.3f}")
            print(f"  Precision: {result.essentiality_precision:.3f}")
            print(f"  Recall:    {result.essentiality_recall:.3f}")
            print(f"  F1-Score:  {result.essentiality_f1:.3f}")
            print(f"  MCC:       {result.essentiality_mcc:.3f}")

    # Theoretical yields summary
    print("\n" + "=" * 80)
    print("THEORETICAL YIELDS")
    print("=" * 80)
    for method_name, result in comparison.results.items():
        if result.success and result.yield_results:
            print(f"\n{method_name.upper()}:")
            print(f"  Avg Aerobic Yield:   {result.avg_aerobic_yield:.3f}")
            print(f"  Avg Anaerobic Yield: {result.avg_anaerobic_yield:.3f}")
            print(f"  Feasible Conditions: {result.num_feasible_conditions}/{len(result.yield_results)}")

    # Detailed yield table
    print("\n" + "-" * 80)
    print("Detailed Yields by Carbon Source:")
    print("-" * 80)
    yield_data = []
    for method_name, result in comparison.results.items():
        if result.success and result.yield_results:
            for yr in result.yield_results:
                yield_data.append({
                    'Method': method_name,
                    'Carbon Source': yr.carbon_source,
                    'Condition': 'Aerobic' if yr.aerobic else 'Anaerobic',
                    'Yield': yr.theoretical_yield if yr.feasible else 0.0,
                    'Feasible': yr.feasible
                })

    yield_df = pd.DataFrame(yield_data)
    yield_pivot = yield_df.pivot_table(
        index=['Carbon Source', 'Condition'],
        columns='Method',
        values='Yield',
        aggfunc='first'
    )
    print(yield_pivot)

    # ========================================================================
    # 6. Visualizations
    # ========================================================================
    print("\n[6/6] Generating visualizations...")

    try:
        # Essentiality comparison
        print("  - Gene essentiality comparison...")
        plot_essentiality_comparison(
            comparison,
            save_path='essentiality_comparison.png'
        )

        # Confusion matrices
        print("  - Confusion matrices...")
        plot_confusion_matrices(
            comparison,
            save_path='confusion_matrices.png'
        )

        # Yield comparison
        print("  - Theoretical yields comparison...")
        plot_yield_comparison(
            comparison,
            save_path='yield_comparison.png'
        )

        # Yield heatmaps
        print("  - Aerobic yield heatmap...")
        plot_yield_heatmap(
            comparison,
            condition='aerobic',
            save_path='yield_heatmap_aerobic.png'
        )

        print("  - Anaerobic yield heatmap...")
        plot_yield_heatmap(
            comparison,
            condition='anaerobic',
            save_path='yield_heatmap_anaerobic.png'
        )

        # Comprehensive validation
        print("  - Comprehensive validation plot...")
        plot_comprehensive_validation(
            comparison,
            save_path='comprehensive_validation.png'
        )

    except Exception as e:
        print(f"  Warning: Some visualizations failed: {str(e)}")

    # ========================================================================
    # Save Results
    # ========================================================================
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    comparison.save_to_json('benchmark_results_with_validation.json')
    print("  ✓ Results saved to benchmark_results_with_validation.json")

    comparison.get_summary_dataframe().to_csv('benchmark_summary.csv', index=False)
    print("  ✓ Summary saved to benchmark_summary.csv")

    yield_pivot.to_csv('theoretical_yields.csv')
    print("  ✓ Yields saved to theoretical_yields.csv")

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE!")
    print("=" * 80)


# ============================================================================
# Standalone validation examples
# ============================================================================

def example_gene_essentiality_only():
    """
    Example of using gene essentiality validator standalone
    """
    print("\n" + "=" * 80)
    print("Standalone Gene Essentiality Validation Example")
    print("=" * 80)

    # Load model
    model = cobra.io.read_sbml_model('data/HumanGEM_Consistent_COVID19_HAM.xml')

    # Create validator
    validator = GeneEssentialityValidator(
        essential_genes=['ENSG00000111640', 'ENSG00000100292'],
        non_essential_genes=['ENSG00000139618', 'ENSG00000141510'],
        growth_threshold=0.01
    )

    # Validate
    result = validator.validate(model)

    print(f"\nResults:")
    print(f"  Accuracy:  {result.accuracy:.3f}")
    print(f"  Precision: {result.precision:.3f}")
    print(f"  Recall:    {result.recall:.3f}")
    print(f"  F1-Score:  {result.f1_score:.3f}")
    print(f"  MCC:       {result.matthews_correlation:.3f}")

    print(f"\nConfusion Matrix:")
    print(f"  True Positives:  {len(result.true_positives)}")
    print(f"  False Positives: {len(result.false_positives)}")
    print(f"  True Negatives:  {len(result.true_negatives)}")
    print(f"  False Negatives: {len(result.false_negatives)}")


def example_theoretical_yield_only():
    """
    Example of using theoretical yield calculator standalone
    """
    print("\n" + "=" * 80)
    print("Standalone Theoretical Yield Calculation Example")
    print("=" * 80)

    # Load model
    model = cobra.io.read_sbml_model('data/HumanGEM_Consistent_COVID19_HAM.xml')

    # Create calculator
    calculator = TheoreticalYieldCalculator(
        biomass_reaction='biomass_human'
    )

    # Calculate yields
    results = calculator.calculate_yields(
        model,
        carbon_sources=['glucose', 'acetate', 'glycerol'],
        aerobic=True,
        anaerobic=True
    )

    print("\nTheoretical Yields:")
    print("-" * 80)
    for result in results:
        condition = "Aerobic" if result.aerobic else "Anaerobic"
        status = "✓" if result.feasible else "✗"
        print(f"{status} {result.carbon_source:10} ({condition:10}): "
              f"{result.theoretical_yield:.4f} gDW/mmol")

    # Convert to DataFrame
    df = calculator.results_to_dataframe(results)
    print("\nDataFrame:")
    print(df[['carbon_source', 'aerobic', 'theoretical_yield', 'feasible']])


if __name__ == '__main__':
    # Run main example
    # Uncomment to run:
    # main()

    # Or run standalone examples:
    # example_gene_essentiality_only()
    # example_theoretical_yield_only()

    print("\nTo run this example:")
    print("  1. Ensure you have the required data files")
    print("  2. Uncomment the function calls at the end of the script")
    print("  3. Run: python benchmark_with_validation_example.py")
