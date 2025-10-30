"""
Example: Benchmark with Entrez Gene IDs

This example demonstrates how to use the benchmark framework when your
gene essentiality data uses Entrez IDs instead of Ensembl IDs or gene symbols.

Entrez ID를 사용한 벤치마크 예제
"""

import pandas as pd
import cobra
import re
import numpy as np

# Troppo imports
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper
from troppo.benchmark import BenchmarkRunner
from troppo.benchmark.gene_id_utils import (
    detect_id_type,
    convert_gene_ids,
    FlexibleGeneIDHandler
)


def example_with_entrez_ids():
    """
    Complete example using Entrez IDs for gene essentiality data
    """

    print("=" * 80)
    print("Troppo Benchmark with Entrez Gene IDs")
    print("=" * 80)

    # ========================================================================
    # 1. Load Model and Data
    # ========================================================================
    print("\n[1/5] Loading model and omics data...")

    # GPR parsing
    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    # Load model
    model = cobra.io.read_sbml_model('data/HumanGEM_Consistent_COVID19_HAM.xml')
    print(f"  Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    # Check model gene ID type
    sample_genes = [g.id for g in list(model.genes)[:10]]
    print(f"  Sample model genes: {sample_genes}")

    # Load omics data
    omics_data = pd.read_csv('data/Desai-GTEx_ensembl.csv', index_col=0)
    print(f"  Omics: {omics_data.shape[0]} samples, {omics_data.shape[1]} genes")

    # ========================================================================
    # 2. Prepare Omics Data
    # ========================================================================
    print("\n[2/5] Preparing omics data...")

    reader = TabularReader(
        path_or_df=omics_data,
        nomenclature='ensemble_gene_id',
        omics_type='transcriptomics'
    )
    omics_container = reader.to_containers()[0]

    model_wrapper = ModelBasedWrapper(
        model=model,
        ttg_ratio=9999,
        gpr_gene_parse_function=replace_alt_transcripts
    )

    data_map = omics_container.get_integrated_data_map(
        model_reader=model_wrapper.model_reader,
        and_func=min,
        or_func=sum
    )

    print(f"  Gene-reaction mapping complete: {len(data_map.get_scores())} scores")

    # ========================================================================
    # 3. Prepare Gene Essentiality Data with Entrez IDs
    # ========================================================================
    print("\n[3/5] Preparing gene essentiality data (Entrez IDs)...")

    # Example: Essential genes in Entrez ID format
    # These are actual essential genes (examples)
    essential_genes_entrez = [
        '2597',   # GAPDH
        '3156',   # HMGCR
        '5230',   # PGK1
        '230',    # ALDOC
        '226',    # ALDOA
        '2821',   # GPI
        '5105',   # PCK1
        '5315',   # PKM
        '10327',  # AKR1A1
        '229',    # ALDOB
    ]

    # Example: Non-essential genes in Entrez ID format
    non_essential_genes_entrez = [
        '672',    # BRCA1
        '675',    # BRCA2
        '7157',   # TP53
        '5728',   # PTEN
        '999',    # CDH1
        '1029',   # CDKN2A
        '2064',   # ERBB2
        '5290',   # PIK3CA
        '3845',   # KRAS
        '7490',   # WT1
    ]

    print(f"  Essential genes (Entrez IDs): {len(essential_genes_entrez)}")
    print(f"  Sample: {essential_genes_entrez[:3]}")
    print(f"  Non-essential genes (Entrez IDs): {len(non_essential_genes_entrez)}")
    print(f"  Sample: {non_essential_genes_entrez[:3]}")

    # ========================================================================
    # 4. Automatic ID Detection and Conversion
    # ========================================================================
    print("\n[4/5] Demonstrating automatic ID detection...")

    # Detect ID types
    essential_type = detect_id_type(essential_genes_entrez)
    non_essential_type = detect_id_type(non_essential_genes_entrez)

    print(f"  Detected essential genes type: {essential_type}")
    print(f"  Detected non-essential genes type: {non_essential_type}")

    # Get model genes
    model_genes = [g.id for g in model.genes]
    model_gene_type = detect_id_type(model_genes[:100])
    print(f"  Detected model gene type: {model_gene_type}")

    # Manual conversion example
    if essential_type != model_gene_type:
        print(f"\n  Converting essential genes from {essential_type} to {model_gene_type}...")
        converted = convert_gene_ids(
            essential_genes_entrez,
            from_type=essential_type,
            to_type=model_gene_type
        )
        print(f"  Converted {len(converted)} genes")
        print(f"  Sample conversions:")
        for orig, conv in list(converted.items())[:3]:
            print(f"    {orig} ({essential_type}) → {conv} ({model_gene_type})")

    # ========================================================================
    # 5. Run Benchmark with Automatic ID Conversion
    # ========================================================================
    print("\n[5/5] Running benchmark with automatic ID conversion...")

    # The BenchmarkRunner will automatically convert IDs!
    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'tinit', 'imat'],
        biomass_reaction='biomass_human',
        # Essential genes in Entrez ID format
        essential_genes=essential_genes_entrez,
        essential_genes_id_type='entrez_id',  # Specify ID type
        # Non-essential genes in Entrez ID format
        non_essential_genes=non_essential_genes_entrez,
        non_essential_genes_id_type='entrez',  # Can use alias
        # Theoretical yields
        carbon_sources=['glucose', 'acetate'],
        verbose=True
    )

    print("\nBenchmark will automatically:")
    print("  1. Detect model gene ID type")
    print("  2. Convert Entrez IDs to model gene ID type")
    print("  3. Run validation with converted IDs")
    print()

    comparison = runner.run_benchmark(
        validate_biology=True,
        validate_essentiality=True,
        validate_yields=True
    )

    # ========================================================================
    # 6. View Results
    # ========================================================================
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)

    summary = comparison.get_summary_dataframe()
    print(summary)

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

    print("\n" + "=" * 80)
    print("Benchmark completed successfully with Entrez IDs!")
    print("=" * 80)


def example_flexible_id_handler():
    """
    Example using FlexibleGeneIDHandler for mixed ID types
    """
    print("\n" + "=" * 80)
    print("FlexibleGeneIDHandler Example")
    print("=" * 80)

    # Create handler
    handler = FlexibleGeneIDHandler(target_nomenclature='symbol')

    # Add gene lists with different ID types
    print("\nAdding gene lists...")

    # Essential genes in Entrez IDs
    handler.add_gene_list(
        'essential',
        ['2597', '3156', '5230'],
        id_type='entrez_id'
    )

    # Non-essential genes in Ensembl IDs
    handler.add_gene_list(
        'non_essential',
        ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000171862'],
        id_type='ensembl_gene_id'
    )

    # Model genes in symbols
    handler.add_gene_list(
        'model',
        ['GAPDH', 'HMGCR', 'PGK1', 'BRCA1', 'BRCA2'],
        id_type='symbol'
    )

    # Get converted lists
    print("\nConverted lists (all in symbol format):")
    for name in ['essential', 'non_essential', 'model']:
        genes = handler.get_converted_list(name)
        print(f"  {name}: {genes}")

    # Match lists
    matched = handler.match_lists('essential', 'model')
    print(f"\nMatched genes between essential and model: {matched}")


def example_mixed_id_types():
    """
    Example with mixed ID types in experimental data
    """
    print("\n" + "=" * 80)
    print("Mixed ID Types Example")
    print("=" * 80)

    # Simulated experimental data with mixed ID types
    experimental_data = pd.DataFrame({
        'gene_id': [
            '2597',     # Entrez
            'ENSG00000111640',  # Ensembl
            'HMGCR',    # Symbol
            '5230',     # Entrez
            'ENSG00000100320',  # Ensembl
        ],
        'essential': [True, True, True, False, False]
    })

    print("\nExperimental data with mixed ID types:")
    print(experimental_data)

    # Separate by essentiality
    essential = experimental_data[experimental_data['essential']]['gene_id'].tolist()
    non_essential = experimental_data[~experimental_data['essential']]['gene_id'].tolist()

    print(f"\nEssential genes (mixed IDs): {essential}")
    print(f"Non-essential genes (mixed IDs): {non_essential}")

    # Strategy: Let auto-detection handle it
    # The framework will try to detect and convert automatically
    print("\nNote: For mixed ID types, it's recommended to:")
    print("  1. Pre-process data to use consistent ID type")
    print("  2. Or use FlexibleGeneIDHandler to manage conversions")


if __name__ == '__main__':
    print("\n" + "=" * 80)
    print("Troppo Benchmark: Entrez ID Support Examples")
    print("=" * 80)

    # Uncomment to run examples:

    # Main example with Entrez IDs
    # example_with_entrez_ids()

    # Flexible ID handler
    # example_flexible_id_handler()

    # Mixed ID types
    # example_mixed_id_types()

    print("\n" + "=" * 80)
    print("Examples defined. Uncomment function calls to run.")
    print("=" * 80)

    print("\nSupported ID types:")
    print("  - entrez_id (NCBI Entrez Gene IDs)")
    print("  - ensembl_gene_id (Ensembl Gene IDs)")
    print("  - symbol (HGNC Gene Symbols)")
    print("  - hgnc_id (HGNC IDs)")
    print("  - uniprot_ids (UniProt IDs)")
    print("  - And more via HGNC database")

    print("\nKey features:")
    print("  ✓ Automatic ID type detection")
    print("  ✓ Automatic ID conversion using HGNC database")
    print("  ✓ Support for mixed ID types")
    print("  ✓ Flexible ID handler for complex scenarios")
    print("  ✓ Verbose logging of conversion process")
