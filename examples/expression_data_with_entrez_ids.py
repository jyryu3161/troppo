"""
Example: Using Expression Data with Entrez IDs for Omics Integration and Benchmarking

Entrez ID 형식의 발현 데이터를 사용한 오믹스 통합 및 벤치마킹 예제

This example demonstrates how to use expression data with Entrez gene IDs
(or any other nomenclature) for tissue-specific model reconstruction and
benchmarking, with automatic ID conversion to match the metabolic model.

이 예제는 Entrez 유전자 ID(또는 다른 명명법)를 사용하는 발현 데이터로
조직 특이적 모델 재구성 및 벤치마킹을 수행하는 방법을 보여줍니다.
모델의 유전자 ID에 맞춰 자동 ID 변환이 수행됩니다.

Key Features Demonstrated:
- Automatic gene ID type detection
- Automatic ID conversion to match model
- Integration with benchmark framework
- Support for multiple ID nomenclatures
"""

import numpy as np
from troppo.methods_wrappers import ModelBasedWrapper
from troppo.omics import OmicsContainer, create_data_map_from_dict
from troppo.benchmark import BenchmarkRunner
import warnings

warnings.filterwarnings('ignore')


def load_example_model():
    """
    Load an example metabolic model

    예제 대사 모델을 불러옵니다

    Replace this with your own model loading code.
    예제 모델 로딩 코드를 실제 모델로 교체하세요.
    """
    import cobra
    from cobra.test import create_test_model

    print("Loading metabolic model...")
    print("대사 모델 로딩 중...")

    # Use E. coli core model as example
    # 예제로 E. coli core 모델 사용
    model = create_test_model('textbook')

    print(f"  Model: {model.id}")
    print(f"  Reactions: {len(model.reactions)}")
    print(f"  Metabolites: {len(model.metabolites)}")
    print(f"  Genes: {len(model.genes)}")

    # Wrap the model
    model_wrapper = ModelBasedWrapper(model)

    return model_wrapper


def create_example_expression_data_entrez():
    """
    Create example expression data with Entrez IDs

    Entrez ID를 사용하는 예제 발현 데이터를 생성합니다

    Returns
    -------
    dict
        Expression data as {entrez_id: expression_value}
    """
    print("\nCreating example expression data with Entrez IDs...")
    print("Entrez ID 형식의 예제 발현 데이터 생성 중...")

    # Example: Human glycolysis and central metabolism genes (Entrez IDs)
    # 예제: 인간 해당과정 및 중심 대사 유전자들 (Entrez ID)
    expression_data_entrez = {
        # Glycolysis genes (높은 발현)
        '2597': 100.5,    # GAPDH - Glyceraldehyde-3-phosphate dehydrogenase
        '5230': 95.2,     # PGK1 - Phosphoglycerate kinase 1
        '5315': 88.7,     # PKM - Pyruvate kinase M
        '226': 92.3,      # ALDOA - Aldolase A
        '2821': 85.1,     # GPI - Glucose-6-phosphate isomerase
        '5213': 78.9,     # PFKM - Phosphofructokinase, muscle
        '5223': 81.4,     # PGAM1 - Phosphoglycerate mutase 1
        '2027': 87.6,     # ENO3 - Enolase 3

        # TCA cycle genes (중간 발현)
        '1431': 65.3,     # CS - Citrate synthase
        '47': 58.7,       # ACLY - ATP citrate lyase
        '3417': 62.1,     # IDH1 - Isocitrate dehydrogenase 1
        '4967': 55.9,     # OGDH - Oxoglutarate dehydrogenase
        '6389': 60.4,     # SDHA - Succinate dehydrogenase complex subunit A

        # Cholesterol biosynthesis genes (낮은 발현)
        '3156': 25.8,     # HMGCR - HMG-CoA reductase
        '3157': 22.3,     # HMGCS1 - HMG-CoA synthase 1
        '4597': 18.9,     # MVD - Mevalonate decarboxylase

        # Fatty acid metabolism genes (중간 발현)
        '31': 45.2,       # ACACA - Acetyl-CoA carboxylase alpha
        '32': 42.8,       # ACACB - Acetyl-CoA carboxylase beta
        '2194': 48.7,     # FASN - Fatty acid synthase

        # Pentose phosphate pathway (중간 발현)
        '5236': 52.3,     # PGM1 - Phosphoglucomutase 1
        '5226': 49.1,     # PGD - Phosphogluconate dehydrogenase

        # Additional genes (다양한 발현 수준)
        '6714': 38.5,     # SRC - SRC proto-oncogene
        '7124': 72.9,     # TNF - Tumor necrosis factor
        '7157': 68.3,     # TP53 - Tumor protein p53
    }

    print(f"  Created expression data for {len(expression_data_entrez)} genes (Entrez IDs)")
    print(f"  {len(expression_data_entrez)}개 유전자의 발현 데이터 생성 완료 (Entrez ID)")

    # Show some examples
    print("\n  Sample data:")
    for i, (gene_id, expr_val) in enumerate(list(expression_data_entrez.items())[:5]):
        print(f"    Entrez ID {gene_id}: {expr_val:.1f}")
    print("    ...")

    return expression_data_entrez


def example1_direct_dict_to_datamap():
    """
    Example 1: Create data map directly from expression data dictionary

    예제 1: 발현 데이터 딕셔너리로부터 직접 데이터 맵 생성

    This is the simplest approach - just provide the expression data dict
    and the function handles everything else.

    가장 간단한 방법 - 발현 데이터 딕셔너리만 제공하면 나머지는 자동 처리됩니다.
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Direct dictionary to data map (simplest method)")
    print("예제 1: 딕셔너리에서 데이터 맵으로 직접 변환 (가장 간단한 방법)")
    print("="*80)

    # Load model
    model_wrapper = load_example_model()

    # Create expression data with Entrez IDs
    expression_data_entrez = create_example_expression_data_entrez()

    # Create data map with automatic ID conversion
    # ID 자동 변환을 통한 데이터 맵 생성
    print("\nCreating data map with automatic ID conversion...")
    print("자동 ID 변환을 통한 데이터 맵 생성 중...")

    data_map = create_data_map_from_dict(
        expression_data=expression_data_entrez,
        model_wrapper=model_wrapper,
        gene_id_type='entrez_id',  # Specify or let it auto-detect
        auto_convert=True,
        verbose=True
    )

    print(f"\nData map created successfully!")
    print(f"데이터 맵 생성 완료!")
    print(f"  Number of reactions scored: {len(data_map.get_scores())}")
    print(f"  반응 점수 개수: {len(data_map.get_scores())}")

    # Show sample reaction scores
    reaction_scores = data_map.get_scores()
    print("\n  Sample reaction scores:")
    for i, (rxn_id, score) in enumerate(list(reaction_scores.items())[:5]):
        print(f"    {rxn_id}: {score}")

    return model_wrapper, data_map


def example2_omics_container_approach():
    """
    Example 2: Using OmicsContainer with automatic preparation

    예제 2: OmicsContainer를 사용한 자동 준비

    This approach gives more control over the OmicsContainer while still
    benefiting from automatic ID conversion.

    OmicsContainer를 더 세밀하게 제어하면서도 자동 ID 변환의 이점을 활용합니다.
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: OmicsContainer with automatic ID preparation")
    print("예제 2: 자동 ID 준비를 사용하는 OmicsContainer")
    print("="*80)

    from troppo.omics import create_compatible_data_map

    # Load model
    model_wrapper = load_example_model()

    # Create expression data with Entrez IDs
    expression_data_entrez = create_example_expression_data_entrez()

    # Create OmicsContainer
    print("\nCreating OmicsContainer...")
    print("OmicsContainer 생성 중...")

    omics_container = OmicsContainer(
        omicstype='transcriptomics',
        condition='cancer_tissue',
        data=expression_data_entrez,
        nomenclature='entrez_id'  # Explicitly specify the nomenclature
    )

    print(f"  OmicsContainer created: {omics_container}")

    # Create compatible data map (with automatic ID conversion)
    # 호환 가능한 데이터 맵 생성 (자동 ID 변환 포함)
    print("\nCreating compatible data map...")
    print("호환 가능한 데이터 맵 생성 중...")

    data_map = create_compatible_data_map(
        omics_container=omics_container,
        model_wrapper=model_wrapper,
        auto_convert=True,
        verbose=True
    )

    print(f"\nData map created successfully!")
    print(f"  Number of reactions scored: {len(data_map.get_scores())}")

    return model_wrapper, data_map


def example3_benchmark_with_entrez_expression_data():
    """
    Example 3: Complete benchmark with Entrez ID expression data

    예제 3: Entrez ID 발현 데이터를 사용한 완전한 벤치마크

    This example shows the complete workflow: expression data with Entrez IDs,
    essentiality data with different IDs, and full benchmarking.

    이 예제는 완전한 워크플로우를 보여줍니다: Entrez ID 발현 데이터,
    다른 ID 형식의 필수성 데이터, 그리고 완전한 벤치마킹.
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Complete benchmark with mixed gene ID types")
    print("예제 3: 혼합 유전자 ID 타입을 사용한 완전한 벤치마크")
    print("="*80)

    # Load model
    model_wrapper = load_example_model()

    # Create expression data with Entrez IDs
    expression_data_entrez = create_example_expression_data_entrez()

    # Create data map
    print("\nCreating data map from Entrez ID expression data...")
    data_map = create_data_map_from_dict(
        expression_data=expression_data_entrez,
        model_wrapper=model_wrapper,
        gene_id_type='entrez_id',
        verbose=True
    )

    # Define essential and non-essential genes using gene symbols
    # (These will be auto-converted to match the model)
    # 유전자 심볼을 사용한 필수/비필수 유전자 정의
    # (모델에 맞게 자동 변환됩니다)
    print("\nDefining validation genes...")
    print("검증용 유전자 정의 중...")

    # Essential genes in Symbol format
    essential_genes_symbols = [
        'GAPDH', 'PGK1', 'PKM', 'ALDOA',
        'ENO1', 'PFKM', 'CS', 'IDH1'
    ]

    # Non-essential genes in Symbol format
    non_essential_genes_symbols = [
        'HMGCR', 'HMGCS1', 'MVD',
        'TP53', 'TNF', 'SRC'
    ]

    print(f"  Essential genes (symbols): {len(essential_genes_symbols)}")
    print(f"  Non-essential genes (symbols): {len(non_essential_genes_symbols)}")

    # Run benchmark with mixed ID types
    # 혼합 ID 타입으로 벤치마크 실행
    print("\nRunning benchmark...")
    print("벤치마크 실행 중...")
    print("(This may take a few minutes / 수 분 소요될 수 있습니다)")

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme', 'fastcore'],  # Limit to 2 methods for speed
        # Essential genes in Symbol format
        essential_genes=essential_genes_symbols,
        essential_genes_id_type='symbol',
        # Non-essential genes in Symbol format
        non_essential_genes=non_essential_genes_symbols,
        non_essential_genes_id_type='symbol',
        verbose=True
    )

    # Run with validation
    comparison = runner.run_benchmark(
        validate_essentiality=True,
        validate_growth=True
    )

    # Display results
    print("\n" + "="*80)
    print("BENCHMARK RESULTS / 벤치마크 결과")
    print("="*80)

    summary_df = comparison.get_summary_dataframe()
    print("\nPerformance Summary:")
    print(summary_df[['method_name', 'execution_time', 'peak_memory',
                      'num_reactions_selected', 'percentage_retained']])

    if comparison.results[list(comparison.results.keys())[0]].essentiality_accuracy is not None:
        print("\nEssentiality Validation:")
        print(summary_df[['method_name', 'essentiality_accuracy',
                         'essentiality_precision', 'essentiality_recall',
                         'essentiality_f1']])

    return comparison


def example4_multiple_id_types():
    """
    Example 4: Mixing different ID types in a single workflow

    예제 4: 하나의 워크플로우에서 여러 ID 타입 혼합 사용

    Demonstrates flexibility in handling:
    - Expression data: Entrez IDs
    - Essential genes: Ensembl IDs
    - Non-essential genes: Gene Symbols

    다양한 ID 타입 처리의 유연성 시연:
    - 발현 데이터: Entrez ID
    - 필수 유전자: Ensembl ID
    - 비필수 유전자: 유전자 심볼
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Multiple ID types in one workflow")
    print("예제 4: 하나의 워크플로우에서 여러 ID 타입 사용")
    print("="*80)

    # Load model
    model_wrapper = load_example_model()

    # Expression data with Entrez IDs
    print("\n1. Expression data: Entrez IDs")
    expression_data_entrez = {
        '2597': 100.5,  # GAPDH
        '5230': 95.2,   # PGK1
        '5315': 88.7,   # PKM
        '226': 92.3,    # ALDOA
    }

    # Essential genes with Ensembl IDs
    print("2. Essential genes: Ensembl IDs")
    essential_genes_ensembl = [
        'ENSG00000111640',  # GAPDH
        'ENSG00000102144',  # PGK1
        'ENSG00000067225',  # PKM
    ]

    # Non-essential genes with Symbols
    print("3. Non-essential genes: Gene Symbols")
    non_essential_genes_symbols = [
        'HMGCR',
        'TP53',
        'TNF'
    ]

    # Create data map from Entrez IDs
    print("\nCreating data map (auto-converting Entrez IDs)...")
    data_map = create_data_map_from_dict(
        expression_data=expression_data_entrez,
        model_wrapper=model_wrapper,
        gene_id_type='entrez_id',
        verbose=True
    )

    # Run benchmark (all ID types will be auto-converted)
    print("\nRunning benchmark with mixed ID types...")
    print("(All ID types will be automatically converted to match the model)")
    print("(모든 ID 타입이 모델에 맞게 자동 변환됩니다)")

    runner = BenchmarkRunner(
        model_wrapper=model_wrapper,
        data_map=data_map,
        methods=['gimme'],
        essential_genes=essential_genes_ensembl,
        essential_genes_id_type='ensembl_gene_id',
        non_essential_genes=non_essential_genes_symbols,
        non_essential_genes_id_type='symbol',
        verbose=True
    )

    comparison = runner.run_benchmark(validate_essentiality=True)

    print("\n✓ Success! All ID types were automatically converted and integrated.")
    print("✓ 성공! 모든 ID 타입이 자동으로 변환되고 통합되었습니다.")

    return comparison


def main():
    """
    Run all examples

    모든 예제 실행
    """
    print("\n" + "="*80)
    print("TROPPO: Expression Data with Different Gene ID Types")
    print("TROPPO: 다양한 유전자 ID 타입을 사용한 발현 데이터 처리")
    print("="*80)
    print("\nThis script demonstrates how to use expression data with Entrez IDs")
    print("(or any other gene ID nomenclature) with automatic conversion to match")
    print("the metabolic model.")
    print("\n이 스크립트는 Entrez ID(또는 다른 유전자 ID 명명법)를 사용하는")
    print("발현 데이터를 대사 모델에 맞게 자동 변환하여 사용하는 방법을 보여줍니다.")

    try:
        # Run examples
        print("\nRunning examples...")

        # Example 1: Direct approach
        model_wrapper, data_map = example1_direct_dict_to_datamap()

        # Example 2: OmicsContainer approach
        model_wrapper, data_map = example2_omics_container_approach()

        # Example 3: Complete benchmark
        comparison = example3_benchmark_with_entrez_expression_data()

        # Example 4: Multiple ID types
        comparison = example4_multiple_id_types()

        print("\n" + "="*80)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("모든 예제가 성공적으로 완료되었습니다!")
        print("="*80)

    except Exception as e:
        print(f"\nError running examples: {e}")
        print(f"예제 실행 중 오류 발생: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
