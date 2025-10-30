"""
Troppo Benchmark Module

Tools for comparing and evaluating omics integration methods.

벤치마크 및 방법론 비교 도구
"""

from .framework import (
    BenchmarkResult,
    BenchmarkComparison,
    BenchmarkRunner,
    quick_benchmark
)

from .validation import (
    GeneEssentialityValidator,
    GeneEssentialityResult,
    TheoreticalYieldCalculator,
    TheoreticalYieldResult,
    compare_yields_across_methods,
    create_essentiality_confusion_matrix
)

from .gene_id_utils import (
    detect_id_type,
    convert_gene_ids,
    standardize_gene_ids,
    match_gene_lists,
    FlexibleGeneIDHandler,
    prepare_validation_genes,
    normalize_id_type
)

__all__ = [
    'BenchmarkResult',
    'BenchmarkComparison',
    'BenchmarkRunner',
    'quick_benchmark',
    'GeneEssentialityValidator',
    'GeneEssentialityResult',
    'TheoreticalYieldCalculator',
    'TheoreticalYieldResult',
    'compare_yields_across_methods',
    'create_essentiality_confusion_matrix',
    'detect_id_type',
    'convert_gene_ids',
    'standardize_gene_ids',
    'match_gene_lists',
    'FlexibleGeneIDHandler',
    'prepare_validation_genes',
    'normalize_id_type'
]
