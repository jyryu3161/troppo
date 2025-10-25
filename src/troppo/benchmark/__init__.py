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

__all__ = [
    'BenchmarkResult',
    'BenchmarkComparison',
    'BenchmarkRunner',
    'quick_benchmark'
]
