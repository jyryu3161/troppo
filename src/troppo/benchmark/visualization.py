"""
Visualization Tools for Benchmark Results

This module provides visualization utilities for comparing benchmark results.

벤치마크 결과 시각화 도구
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Any
import warnings


def plot_performance_comparison(
    comparison,
    metrics: Optional[List[str]] = None,
    figsize: tuple = (15, 10),
    save_path: Optional[str] = None
):
    """
    Plot comprehensive performance comparison

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    metrics : List[str], optional
        Metrics to plot (default: all available)
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure

    Examples
    --------
    >>> plot_performance_comparison(comparison, save_path='benchmark_results.png')
    """
    if metrics is None:
        metrics = ['execution_time', 'peak_memory', 'percentage_retained', 'biomass_flux']

    # Filter successful results
    successful_results = {
        name: result for name, result in comparison.results.items()
        if result.success
    }

    if not successful_results:
        print("No successful results to plot")
        return

    # Create subplots
    n_metrics = len(metrics)
    n_cols = 2
    n_rows = (n_metrics + 1) // 2

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_metrics > 1 else [axes]

    # Plot each metric
    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        # Collect data
        methods = []
        values = []

        for method_name, result in successful_results.items():
            if hasattr(result, metric):
                value = getattr(result, metric)
                if value is not None:
                    methods.append(method_name)
                    values.append(value)

        if not values:
            ax.text(0.5, 0.5, f'No data for {metric}',
                   ha='center', va='center', transform=ax.transAxes)
            continue

        # Create bar plot
        colors = sns.color_palette("husl", len(methods))
        bars = ax.bar(methods, values, color=colors)

        # Customize
        ax.set_title(metric.replace('_', ' ').title(), fontsize=12, fontweight='bold')
        ax.set_ylabel('Value')
        ax.tick_params(axis='x', rotation=45)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.2f}',
                   ha='center', va='bottom', fontsize=9)

        ax.grid(axis='y', alpha=0.3)

    # Remove empty subplots
    for idx in range(len(metrics), len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()


def plot_overlap_heatmap(
    comparison,
    figsize: tuple = (10, 8),
    save_path: Optional[str] = None
):
    """
    Plot reaction overlap heatmap between methods

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    overlap_matrix = comparison.get_overlap_matrix()

    plt.figure(figsize=figsize)
    sns.heatmap(
        overlap_matrix,
        annot=True,
        fmt='.3f',
        cmap='YlOrRd',
        square=True,
        cbar_kws={'label': 'Jaccard Similarity'},
        vmin=0,
        vmax=1
    )

    plt.title('Reaction Overlap Between Methods (Jaccard Similarity)',
             fontsize=14, fontweight='bold', pad=20)
    plt.xlabel('Method', fontsize=12)
    plt.ylabel('Method', fontsize=12)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved to {save_path}")

    plt.show()


def plot_pareto_front(
    comparison,
    x_metric: str = 'execution_time',
    y_metric: str = 'percentage_retained',
    figsize: tuple = (10, 6),
    save_path: Optional[str] = None
):
    """
    Plot Pareto front for two metrics

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    x_metric : str
        X-axis metric
    y_metric : str
        Y-axis metric
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    # Collect data
    data = []
    for method_name, result in comparison.results.items():
        if result.success:
            x_val = getattr(result, x_metric, None)
            y_val = getattr(result, y_metric, None)

            if x_val is not None and y_val is not None:
                data.append({
                    'method': method_name,
                    'x': x_val,
                    'y': y_val
                })

    if not data:
        print("No data to plot")
        return

    df = pd.DataFrame(data)

    # Create plot
    plt.figure(figsize=figsize)

    # Scatter plot
    colors = sns.color_palette("husl", len(df))
    for idx, row in df.iterrows():
        plt.scatter(row['x'], row['y'], s=200, c=[colors[idx]], alpha=0.7)
        plt.annotate(
            row['method'],
            (row['x'], row['y']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=10,
            fontweight='bold'
        )

    plt.xlabel(x_metric.replace('_', ' ').title(), fontsize=12)
    plt.ylabel(y_metric.replace('_', ' ').title(), fontsize=12)
    plt.title(f'{y_metric.replace("_", " ").title()} vs {x_metric.replace("_", " ").title()}',
             fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")

    plt.show()


def plot_radar_chart(
    comparison,
    metrics: Optional[List[str]] = None,
    normalize: bool = True,
    figsize: tuple = (10, 10),
    save_path: Optional[str] = None
):
    """
    Plot radar chart comparing methods across multiple metrics

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    metrics : List[str], optional
        Metrics to include
    normalize : bool
        Whether to normalize metrics to 0-1 range
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    if metrics is None:
        metrics = ['execution_time', 'peak_memory', 'percentage_retained',
                  'biomass_flux', 'task_completion_rate']

    # Filter successful results
    successful_results = {
        name: result for name, result in comparison.results.items()
        if result.success
    }

    if not successful_results:
        print("No successful results to plot")
        return

    # Collect data
    data = {}
    for metric in metrics:
        data[metric] = []
        for method_name in successful_results.keys():
            result = successful_results[method_name]
            value = getattr(result, metric, None)
            data[metric].append(value if value is not None else 0)

    # Normalize if requested
    if normalize:
        for metric in metrics:
            values = data[metric]
            if any(v != 0 for v in values):
                min_val = min(values)
                max_val = max(values)
                if max_val > min_val:
                    # Invert for time/memory (lower is better)
                    if metric in ['execution_time', 'peak_memory']:
                        data[metric] = [1 - (v - min_val) / (max_val - min_val) for v in values]
                    else:
                        data[metric] = [(v - min_val) / (max_val - min_val) for v in values]

    # Number of variables
    num_vars = len(metrics)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]  # Complete the circle

    # Create plot
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))

    # Plot each method
    colors = sns.color_palette("husl", len(successful_results))

    for idx, (method_name, result) in enumerate(successful_results.items()):
        values = [data[metric][idx] for metric in metrics]
        values += values[:1]  # Complete the circle

        ax.plot(angles, values, 'o-', linewidth=2, label=method_name, color=colors[idx])
        ax.fill(angles, values, alpha=0.15, color=colors[idx])

    # Customize
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([m.replace('_', ' ').title() for m in metrics])
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'])
    ax.grid(True)

    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    plt.title('Multi-Metric Comparison Radar Chart',
             fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Radar chart saved to {save_path}")

    plt.show()


def create_comparison_report(
    comparison,
    output_dir: str = '.',
    include_plots: bool = True
):
    """
    Create comprehensive comparison report

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    output_dir : str
        Output directory for report files
    include_plots : bool
        Whether to include visualizations
    """
    import os
    from datetime import datetime

    os.makedirs(output_dir, exist_ok=True)

    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save summary CSV
    summary_df = comparison.get_summary_dataframe()
    summary_path = os.path.join(output_dir, f'benchmark_summary_{timestamp}.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to {summary_path}")

    # Save overlap matrix
    overlap_df = comparison.get_overlap_matrix()
    overlap_path = os.path.join(output_dir, f'overlap_matrix_{timestamp}.csv')
    overlap_df.to_csv(overlap_path)
    print(f"Overlap matrix saved to {overlap_path}")

    # Save JSON
    json_path = os.path.join(output_dir, f'benchmark_results_{timestamp}.json')
    comparison.save_to_json(json_path)
    print(f"Full results saved to {json_path}")

    # Create plots
    if include_plots:
        try:
            # Performance comparison
            perf_path = os.path.join(output_dir, f'performance_comparison_{timestamp}.png')
            plot_performance_comparison(comparison, save_path=perf_path)

            # Overlap heatmap
            heatmap_path = os.path.join(output_dir, f'overlap_heatmap_{timestamp}.png')
            plot_overlap_heatmap(comparison, save_path=heatmap_path)

            # Pareto front
            pareto_path = os.path.join(output_dir, f'pareto_front_{timestamp}.png')
            plot_pareto_front(comparison, save_path=pareto_path)

            # Radar chart
            radar_path = os.path.join(output_dir, f'radar_chart_{timestamp}.png')
            plot_radar_chart(comparison, save_path=radar_path)

        except Exception as e:
            print(f"Warning: Some plots could not be created: {str(e)}")

    # Create markdown report
    report_path = os.path.join(output_dir, f'benchmark_report_{timestamp}.md')
    _create_markdown_report(comparison, report_path, timestamp)
    print(f"Markdown report saved to {report_path}")

    print(f"\n✓ Comparison report created in: {output_dir}")


def _create_markdown_report(comparison, filepath: str, timestamp: str):
    """Create markdown report"""
    with open(filepath, 'w') as f:
        f.write(f"# Omics Integration Methods Benchmark Report\n\n")
        f.write(f"**Generated:** {timestamp}\n\n")
        f.write(f"---\n\n")

        # Dataset info
        f.write("## Dataset Information\n\n")
        for key, value in comparison.dataset_info.items():
            f.write(f"- **{key.replace('_', ' ').title()}:** {value}\n")
        f.write("\n---\n\n")

        # Summary table
        f.write("## Performance Summary\n\n")
        summary_df = comparison.get_summary_dataframe()
        f.write(summary_df.to_markdown(index=False))
        f.write("\n\n---\n\n")

        # Rankings
        f.write("## Method Rankings\n\n")

        metrics = ['execution_time', 'peak_memory', 'percentage_retained']
        for metric in metrics:
            rankings = comparison.rank_methods(metric)
            if rankings:
                f.write(f"### {metric.replace('_', ' ').title()}\n\n")
                for rank, (method, value) in enumerate(rankings, 1):
                    f.write(f"{rank}. **{method}**: {value:.3f}\n")
                f.write("\n")

        f.write("---\n\n")

        # Common reactions
        common = comparison.get_common_reactions()
        f.write(f"## Common Reactions\n\n")
        f.write(f"**Number of reactions common to all methods:** {len(common)}\n\n")

        # Unique reactions
        f.write("## Unique Reactions by Method\n\n")
        for method_name in comparison.results.keys():
            unique = comparison.get_unique_reactions(method_name)
            f.write(f"- **{method_name}:** {len(unique)} unique reactions\n")

        f.write("\n---\n\n")

        # Detailed results
        f.write("## Detailed Results\n\n")
        for method_name, result in comparison.results.items():
            f.write(f"### {method_name.upper()}\n\n")

            if result.success:
                f.write(f"- **Execution Time:** {result.execution_time:.2f}s\n")
                f.write(f"- **Peak Memory:** {result.peak_memory:.2f} MB\n")
                f.write(f"- **Reactions Selected:** {result.num_reactions_selected}\n")
                f.write(f"- **Percentage Retained:** {result.percentage_retained:.2f}%\n")

                if result.biomass_flux is not None:
                    f.write(f"- **Biomass Flux:** {result.biomass_flux:.4f}\n")

                if result.task_completion_rate is not None:
                    f.write(f"- **Task Success Rate:** {result.task_completion_rate:.2f}%\n")

                if result.num_blocked_reactions is not None:
                    f.write(f"- **Blocked Reactions:** {result.num_blocked_reactions}\n")

                if result.network_consistency is not None:
                    f.write(f"- **Network Consistent:** {'Yes' if result.network_consistency else 'No'}\n")
            else:
                f.write(f"- **Status:** Failed\n")
                f.write(f"- **Error:** {result.error_message}\n")

            f.write("\n")

        f.write("---\n\n")
        f.write("*Report generated by Troppo Benchmark Framework*\n")
