"""
Extended Visualization Tools for Benchmark Results

Additional visualizations for gene essentiality and theoretical yields.

Gene essentiality와 이론적 수율을 위한 추가 시각화 도구
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Dict
import warnings


def plot_essentiality_comparison(
    comparison,
    metrics: Optional[List[str]] = None,
    figsize: tuple = (12, 8),
    save_path: Optional[str] = None
):
    """
    Plot gene essentiality metrics comparison across methods

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    metrics : List[str], optional
        Metrics to plot (default: accuracy, precision, recall, f1, mcc)
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    if metrics is None:
        metrics = ['essentiality_accuracy', 'essentiality_precision',
                  'essentiality_recall', 'essentiality_f1', 'essentiality_mcc']

    # Collect data
    data = []
    for method_name, result in comparison.results.items():
        if result.success:
            row = {'method': method_name}
            for metric in metrics:
                value = getattr(result, metric, None)
                if value is not None:
                    row[metric] = value
            if len(row) > 1:  # Has at least one metric
                data.append(row)

    if not data:
        print("No essentiality data to plot")
        return

    df = pd.DataFrame(data)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for grouped bar chart
    x = np.arange(len(df))
    width = 0.15
    metric_labels = {
        'essentiality_accuracy': 'Accuracy',
        'essentiality_precision': 'Precision',
        'essentiality_recall': 'Recall',
        'essentiality_f1': 'F1-Score',
        'essentiality_mcc': 'MCC'
    }

    colors = sns.color_palette("Set2", len(metrics))

    for idx, metric in enumerate(metrics):
        if metric in df.columns:
            offset = width * (idx - len(metrics)/2 + 0.5)
            ax.bar(x + offset, df[metric], width, label=metric_labels.get(metric, metric),
                  color=colors[idx], alpha=0.8)

    # Customize
    ax.set_xlabel('Method', fontsize=12, fontweight='bold')
    ax.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax.set_title('Gene Essentiality Prediction Performance', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(df['method'], rotation=45, ha='right')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim([0, 1.05])

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Essentiality comparison saved to {save_path}")

    plt.show()


def plot_confusion_matrices(
    comparison,
    figsize: tuple = (15, 5),
    save_path: Optional[str] = None
):
    """
    Plot confusion matrices for gene essentiality predictions

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    # Collect essentiality results
    results_with_essentiality = {}
    for method_name, result in comparison.results.items():
        if result.success and result.essentiality_result is not None:
            results_with_essentiality[method_name] = result.essentiality_result

    if not results_with_essentiality:
        print("No essentiality results to plot")
        return

    n_methods = len(results_with_essentiality)
    n_cols = min(3, n_methods)
    n_rows = (n_methods + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_methods == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, (method_name, ess_result) in enumerate(results_with_essentiality.items()):
        ax = axes[idx]

        # Create confusion matrix
        cm = np.array([
            [len(ess_result.true_positives), len(ess_result.false_negatives)],
            [len(ess_result.false_positives), len(ess_result.true_negatives)]
        ])

        # Plot heatmap
        sns.heatmap(
            cm,
            annot=True,
            fmt='d',
            cmap='Blues',
            ax=ax,
            cbar=False,
            xticklabels=['Essential', 'Non-Essential'],
            yticklabels=['Pred. Essential', 'Pred. Non-Essential']
        )

        ax.set_title(f'{method_name.upper()}\nAcc: {ess_result.accuracy:.3f}, '
                    f'F1: {ess_result.f1_score:.3f}',
                    fontweight='bold')
        ax.set_xlabel('Actual', fontweight='bold')
        ax.set_ylabel('Predicted', fontweight='bold')

    # Remove empty subplots
    for idx in range(n_methods, len(axes)):
        fig.delaxes(axes[idx])

    plt.suptitle('Gene Essentiality Confusion Matrices', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Confusion matrices saved to {save_path}")

    plt.show()


def plot_yield_comparison(
    comparison,
    metric: str = 'theoretical_yield',
    figsize: tuple = (14, 6),
    save_path: Optional[str] = None
):
    """
    Plot theoretical yields comparison across methods and carbon sources

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    metric : str
        Metric to plot ('theoretical_yield', 'biomass_yield', etc.)
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    # Collect yield data
    data = []
    for method_name, result in comparison.results.items():
        if result.success and result.yield_results:
            for yield_result in result.yield_results:
                data.append({
                    'method': method_name,
                    'carbon_source': yield_result.carbon_source,
                    'condition': 'Aerobic' if yield_result.aerobic else 'Anaerobic',
                    metric: getattr(yield_result, metric),
                    'feasible': yield_result.feasible
                })

    if not data:
        print("No yield data to plot")
        return

    df = pd.DataFrame(data)

    # Create subplots for aerobic and anaerobic
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Aerobic plot
    aerobic_data = df[df['condition'] == 'Aerobic']
    if not aerobic_data.empty:
        pivot_aerobic = aerobic_data.pivot_table(
            index='carbon_source',
            columns='method',
            values=metric,
            aggfunc='first'
        )

        pivot_aerobic.plot(kind='bar', ax=ax1, width=0.8, alpha=0.8)
        ax1.set_title('Aerobic Conditions (O₂ Available)', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Carbon Source', fontsize=10, fontweight='bold')
        ax1.set_ylabel(metric.replace('_', ' ').title(), fontsize=10, fontweight='bold')
        ax1.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(axis='y', alpha=0.3)
        ax1.tick_params(axis='x', rotation=45)

    # Anaerobic plot
    anaerobic_data = df[df['condition'] == 'Anaerobic']
    if not anaerobic_data.empty:
        pivot_anaerobic = anaerobic_data.pivot_table(
            index='carbon_source',
            columns='method',
            values=metric,
            aggfunc='first'
        )

        pivot_anaerobic.plot(kind='bar', ax=ax2, width=0.8, alpha=0.8)
        ax2.set_title('Anaerobic Conditions (No O₂)', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Carbon Source', fontsize=10, fontweight='bold')
        ax2.set_ylabel(metric.replace('_', ' ').title(), fontsize=10, fontweight='bold')
        ax2.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(axis='y', alpha=0.3)
        ax2.tick_params(axis='x', rotation=45)

    plt.suptitle(f'Theoretical Yields Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Yield comparison saved to {save_path}")

    plt.show()


def plot_yield_heatmap(
    comparison,
    condition: str = 'aerobic',
    metric: str = 'theoretical_yield',
    figsize: tuple = (12, 8),
    save_path: Optional[str] = None
):
    """
    Plot heatmap of yields across methods and carbon sources

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    condition : str
        'aerobic' or 'anaerobic'
    metric : str
        Metric to plot
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    # Collect data
    data = []
    for method_name, result in comparison.results.items():
        if result.success and result.yield_results:
            for yield_result in result.yield_results:
                is_aerobic = yield_result.aerobic
                if (condition.lower() == 'aerobic' and is_aerobic) or \
                   (condition.lower() == 'anaerobic' and not is_aerobic):
                    data.append({
                        'method': method_name,
                        'carbon_source': yield_result.carbon_source,
                        metric: getattr(yield_result, metric) if yield_result.feasible else 0
                    })

    if not data:
        print("No yield data to plot")
        return

    df = pd.DataFrame(data)

    # Create pivot table
    pivot = df.pivot_table(
        index='carbon_source',
        columns='method',
        values=metric,
        aggfunc='first'
    )

    # Create heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(
        pivot,
        annot=True,
        fmt='.3f',
        cmap='YlGnBu',
        cbar_kws={'label': metric.replace('_', ' ').title()},
        linewidths=0.5
    )

    plt.title(f'{condition.title()} Yields Heatmap', fontsize=14, fontweight='bold', pad=20)
    plt.xlabel('Method', fontsize=12, fontweight='bold')
    plt.ylabel('Carbon Source', fontsize=12, fontweight='bold')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Yield heatmap saved to {save_path}")

    plt.show()


def plot_comprehensive_validation(
    comparison,
    figsize: tuple = (16, 12),
    save_path: Optional[str] = None
):
    """
    Create comprehensive validation plot with all metrics

    Parameters
    ----------
    comparison : BenchmarkComparison
        Comparison results
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # 1. Performance metrics (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    _plot_performance_bars(comparison, ax1, ['execution_time', 'peak_memory'])

    # 2. Model size metrics (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    _plot_performance_bars(comparison, ax2, ['percentage_retained', 'biomass_flux'])

    # 3. Essentiality metrics (middle row)
    ax3 = fig.add_subplot(gs[1, :])
    _plot_essentiality_bars(comparison, ax3)

    # 4. Aerobic yields (bottom left)
    ax4 = fig.add_subplot(gs[2, 0])
    _plot_yield_bars(comparison, ax4, aerobic=True)

    # 5. Anaerobic yields (bottom right)
    ax5 = fig.add_subplot(gs[2, 1])
    _plot_yield_bars(comparison, ax5, aerobic=False)

    plt.suptitle('Comprehensive Benchmark Validation Results',
                fontsize=16, fontweight='bold', y=0.995)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Comprehensive validation plot saved to {save_path}")

    plt.show()


def _plot_performance_bars(comparison, ax, metrics):
    """Helper function for performance bars"""
    data = []
    for method_name, result in comparison.results.items():
        if result.success:
            row = {'method': method_name}
            for metric in metrics:
                row[metric] = getattr(result, metric, None)
            data.append(row)

    if not data:
        return

    df = pd.DataFrame(data)
    x = np.arange(len(df))
    width = 0.35

    for idx, metric in enumerate(metrics):
        if metric in df.columns:
            offset = width * (idx - 0.5)
            ax.bar(x + offset, df[metric], width, label=metric.replace('_', ' ').title(),
                  alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(df['method'], rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.set_title('Performance Metrics', fontweight='bold')


def _plot_essentiality_bars(comparison, ax):
    """Helper function for essentiality bars"""
    metrics = ['essentiality_accuracy', 'essentiality_precision',
              'essentiality_recall', 'essentiality_f1']

    data = []
    for method_name, result in comparison.results.items():
        if result.success and result.essentiality_accuracy is not None:
            row = {'method': method_name}
            for metric in metrics:
                row[metric] = getattr(result, metric, None)
            data.append(row)

    if not data:
        ax.text(0.5, 0.5, 'No essentiality data', ha='center', va='center',
               transform=ax.transAxes)
        return

    df = pd.DataFrame(data)
    x = np.arange(len(df))
    width = 0.2

    colors = sns.color_palette("Set2", len(metrics))
    for idx, metric in enumerate(metrics):
        if metric in df.columns:
            offset = width * (idx - len(metrics)/2 + 0.5)
            label = metric.replace('essentiality_', '').title()
            ax.bar(x + offset, df[metric], width, label=label,
                  color=colors[idx], alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(df['method'], rotation=45, ha='right')
    ax.legend(loc='upper left')
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim([0, 1.05])
    ax.set_title('Gene Essentiality Prediction', fontweight='bold')
    ax.set_ylabel('Score')


def _plot_yield_bars(comparison, ax, aerobic=True):
    """Helper function for yield bars"""
    data = []
    for method_name, result in comparison.results.items():
        if result.success and result.yield_results:
            yields = [
                yr.theoretical_yield for yr in result.yield_results
                if yr.aerobic == aerobic and yr.feasible
            ]
            if yields:
                data.append({
                    'method': method_name,
                    'avg_yield': np.mean(yields),
                    'std_yield': np.std(yields)
                })

    if not data:
        ax.text(0.5, 0.5, 'No yield data', ha='center', va='center',
               transform=ax.transAxes)
        return

    df = pd.DataFrame(data)
    x = np.arange(len(df))

    ax.bar(x, df['avg_yield'], yerr=df['std_yield'], capsize=5, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(df['method'], rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)
    condition = 'Aerobic' if aerobic else 'Anaerobic'
    ax.set_title(f'{condition} Yields', fontweight='bold')
    ax.set_ylabel('Avg Theoretical Yield')
