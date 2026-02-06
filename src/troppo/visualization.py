"""
Publication-quality visualization for tINIT reconstruction results.

All plots use 300 DPI, clean scientific styling suitable for journal submission.
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import Optional, List, TYPE_CHECKING

if TYPE_CHECKING:
    from troppo.results import ReconstructionResults

# Publication style defaults
STYLE_DEFAULTS = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
}

# Nature-style color palette
COLORS = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F',
          '#EDC948', '#B07AA1', '#FF9DA7', '#9C755F', '#BAB0AC']


def _apply_style():
    """Apply publication style settings."""
    plt.rcParams.update(STYLE_DEFAULTS)


def _save_if_needed(fig, save_path):
    """Save figure if path is provided."""
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')


def plot_model_sizes(results: 'ReconstructionResults', ax=None, **kwargs) -> Figure:
    """Bar plot of active reaction counts per sample.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the plot.
    """
    _apply_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (max(6, len(results.sample_names) * 0.5), 4)))
    else:
        fig = ax.get_figure()

    counts = results.n_reactions_per_sample.sort_values(ascending=False)
    bars = ax.bar(range(len(counts)), counts.values, color=COLORS[0], edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha='right')
    ax.set_ylabel('Active Reactions')
    ax.set_title('Model Size per Sample')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add mean line
    mean_val = counts.mean()
    ax.axhline(y=mean_val, color=COLORS[1], linestyle='--', linewidth=1, label=f'Mean: {mean_val:.0f}')
    ax.legend()

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_reaction_overlap_heatmap(results: 'ReconstructionResults', ax=None, **kwargs) -> Figure:
    """Jaccard similarity heatmap with hierarchical clustering.

    Computes pairwise Jaccard similarity between all sample reaction sets
    and displays the result as a clustered heatmap. Hierarchical clustering
    is applied when scipy is available; otherwise the original sample order
    is preserved.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the heatmap.
    """
    _apply_style()

    jac = results.jaccard_similarity()
    n = len(jac)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (max(6, n * 0.6), max(5, n * 0.5))))
    else:
        fig = ax.get_figure()

    # Cluster if scipy available
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        dist = 1 - jac.values
        np.fill_diagonal(dist, 0)
        condensed = squareform(dist)
        Z = linkage(condensed, method='average')
        order = leaves_list(Z)
        jac = jac.iloc[order, order]
    except ImportError:
        pass

    im = ax.imshow(jac.values, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(jac.columns, rotation=45, ha='right')
    ax.set_yticklabels(jac.index)
    ax.set_title('Jaccard Similarity')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Jaccard Index')

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_subsystem_coverage(results: 'ReconstructionResults', ax=None, top_n: int = 20, **kwargs) -> Figure:
    """Subsystem-level reaction activation heatmap.

    Displays the fraction of reactions active in each metabolic subsystem
    across all samples. Subsystems are ranked by variance to highlight
    the most differentially active pathways.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    top_n : int, optional
        Number of most variable subsystems to display. Default is 20.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the heatmap.
    """
    _apply_style()

    coverage = results.subsystem_coverage()
    # Select top_n most variable subsystems
    variance = coverage.var()
    top_subsystems = variance.nlargest(top_n).index
    coverage = coverage[top_subsystems]

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (max(8, len(results.sample_names) * 0.5), max(6, top_n * 0.3))))
    else:
        fig = ax.get_figure()

    im = ax.imshow(coverage.T.values, cmap='viridis', aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(len(coverage)))
    ax.set_yticks(range(len(top_subsystems)))
    ax.set_xticklabels(coverage.index, rotation=45, ha='right')
    ax.set_yticklabels(top_subsystems, fontsize=8)
    ax.set_title(f'Subsystem Coverage (Top {top_n})')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Fraction Active')

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_biomass_distribution(results: 'ReconstructionResults', ax=None, **kwargs) -> Figure:
    """Biomass flux distribution as a violin plot with overlaid data points.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the violin plot.
    """
    _apply_style()

    fluxes = results.biomass_fluxes()

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (5, 4)))
    else:
        fig = ax.get_figure()

    parts = ax.violinplot([fluxes.values], positions=[0], showmeans=True, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor(COLORS[0])
        pc.set_alpha(0.7)
    parts['cmeans'].set_color(COLORS[1])
    parts['cmedians'].set_color(COLORS[2])

    ax.scatter(np.zeros(len(fluxes)), fluxes.values, color=COLORS[0], alpha=0.6, s=20, zorder=3)
    ax.set_ylabel('Biomass Flux')
    ax.set_title('Biomass Flux Distribution')
    ax.set_xticks([0])
    ax.set_xticklabels(['Samples'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_pca_models(results: 'ReconstructionResults', ax=None, labels: list = None, **kwargs) -> Figure:
    """PCA scatter plot based on reaction composition.

    Projects the binary reaction-presence matrix into two principal
    components and plots each sample as a labeled point. Requires
    scikit-learn for PCA computation.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    labels : list, optional
        Custom labels for each sample point. Defaults to sample names.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the PCA scatter plot.

    Raises
    ------
    ImportError
        If scikit-learn is not installed.
    """
    _apply_style()

    from sklearn.decomposition import PCA

    matrix = results.reaction_matrix.values.astype(float)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(matrix)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (6, 5)))
    else:
        fig = ax.get_figure()

    scatter = ax.scatter(
        coords[:, 0], coords[:, 1],
        c=COLORS[0], s=50,
        edgecolors='white', linewidth=0.5, zorder=3,
    )

    if labels is None:
        labels = results.sample_names
    for i, label in enumerate(labels):
        ax.annotate(
            label, (coords[i, 0], coords[i, 1]),
            fontsize=7, ha='left', va='bottom',
            xytext=(4, 4), textcoords='offset points',
        )

    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)')
    ax.set_title('PCA of Reaction Composition')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_essentiality_metrics(results: 'ReconstructionResults', validation_result: pd.DataFrame,
                              ax=None, **kwargs) -> Figure:
    """Gene essentiality validation bar chart.

    Displays mean values of accuracy, precision, recall, F1-score, and
    Matthews correlation coefficient from gene essentiality validation.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    validation_result : pd.DataFrame
        DataFrame with essentiality validation metrics. Expected columns
        include 'accuracy', 'precision', 'recall', 'f1_score', and
        'matthews_correlation'.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the bar chart.
    """
    _apply_style()

    metrics = ['accuracy', 'precision', 'recall', 'f1_score', 'matthews_correlation']
    available = [m for m in metrics if m in validation_result.columns]

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (8, 4)))
    else:
        fig = ax.get_figure()

    mean_metrics = validation_result[available].mean()
    x = range(len(available))
    bars = ax.bar(x, mean_metrics.values, color=COLORS[:len(available)], edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('_', ' ').title() for m in available], rotation=30, ha='right')
    ax.set_ylabel('Score')
    ax.set_ylim(0, 1.05)
    ax.set_title('Gene Essentiality Validation')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for bar, val in zip(bars, mean_metrics.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
            f'{val:.3f}', ha='center', va='bottom', fontsize=8,
        )

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_yield_comparison(results: 'ReconstructionResults', yield_result: pd.DataFrame,
                          ax=None, **kwargs) -> Figure:
    """Theoretical yield comparison bar chart.

    Plots mean theoretical yield grouped by carbon source. Falls back
    to a placeholder message when no yield data is available.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    yield_result : pd.DataFrame
        DataFrame with yield analysis results. Expected to contain
        'theoretical_yield' and 'carbon_source' columns.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the bar chart.
    """
    _apply_style()

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (7, 4)))
    else:
        fig = ax.get_figure()

    if 'theoretical_yield' in yield_result.columns:
        data = yield_result.groupby('carbon_source')['theoretical_yield'].mean()
        bars = ax.bar(range(len(data)), data.values, color=COLORS[0], edgecolor='white', linewidth=0.5)
        ax.set_xticks(range(len(data)))
        ax.set_xticklabels(data.index, rotation=30, ha='right')
        ax.set_ylabel('Theoretical Yield')
        ax.set_title('Mean Theoretical Yield by Carbon Source')
    else:
        ax.text(0.5, 0.5, 'No yield data available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Theoretical Yield')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_reaction_frequency(results: 'ReconstructionResults', ax=None, top_n: int = 30, **kwargs) -> Figure:
    """Reaction activation frequency across samples.

    Displays a horizontal bar chart showing the percentage of samples
    in which each reaction is active, sorted by frequency.

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new Figure and Axes are created.
    top_n : int, optional
        Number of most frequent reactions to display. Default is 30.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height).
        save_path : str, optional
            File path to save the figure.

    Returns
    -------
    Figure
        The matplotlib Figure containing the horizontal bar chart.
    """
    _apply_style()

    freq = results.reaction_matrix.sum(axis=0).sort_values(ascending=False)
    freq_pct = (freq / results.n_samples * 100)

    if top_n:
        freq_pct = freq_pct.head(top_n)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, max(4, top_n * 0.2))))
    else:
        fig = ax.get_figure()

    y_pos = range(len(freq_pct))
    ax.barh(y_pos, freq_pct.values, color=COLORS[0], edgecolor='white', linewidth=0.3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(freq_pct.index, fontsize=7)
    ax.set_xlabel('Samples with Reaction Active (%)')
    ax.set_title(f'Reaction Frequency (Top {top_n})')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    _save_if_needed(fig, kwargs.get('save_path'))
    return fig


def plot_overview(results: 'ReconstructionResults', save_dir: str = None,
                  format: str = 'pdf', **kwargs) -> Figure:
    """Generate comprehensive 2x3 overview figure.

    Creates a multi-panel figure combining the most informative plots
    into a single publication-ready overview:

      - Row 1: Model sizes, Jaccard similarity heatmap, reaction frequency
      - Row 2: Subsystem coverage, PCA scatter, summary statistics

    Parameters
    ----------
    results : ReconstructionResults
        Reconstruction results object containing reaction data.
    save_dir : str, optional
        Directory to save the overview figure. The directory is created
        if it does not exist.
    format : str, optional
        Output format for the saved figure (e.g. 'pdf', 'png', 'svg').
        Default is 'pdf'.
    **kwargs
        figsize : tuple, optional
            Figure size in inches (width, height). Default is (16, 10).

    Returns
    -------
    Figure
        The matplotlib Figure containing the six-panel overview.
    """
    _apply_style()

    fig, axes = plt.subplots(2, 3, figsize=kwargs.get('figsize', (16, 10)))
    fig.suptitle('tINIT Reconstruction Overview', fontsize=14, fontweight='bold')

    # Row 1
    plot_model_sizes(results, ax=axes[0, 0])
    plot_reaction_overlap_heatmap(results, ax=axes[0, 1])
    plot_reaction_frequency(results, ax=axes[0, 2], top_n=15)

    # Row 2
    plot_subsystem_coverage(results, ax=axes[1, 0], top_n=10)

    # PCA - may fail if sklearn not installed
    try:
        plot_pca_models(results, ax=axes[1, 1])
    except ImportError:
        axes[1, 1].text(
            0.5, 0.5, 'Install scikit-learn\nfor PCA plot',
            ha='center', va='center',
            transform=axes[1, 1].transAxes, fontsize=10, color='gray',
        )
        axes[1, 1].set_title('PCA of Reaction Composition')

    # Summary stats text
    ax_text = axes[1, 2]
    ax_text.axis('off')
    summary = results.summary()
    text_lines = [
        f"Samples: {results.n_samples}",
        f"Total reactions: {len(results.reaction_matrix.columns)}",
        f"Mean active: {summary['n_active_reactions'].mean():.0f}",
        f"Std active: {summary['n_active_reactions'].std():.0f}",
        f"Common reactions: {len(results.get_common_reactions())}",
        f"Mean retained: {summary['pct_retained'].mean():.1f}%",
    ]
    ax_text.text(
        0.1, 0.8, '\n'.join(text_lines),
        transform=ax_text.transAxes,
        fontsize=11, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.3),
    )
    ax_text.set_title('Summary Statistics')

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save_dir:
        import os
        os.makedirs(save_dir, exist_ok=True)
        fig.savefig(os.path.join(save_dir, f'overview.{format}'), dpi=300, bbox_inches='tight')

    return fig
