"""
Plotting functions for cfDNA fragmentomics visualization.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from fragmentomics.features.sizes import SizeDistribution

# Color scheme
COLORS = {
    "primary": "#2563eb",  # Blue
    "secondary": "#f97316",  # Orange
    "short": "#ef4444",  # Red (tumor-associated)
    "mono": "#22c55e",  # Green
    "di": "#8b5cf6",  # Purple
    "background": "#f8fafc",
    "grid": "#e2e8f0",
}


def plot_size_distribution(
    dist: SizeDistribution,
    ax: Axes | None = None,
    title: str | None = None,
    show_features: bool = True,
    color: str = COLORS["primary"],
    alpha: float = 0.7,
    figsize: tuple[float, float] = (10, 6),
) -> tuple[Figure, Axes]:
    """
    Plot fragment size distribution with optional feature annotations.

    Parameters
    ----------
    dist : SizeDistribution
        Size distribution object from analysis
    ax : Axes, optional
        Matplotlib axes to plot on
    title : str, optional
        Plot title
    show_features : bool, default True
        Whether to annotate peaks and regions
    color : str, default blue
        Line/fill color
    alpha : float, default 0.7
        Fill transparency
    figsize : tuple, default (10, 6)
        Figure size if creating new figure

    Returns
    -------
    tuple[Figure, Axes]
        Matplotlib figure and axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Normalize counts to density
    density = dist.counts / dist.counts.sum()

    # Plot distribution
    ax.fill_between(
        dist.bin_centers,
        density,
        alpha=alpha,
        color=color,
        label="Fragment sizes",
    )
    ax.plot(dist.bin_centers, density, color=color, linewidth=1)

    if show_features:
        # Highlight regions
        _add_region_highlights(ax, density, dist.bin_centers)

        # Mark peaks
        if dist.peak_mono:
            ax.axvline(
                dist.peak_mono,
                color=COLORS["mono"],
                linestyle="--",
                linewidth=1.5,
                label=f"Mono peak: {dist.peak_mono}bp",
            )

        if dist.peak_di:
            ax.axvline(
                dist.peak_di,
                color=COLORS["di"],
                linestyle="--",
                linewidth=1.5,
                label=f"Di peak: {dist.peak_di}bp",
            )

        # Add statistics box
        _add_stats_box(ax, dist)

    # Styling
    ax.set_xlabel("Fragment Size (bp)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(title or "Fragment Size Distribution", fontsize=14, fontweight="bold")
    ax.set_xlim(dist.bin_centers.min(), dist.bin_centers.max())
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3, color=COLORS["grid"])
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    return fig, ax


def _add_region_highlights(
    ax: Axes,
    density: NDArray,
    bin_centers: NDArray,
) -> None:
    """Add shaded regions for short, mono, and di fragments."""
    density.max() * 1.1

    # Short fragments (<150bp) - tumor-associated
    ax.axvspan(
        bin_centers.min(),
        150,
        alpha=0.1,
        color=COLORS["short"],
        label="Short (<150bp)",
    )

    # Mononucleosome region (140-180bp)
    ax.axvspan(140, 180, alpha=0.1, color=COLORS["mono"])

    # Dinucleosome region (280-360bp)
    ax.axvspan(280, 360, alpha=0.1, color=COLORS["di"])


def _add_stats_box(ax: Axes, dist: SizeDistribution) -> None:
    """Add a text box with key statistics."""
    stats_text = (
        f"n = {dist.n_fragments:,}\n"
        f"Median: {dist.median:.0f}bp\n"
        f"Mode: {dist.mode}bp\n"
        f"Short: {dist.ratio_short:.1%}\n"
        f"Mono: {dist.ratio_mono:.1%}"
    )

    props = dict(
        boxstyle="round", facecolor="white", alpha=0.8, edgecolor=COLORS["grid"]
    )
    ax.text(
        0.98,
        0.72,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=props,
        family="monospace",
    )


def plot_size_comparison(
    distributions: Sequence[SizeDistribution],
    labels: Sequence[str],
    ax: Axes | None = None,
    title: str | None = None,
    colors: Sequence[str] | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> tuple[Figure, Axes]:
    """
    Compare multiple fragment size distributions.

    Parameters
    ----------
    distributions : Sequence[SizeDistribution]
        List of distributions to compare
    labels : Sequence[str]
        Labels for each distribution
    ax : Axes, optional
        Matplotlib axes
    title : str, optional
        Plot title
    colors : Sequence[str], optional
        Colors for each distribution
    figsize : tuple, default (10, 6)
        Figure size

    Returns
    -------
    tuple[Figure, Axes]
    """
    if len(distributions) != len(labels):
        raise ValueError("Number of distributions must match number of labels")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    if colors is None:
        cmap = plt.cm.tab10
        colors = [cmap(i) for i in range(len(distributions))]

    for dist, label, color in zip(distributions, labels, colors):
        density = dist.counts / dist.counts.sum()
        ax.plot(
            dist.bin_centers,
            density,
            label=f"{label} (n={dist.n_fragments:,})",
            color=color,
            linewidth=1.5,
        )

    ax.set_xlabel("Fragment Size (bp)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(
        title or "Fragment Size Distribution Comparison", fontsize=14, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    ax.set_ylim(0, None)

    plt.tight_layout()
    return fig, ax


def plot_fragment_profile(
    sizes: NDArray,
    positions: NDArray,
    chrom: str,
    start: int,
    end: int,
    ax: Axes | None = None,
    window_size: int = 1000,
    figsize: tuple[float, float] = (12, 4),
) -> tuple[Figure, Axes]:
    """
    Plot fragment size profile across a genomic region.

    Shows how fragment sizes vary along a genomic region,
    useful for identifying nucleosome positioning.

    Parameters
    ----------
    sizes : NDArray
        Fragment sizes
    positions : NDArray
        Fragment midpoint positions
    chrom : str
        Chromosome name
    start : int
        Region start
    end : int
        Region end
    ax : Axes, optional
        Matplotlib axes
    window_size : int, default 1000
        Smoothing window size in bp
    figsize : tuple, default (12, 4)
        Figure size

    Returns
    -------
    tuple[Figure, Axes]
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Bin fragments by position
    n_bins = (end - start) // window_size
    bins = np.linspace(start, end, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    mean_sizes = []
    for i in range(n_bins):
        mask = (positions >= bins[i]) & (positions < bins[i + 1])
        if mask.sum() > 0:
            mean_sizes.append(np.mean(sizes[mask]))
        else:
            mean_sizes.append(np.nan)

    mean_sizes = np.array(mean_sizes)

    ax.plot(bin_centers, mean_sizes, color=COLORS["primary"], linewidth=1)
    ax.fill_between(
        bin_centers,
        mean_sizes,
        alpha=0.3,
        color=COLORS["primary"],
    )

    ax.set_xlabel(f"Position on {chrom}", fontsize=12)
    ax.set_ylabel("Mean Fragment Size (bp)", fontsize=12)
    ax.set_title(
        f"Fragment Size Profile: {chrom}:{start:,}-{end:,}",
        fontsize=14,
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig, ax


def save_figure(
    fig: Figure,
    path: str | Path,
    dpi: int = 150,
    **kwargs,
) -> None:
    """
    Save a figure with sensible defaults.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure
    path : str or Path
        Output path
    dpi : int, default 150
        Resolution in dots per inch
    **kwargs
        Additional arguments to savefig
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)


def plot_delfi_comparison(
    profiles: list,
    labels: list[str],
    colors: list[str] | None = None,
    ax: Axes | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (14, 5),
) -> tuple[Figure, Axes]:
    """
    Compare DELFI profiles from multiple samples.

    Overlays short/long ratios across the genome for healthy vs cancer
    or other comparisons.

    Parameters
    ----------
    profiles : list[DELFIProfile]
        DELFI profiles to compare
    labels : list[str]
        Sample labels
    colors : list[str], optional
        Colors for each profile
    ax : Axes, optional
        Matplotlib axes
    title : str, optional
        Plot title
    figsize : tuple, default (14, 5)
        Figure size

    Returns
    -------
    tuple[Figure, Axes]
    """
    if len(profiles) != len(labels):
        raise ValueError("Number of profiles must match number of labels")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    if colors is None:
        cmap = plt.cm.tab10
        colors = [cmap(i) for i in range(len(profiles))]

    for profile, label, color in zip(profiles, labels, colors):
        ratios = profile.to_ratio_vector()
        x = np.arange(len(ratios))
        ax.plot(x, ratios, label=label, color=color, linewidth=0.8, alpha=0.8)

    ax.set_xlabel("Genomic Bin", fontsize=12)
    ax.set_ylabel("Short/Long Ratio", fontsize=12)
    ax.set_title(
        title or "DELFI Fragmentation Profile Comparison",
        fontsize=14,
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    plt.tight_layout()
    return fig, ax
