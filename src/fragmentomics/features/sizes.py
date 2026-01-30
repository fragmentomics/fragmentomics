"""
Fragment size distribution analysis.

Analyzes cfDNA fragment size distributions to extract features useful for
cancer detection and tissue-of-origin inference.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy import signal

logger = logging.getLogger(__name__)


@dataclass
class SizeDistribution:
    """
    Fragment size distribution with computed statistics.

    Attributes
    ----------
    sizes : NDArray
        Raw fragment sizes
    counts : NDArray
        Histogram counts
    bin_edges : NDArray
        Histogram bin edges
    bin_centers : NDArray
        Histogram bin centers

    Statistics
    ----------
    mean : float
    median : float
    mode : int
    std : float
    min_size : int
    max_size : int
    n_fragments : int

    Fragmentomics Features
    ----------------------
    ratio_short : float
        Fraction of fragments < 150bp (tumor-enriched)
    ratio_mono : float
        Fraction in mononucleosome range (140-180bp)
    ratio_di : float
        Fraction in dinucleosome range (280-360bp)
    peak_mono : int
        Position of mononucleosome peak
    peak_di : int
        Position of dinucleosome peak (if present)
    amplitude_ratio : float
        Ratio of di/mono peak amplitudes
    periodicity_10bp : float
        Strength of 10bp periodicity (nucleosome positioning)
    """

    # Raw data
    sizes: NDArray = field(repr=False)
    counts: NDArray = field(repr=False)
    bin_edges: NDArray = field(repr=False)
    bin_centers: NDArray = field(repr=False)

    # Basic statistics
    mean: float = 0.0
    median: float = 0.0
    mode: int = 0
    std: float = 0.0
    min_size: int = 0
    max_size: int = 0
    n_fragments: int = 0

    # Fragmentomics features
    ratio_short: float = 0.0  # < 150bp
    ratio_mono: float = 0.0  # 140-180bp
    ratio_di: float = 0.0  # 280-360bp
    peak_mono: int = 0
    peak_di: int | None = None
    amplitude_ratio: float | None = None
    periodicity_10bp: float = 0.0

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "n_fragments": self.n_fragments,
            "mean": self.mean,
            "median": self.median,
            "mode": self.mode,
            "std": self.std,
            "min_size": self.min_size,
            "max_size": self.max_size,
            "ratio_short": self.ratio_short,
            "ratio_mono": self.ratio_mono,
            "ratio_di": self.ratio_di,
            "peak_mono": self.peak_mono,
            "peak_di": self.peak_di,
            "amplitude_ratio": self.amplitude_ratio,
            "periodicity_10bp": self.periodicity_10bp,
        }

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = [
            "Fragment Size Distribution Summary",
            "=" * 40,
            f"Total fragments: {self.n_fragments:,}",
            f"Size range: {self.min_size} - {self.max_size} bp",
            "",
            "Basic Statistics:",
            f"  Mean:   {self.mean:.1f} bp",
            f"  Median: {self.median:.1f} bp",
            f"  Mode:   {self.mode} bp",
            f"  Std:    {self.std:.1f} bp",
            "",
            "Fragmentomics Features:",
            f"  Short fragments (<150bp): {self.ratio_short:.1%}",
            f"  Mononucleosome (140-180bp): {self.ratio_mono:.1%}",
            f"  Dinucleosome (280-360bp): {self.ratio_di:.1%}",
            f"  Mononucleosome peak: {self.peak_mono} bp",
        ]

        if self.peak_di:
            lines.append(f"  Dinucleosome peak: {self.peak_di} bp")
        if self.amplitude_ratio:
            lines.append(f"  Di/Mono amplitude ratio: {self.amplitude_ratio:.3f}")

        lines.append(f"  10bp periodicity score: {self.periodicity_10bp:.3f}")

        return "\n".join(lines)


class FragmentSizeAnalyzer:
    """
    Analyzer for cfDNA fragment size distributions.

    Computes standard fragmentomics metrics including size ratios,
    peak positions, and nucleosome periodicity.

    Parameters
    ----------
    bin_size : int, default 1
        Histogram bin size in bp
    min_size : int, default 50
        Minimum size for analysis
    max_size : int, default 500
        Maximum size for analysis
    smooth_window : int, default 5
        Window size for peak smoothing

    Examples
    --------
    >>> from fragmentomics.io import read_fragments
    >>> sizes = read_fragments("sample.bam")
    >>> analyzer = FragmentSizeAnalyzer()
    >>> dist = analyzer.analyze(sizes)
    >>> print(dist.summary())
    """

    def __init__(
        self,
        bin_size: int = 1,
        min_size: int = 50,
        max_size: int = 500,
        smooth_window: int = 5,
    ):
        self.bin_size = bin_size
        self.min_size = min_size
        self.max_size = max_size
        self.smooth_window = smooth_window

    def analyze(self, sizes: NDArray[np.int32]) -> SizeDistribution:
        """
        Analyze a fragment size distribution.

        Parameters
        ----------
        sizes : NDArray
            Array of fragment sizes in bp

        Returns
        -------
        SizeDistribution
            Distribution object with computed features
        """
        if len(sizes) == 0:
            raise ValueError("No fragments provided for analysis")

        # Filter to analysis range
        mask = (sizes >= self.min_size) & (sizes <= self.max_size)
        sizes_filtered = sizes[mask]

        if len(sizes_filtered) == 0:
            raise ValueError(
                f"No fragments in range [{self.min_size}, {self.max_size}]"
            )

        # Compute histogram
        bins = np.arange(self.min_size, self.max_size + self.bin_size, self.bin_size)
        counts, bin_edges = np.histogram(sizes_filtered, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Normalize to density
        counts_norm = counts / counts.sum()

        # Basic statistics
        mean = float(np.mean(sizes_filtered))
        median = float(np.median(sizes_filtered))
        mode = int(bin_centers[np.argmax(counts)])
        std = float(np.std(sizes_filtered))

        # Size ratios
        ratio_short = np.sum(sizes_filtered < 150) / len(sizes_filtered)
        ratio_mono = np.sum((sizes_filtered >= 140) & (sizes_filtered <= 180)) / len(
            sizes_filtered
        )
        ratio_di = np.sum((sizes_filtered >= 280) & (sizes_filtered <= 360)) / len(
            sizes_filtered
        )

        # Find peaks
        peak_mono, peak_di, amplitude_ratio = self._find_peaks(counts_norm, bin_centers)

        # Compute 10bp periodicity
        periodicity_10bp = self._compute_periodicity(counts_norm, bin_centers)

        return SizeDistribution(
            sizes=sizes_filtered,
            counts=counts,
            bin_edges=bin_edges,
            bin_centers=bin_centers,
            mean=mean,
            median=median,
            mode=mode,
            std=std,
            min_size=int(sizes_filtered.min()),
            max_size=int(sizes_filtered.max()),
            n_fragments=len(sizes_filtered),
            ratio_short=float(ratio_short),
            ratio_mono=float(ratio_mono),
            ratio_di=float(ratio_di),
            peak_mono=peak_mono,
            peak_di=peak_di,
            amplitude_ratio=amplitude_ratio,
            periodicity_10bp=periodicity_10bp,
        )

    def _find_peaks(
        self,
        counts: NDArray,
        bin_centers: NDArray,
    ) -> tuple[int, int | None, float | None]:
        """Find mono- and di-nucleosome peaks."""
        # Smooth the distribution
        if self.smooth_window > 1:
            kernel = np.ones(self.smooth_window) / self.smooth_window
            counts_smooth = np.convolve(counts, kernel, mode="same")
        else:
            counts_smooth = counts

        # Find peaks with minimum prominence
        peaks, properties = signal.find_peaks(
            counts_smooth,
            prominence=0.001,
            distance=20,
        )

        if len(peaks) == 0:
            # Fallback: use maximum in mononucleosome range
            mono_mask = (bin_centers >= 140) & (bin_centers <= 200)
            if mono_mask.any():
                mono_idx = np.argmax(counts_smooth[mono_mask])
                peak_mono = int(bin_centers[mono_mask][mono_idx])
            else:
                peak_mono = int(bin_centers[np.argmax(counts_smooth)])
            return peak_mono, None, None

        peak_positions = bin_centers[peaks].astype(int)
        peak_heights = counts_smooth[peaks]

        # Find mononucleosome peak (largest peak in 140-200 range)
        mono_mask = (peak_positions >= 140) & (peak_positions <= 200)
        if mono_mask.any():
            mono_peaks = peaks[mono_mask]
            mono_heights = peak_heights[mono_mask]
            peak_mono = int(bin_centers[mono_peaks[np.argmax(mono_heights)]])
            mono_height = mono_heights.max()
        else:
            # Use global maximum as fallback
            peak_mono = int(bin_centers[peaks[np.argmax(peak_heights)]])
            mono_height = peak_heights.max()

        # Find dinucleosome peak (largest peak in 280-400 range)
        di_mask = (peak_positions >= 280) & (peak_positions <= 400)
        if di_mask.any():
            di_peaks = peaks[di_mask]
            di_heights = peak_heights[di_mask]
            peak_di = int(bin_centers[di_peaks[np.argmax(di_heights)]])
            di_height = di_heights.max()
            amplitude_ratio = di_height / mono_height if mono_height > 0 else None
        else:
            peak_di = None
            amplitude_ratio = None

        return peak_mono, peak_di, amplitude_ratio

    def _compute_periodicity(
        self,
        counts: NDArray,
        bin_centers: NDArray,
    ) -> float:
        """
        Compute 10bp periodicity score using autocorrelation.

        The 10bp periodicity reflects nucleosome positioning and is
        characteristic of cfDNA fragmentation patterns.
        """
        # Focus on mononucleosome region (100-200bp)
        mask = (bin_centers >= 100) & (bin_centers <= 200)
        if mask.sum() < 20:
            return 0.0

        region = counts[mask]

        # Compute autocorrelation
        autocorr = np.correlate(region, region, mode="full")
        autocorr = autocorr[len(autocorr) // 2 :]  # Take positive lags
        autocorr = autocorr / autocorr[0]  # Normalize

        # Look for peak around lag 10 (10bp periodicity)
        # With 1bp bins, lag 10 = 10bp
        if len(autocorr) > 15:
            # Average autocorrelation around 10bp lag
            periodicity = np.mean(autocorr[8:12])
        else:
            periodicity = 0.0

        return float(periodicity)


def analyze_sizes(
    sizes: NDArray[np.int32],
    bin_size: int = 1,
    min_size: int = 50,
    max_size: int = 500,
) -> SizeDistribution:
    """
    Convenience function to analyze fragment sizes.

    Parameters
    ----------
    sizes : NDArray
        Array of fragment sizes
    bin_size : int, default 1
        Histogram bin size
    min_size : int, default 50
        Minimum size for analysis
    max_size : int, default 500
        Maximum size for analysis

    Returns
    -------
    SizeDistribution
        Distribution with computed features

    Examples
    --------
    >>> from fragmentomics.io import read_fragments
    >>> from fragmentomics.features import analyze_sizes
    >>> sizes = read_fragments("sample.bam")
    >>> dist = analyze_sizes(sizes)
    >>> print(f"Short fragment ratio: {dist.ratio_short:.1%}")
    """
    analyzer = FragmentSizeAnalyzer(
        bin_size=bin_size,
        min_size=min_size,
        max_size=max_size,
    )
    return analyzer.analyze(sizes)
