"""
GC bias correction for cfDNA fragmentomics.

GC content affects sequencing efficiency, creating biases in fragment
coverage and size distributions. This module implements correction methods.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam
from numpy.typing import NDArray
from scipy import interpolate, stats

logger = logging.getLogger(__name__)


@dataclass
class GCBiasProfile:
    """
    GC bias profile for a sample.

    Attributes
    ----------
    gc_bins : NDArray
        GC content bin centers (0-1)
    observed_counts : NDArray
        Observed fragment counts per GC bin
    expected_counts : NDArray
        Expected counts (uniform or reference)
    correction_factors : NDArray
        Multiplicative correction factors
    r_squared : float
        R² of GC vs coverage relationship
    gc_dropout : float
        Fraction of GC bins with < 10% expected coverage
    """

    gc_bins: NDArray
    observed_counts: NDArray
    expected_counts: NDArray
    correction_factors: NDArray
    r_squared: float = 0.0
    gc_dropout: float = 0.0
    mean_gc: float = 0.0
    std_gc: float = 0.0

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "r_squared": self.r_squared,
            "gc_dropout": self.gc_dropout,
            "mean_gc": self.mean_gc,
            "std_gc": self.std_gc,
            "n_bins": len(self.gc_bins),
        }

    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            f"GC Bias Profile\n"
            f"===============\n"
            f"Mean GC: {self.mean_gc:.1%}\n"
            f"Std GC: {self.std_gc:.1%}\n"
            f"R² (GC vs coverage): {self.r_squared:.3f}\n"
            f"GC dropout: {self.gc_dropout:.1%}\n"
        )


class GCCorrector:
    """
    GC bias correction for cfDNA data.

    Implements LOESS-based GC correction similar to methods used in
    copy number analysis and cfDNA studies.

    Parameters
    ----------
    n_bins : int, default 101
        Number of GC bins (0%, 1%, ..., 100%)
    smoothing_span : float, default 0.3
        LOESS smoothing span
    min_count : int, default 10
        Minimum fragments per bin for reliable correction

    Examples
    --------
    >>> corrector = GCCorrector()
    >>> profile = corrector.compute_bias("sample.bam", "reference.fa")
    >>> corrected_counts = corrector.correct(raw_counts, gc_values)
    """

    def __init__(
        self,
        n_bins: int = 101,
        smoothing_span: float = 0.3,
        min_count: int = 10,
    ):
        self.n_bins = n_bins
        self.smoothing_span = smoothing_span
        self.min_count = min_count
        self.gc_bins = np.linspace(0, 1, n_bins)

    def compute_gc_content(
        self,
        sequence: str,
    ) -> float:
        """Compute GC content of a sequence."""
        if not sequence:
            return 0.0
        gc = sum(1 for b in sequence.upper() if b in "GC")
        return gc / len(sequence)

    def compute_bias_from_fragments(
        self,
        gc_values: NDArray[np.float64],
        weights: NDArray | None = None,
    ) -> GCBiasProfile:
        """
        Compute GC bias profile from fragment GC values.

        Parameters
        ----------
        gc_values : NDArray
            GC content for each fragment (0-1)
        weights : NDArray, optional
            Weights for each fragment

        Returns
        -------
        GCBiasProfile
            Computed bias profile
        """
        if weights is None:
            weights = np.ones_like(gc_values)

        # Bin fragments by GC content
        bin_edges = np.linspace(0, 1, self.n_bins + 1)
        observed, _ = np.histogram(gc_values, bins=bin_edges, weights=weights)

        # Expected under uniform distribution
        expected = np.full(self.n_bins, len(gc_values) / self.n_bins)

        # Compute correction factors
        # Avoid division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            correction = np.where(observed > self.min_count, expected / observed, 1.0)

        # Smooth correction factors with LOESS-like approach
        correction_smooth = self._smooth_corrections(correction, observed)

        # Compute statistics
        mean_gc = float(np.average(gc_values, weights=weights))
        std_gc = float(np.sqrt(np.average((gc_values - mean_gc) ** 2, weights=weights)))

        # R² of GC vs coverage
        valid = observed > self.min_count
        if valid.sum() > 2:
            slope, intercept, r, p, se = stats.linregress(
                self.gc_bins[valid], observed[valid]
            )
            r_squared = r**2
        else:
            r_squared = 0.0

        # GC dropout (bins with very low coverage)
        gc_dropout = np.mean(observed < expected * 0.1)

        return GCBiasProfile(
            gc_bins=self.gc_bins,
            observed_counts=observed,
            expected_counts=expected,
            correction_factors=correction_smooth,
            r_squared=r_squared,
            gc_dropout=gc_dropout,
            mean_gc=mean_gc,
            std_gc=std_gc,
        )

    def _smooth_corrections(
        self,
        corrections: NDArray,
        counts: NDArray,
    ) -> NDArray:
        """Apply smoothing to correction factors."""
        # Use weighted moving average
        valid = counts > self.min_count

        if valid.sum() < 5:
            return corrections

        # Interpolate over invalid regions
        x_valid = self.gc_bins[valid]
        y_valid = corrections[valid]

        # Create interpolation function
        f = interpolate.interp1d(
            x_valid,
            y_valid,
            kind="linear",
            bounds_error=False,
            fill_value=(y_valid[0], y_valid[-1]),
        )

        smoothed = f(self.gc_bins)

        # Apply additional smoothing
        window = max(3, int(self.n_bins * self.smoothing_span / 10))
        if window % 2 == 0:
            window += 1

        kernel = np.ones(window) / window
        smoothed = np.convolve(smoothed, kernel, mode="same")

        return smoothed

    def correct_values(
        self,
        values: NDArray,
        gc_values: NDArray,
        profile: GCBiasProfile,
    ) -> NDArray:
        """
        Apply GC correction to values.

        Parameters
        ----------
        values : NDArray
            Values to correct (e.g., counts, coverage)
        gc_values : NDArray
            GC content for each value
        profile : GCBiasProfile
            Pre-computed bias profile

        Returns
        -------
        NDArray
            Corrected values
        """
        # Map GC values to correction factors
        bin_indices = np.digitize(gc_values, np.linspace(0, 1, self.n_bins + 1)) - 1
        bin_indices = np.clip(bin_indices, 0, self.n_bins - 1)

        corrections = profile.correction_factors[bin_indices]

        return values * corrections

    def compute_fragment_gc(
        self,
        bam_path: str | Path,
        reference_path: str | Path,
        region: str | None = None,
        max_fragments: int = 100000,
    ) -> NDArray[np.float64]:
        """
        Compute GC content for fragments in a BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file
        reference_path : str or Path
            Path to reference FASTA
        region : str, optional
            Genomic region to analyze
        max_fragments : int, default 100000
            Maximum fragments to sample

        Returns
        -------
        NDArray
            GC content for each fragment
        """
        bam_path = Path(bam_path)
        reference_path = Path(reference_path)

        gc_values = []
        ref = pysam.FastaFile(str(reference_path))

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                iterator = bam.fetch(region=region) if region else bam.fetch()

                for read in iterator:
                    if not self._is_valid_read(read):
                        continue

                    # Get fragment sequence
                    try:
                        chrom = read.reference_name
                        start = read.reference_start
                        end = start + abs(read.template_length)

                        if end - start > 1000 or end - start < 50:
                            continue

                        seq = ref.fetch(chrom, start, end)
                        gc = self.compute_gc_content(seq)
                        gc_values.append(gc)

                    except Exception:
                        continue

                    if len(gc_values) >= max_fragments:
                        break
        finally:
            ref.close()

        return np.array(gc_values, dtype=np.float64)

    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read is valid for GC analysis."""
        if read.is_unmapped:
            return False
        if not read.is_proper_pair:
            return False
        if read.is_duplicate:
            return False
        if not read.is_read1:
            return False
        if read.mapping_quality < 30:
            return False
        return True


def compute_gc_bias(
    bam_path: str | Path,
    reference_path: str | Path,
    region: str | None = None,
    max_fragments: int = 100000,
) -> GCBiasProfile:
    """
    Convenience function to compute GC bias profile.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    reference_path : str or Path
        Path to reference FASTA
    region : str, optional
        Genomic region
    max_fragments : int, default 100000
        Max fragments to sample

    Returns
    -------
    GCBiasProfile
        Computed bias profile
    """
    corrector = GCCorrector()
    gc_values = corrector.compute_fragment_gc(
        bam_path, reference_path, region, max_fragments
    )
    return corrector.compute_bias_from_fragments(gc_values)
