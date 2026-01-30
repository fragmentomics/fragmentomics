"""
Coverage-based copy number analysis for cfDNA.

Implements binned coverage analysis with GC correction for
inferring copy number from cfDNA fragment data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


@dataclass
class CoverageBin:
    """A genomic bin with coverage information."""

    chrom: str
    start: int
    end: int
    raw_count: int
    gc_content: float = 0.0
    corrected_count: float = 0.0
    log2_ratio: float = 0.0

    @property
    def midpoint(self) -> int:
        """Bin midpoint."""
        return (self.start + self.end) // 2


@dataclass
class CoverageProfile:
    """
    Coverage profile across genomic bins.

    Attributes
    ----------
    bins : list[CoverageBin]
        Coverage bins
    bin_size : int
        Size of each bin in bp
    total_reads : int
        Total reads counted
    median_coverage : float
        Median raw coverage
    mad : float
        Median absolute deviation
    """

    bins: list[CoverageBin]
    bin_size: int
    total_reads: int = 0
    median_coverage: float = 0.0
    mad: float = 0.0
    gc_corrected: bool = False

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "n_bins": len(self.bins),
            "bin_size": self.bin_size,
            "total_reads": self.total_reads,
            "median_coverage": self.median_coverage,
            "mad": self.mad,
            "gc_corrected": self.gc_corrected,
        }

    def to_array(self) -> tuple[NDArray, NDArray, NDArray]:
        """
        Convert to numpy arrays.

        Returns
        -------
        positions : NDArray
            Bin midpoints
        raw_counts : NDArray
            Raw fragment counts
        log2_ratios : NDArray
            Log2 ratios (GC-corrected if available)
        """
        positions = np.array([b.midpoint for b in self.bins])
        raw_counts = np.array([b.raw_count for b in self.bins])
        log2_ratios = np.array([b.log2_ratio for b in self.bins])
        return positions, raw_counts, log2_ratios


class CoverageAnalyzer:
    """
    Analyze read coverage for copy number estimation.

    Computes binned coverage across the genome with optional
    GC bias correction.

    Parameters
    ----------
    bin_size : int, default 100000
        Size of genomic bins in bp (100kb default)
    min_mapq : int, default 30
        Minimum mapping quality
    min_fragment_size : int, default 50
        Minimum fragment size
    max_fragment_size : int, default 500
        Maximum fragment size

    Examples
    --------
    >>> analyzer = CoverageAnalyzer(bin_size=100000)
    >>> profile = analyzer.compute_coverage("sample.bam")
    >>> print(f"Median coverage: {profile.median_coverage:.1f}")
    """

    def __init__(
        self,
        bin_size: int = 100_000,
        min_mapq: int = 30,
        min_fragment_size: int = 50,
        max_fragment_size: int = 500,
    ):
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.min_fragment_size = min_fragment_size
        self.max_fragment_size = max_fragment_size

    def compute_coverage(
        self,
        bam_path: str | Path,
        chrom: str | None = None,
        reference: str | Path | None = None,
    ) -> CoverageProfile:
        """
        Compute binned coverage from BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file
        chrom : str, optional
            Specific chromosome to analyze (None = all)
        reference : str or Path, optional
            Reference FASTA for GC content calculation

        Returns
        -------
        CoverageProfile
            Computed coverage profile
        """
        bam_path = Path(bam_path)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Get chromosomes to process
            if chrom:
                chroms = [(chrom, bam.get_reference_length(chrom))]
            else:
                chroms = [
                    (bam.get_reference_name(i), bam.get_reference_length(bam.get_reference_name(i)))
                    for i in range(bam.nreferences)
                    if not bam.get_reference_name(i).startswith(("chrM", "MT", "chrUn", "random"))
                ]

            bins = []
            total_reads = 0

            for chrom_name, chrom_len in chroms:
                # Create bins for this chromosome
                for bin_start in range(0, chrom_len, self.bin_size):
                    bin_end = min(bin_start + self.bin_size, chrom_len)

                    # Count fragments in bin
                    count = 0
                    for read in bam.fetch(chrom_name, bin_start, bin_end):
                        if not self._is_valid_read(read):
                            continue

                        # Use fragment midpoint for binning
                        frag_start = read.reference_start
                        frag_end = frag_start + abs(read.template_length)
                        midpoint = (frag_start + frag_end) // 2

                        if bin_start <= midpoint < bin_end:
                            count += 1
                            total_reads += 1

                    bins.append(CoverageBin(
                        chrom=chrom_name,
                        start=bin_start,
                        end=bin_end,
                        raw_count=count,
                    ))

        # Compute median and MAD
        counts = np.array([b.raw_count for b in bins])
        median_cov = float(np.median(counts))
        mad = float(np.median(np.abs(counts - median_cov)))

        # Compute log2 ratios
        for b in bins:
            if median_cov > 0:
                ratio = (b.raw_count + 1) / (median_cov + 1)  # Pseudocount
                b.log2_ratio = float(np.log2(ratio))

        # Optional: GC correction
        if reference:
            self._correct_gc_bias(bins, reference)

        return CoverageProfile(
            bins=bins,
            bin_size=self.bin_size,
            total_reads=total_reads,
            median_coverage=median_cov,
            mad=mad,
            gc_corrected=reference is not None,
        )

    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read passes filters."""
        if read.is_unmapped:
            return False
        if not read.is_proper_pair:
            return False
        if read.is_duplicate:
            return False
        if not read.is_read1:
            return False
        if read.mapping_quality < self.min_mapq:
            return False

        frag_size = abs(read.template_length)
        if frag_size < self.min_fragment_size or frag_size > self.max_fragment_size:
            return False

        return True

    def _correct_gc_bias(
        self,
        bins: list[CoverageBin],
        reference: str | Path,
    ) -> None:
        """Apply GC bias correction to bins."""
        reference = Path(reference)

        try:
            with pysam.FastaFile(str(reference)) as fasta:
                for b in bins:
                    seq = fasta.fetch(b.chrom, b.start, b.end).upper()
                    gc_count = seq.count("G") + seq.count("C")
                    total = len(seq) - seq.count("N")
                    b.gc_content = gc_count / total if total > 0 else 0.5

            # Simple LOESS-style correction
            gc_values = np.array([b.gc_content for b in bins])
            counts = np.array([b.raw_count for b in bins])

            # Bin by GC content
            gc_bins = np.linspace(0.2, 0.8, 31)
            gc_medians = []
            for i in range(len(gc_bins) - 1):
                mask = (gc_values >= gc_bins[i]) & (gc_values < gc_bins[i + 1])
                if mask.sum() > 10:
                    gc_medians.append(np.median(counts[mask]))
                else:
                    gc_medians.append(np.nan)

            gc_medians = np.array(gc_medians)
            global_median = np.nanmedian(gc_medians)

            # Apply correction
            for b in bins:
                gc_idx = int((b.gc_content - 0.2) / 0.02)
                gc_idx = max(0, min(gc_idx, len(gc_medians) - 1))
                if not np.isnan(gc_medians[gc_idx]) and gc_medians[gc_idx] > 0:
                    correction = global_median / gc_medians[gc_idx]
                    b.corrected_count = b.raw_count * correction
                    b.log2_ratio = float(np.log2((b.corrected_count + 1) / (global_median + 1)))

        except Exception as e:
            logger.warning(f"GC correction failed: {e}")


def compute_coverage(
    bam_path: str | Path,
    bin_size: int = 100_000,
    chrom: str | None = None,
) -> CoverageProfile:
    """
    Convenience function to compute coverage.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    bin_size : int, default 100000
        Bin size in bp
    chrom : str, optional
        Specific chromosome

    Returns
    -------
    CoverageProfile
        Coverage profile
    """
    analyzer = CoverageAnalyzer(bin_size=bin_size)
    return analyzer.compute_coverage(bam_path, chrom=chrom)
