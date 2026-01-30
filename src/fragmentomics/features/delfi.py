"""
DELFI-style fragmentation profile analysis.

Implements the approach from Cristiano et al. 2019 (Nature):
"Genome-wide cell-free DNA fragmentation in patients with cancer"
DOI: 10.1038/s41586-019-1272-6

DELFI (DNA Evaluation of Fragments for Early Interception) analyzes
genome-wide fragmentation patterns in 5Mb windows, computing the ratio
of short (100-150bp) to long (151-220bp) fragments.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# DELFI default parameters from Cristiano et al.
DEFAULT_BIN_SIZE = 5_000_000  # 5 Mb windows
SHORT_MIN = 100
SHORT_MAX = 150
LONG_MIN = 151
LONG_MAX = 220


@dataclass
class DELFIBin:
    """A genomic bin with DELFI fragmentation features."""

    chrom: str
    start: int
    end: int
    short_count: int = 0
    long_count: int = 0
    total_count: int = 0
    gc_content: float = 0.0

    @property
    def ratio(self) -> float:
        """Short/long fragment ratio (the key DELFI feature)."""
        if self.long_count == 0:
            return 0.0
        return self.short_count / self.long_count

    @property
    def coverage(self) -> float:
        """Normalized coverage (fragments per Mb)."""
        bin_size_mb = (self.end - self.start) / 1_000_000
        return self.total_count / bin_size_mb if bin_size_mb > 0 else 0.0


@dataclass
class DELFIProfile:
    """
    DELFI fragmentation profile for a sample.

    Contains genome-wide fragmentation features in 5Mb windows,
    along with summary statistics and ML-ready feature vectors.

    Attributes
    ----------
    bins : list[DELFIBin]
        All genomic bins with fragment counts
    bin_size : int
        Size of each bin in bp
    total_fragments : int
        Total fragments analyzed
    genome_wide_ratio : float
        Overall short/long ratio across genome
    profile_correlation : float
        Correlation with healthy reference (if computed)
    """

    bins: list[DELFIBin]
    bin_size: int = DEFAULT_BIN_SIZE
    total_fragments: int = 0
    total_short: int = 0
    total_long: int = 0
    genome_wide_ratio: float = 0.0
    profile_correlation: float = 0.0
    sample_name: str = ""

    def to_feature_vector(self) -> NDArray[np.float64]:
        """
        Convert to feature vector for ML.

        Returns array of [ratio_1, ratio_2, ..., ratio_n, coverage_1, ..., coverage_n]
        """
        ratios = np.array([b.ratio for b in self.bins])
        coverages = np.array([b.coverage for b in self.bins])
        return np.concatenate([ratios, coverages])

    def to_ratio_vector(self) -> NDArray[np.float64]:
        """Return just the short/long ratios (primary DELFI feature)."""
        return np.array([b.ratio for b in self.bins])

    def to_coverage_vector(self) -> NDArray[np.float64]:
        """Return coverage values."""
        return np.array([b.coverage for b in self.bins])

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "sample_name": self.sample_name,
            "n_bins": len(self.bins),
            "bin_size": self.bin_size,
            "total_fragments": self.total_fragments,
            "total_short": self.total_short,
            "total_long": self.total_long,
            "genome_wide_ratio": self.genome_wide_ratio,
            "profile_correlation": self.profile_correlation,
        }

    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            f"DELFI Profile: {self.sample_name or 'unnamed'}\n"
            f"{'=' * 40}\n"
            f"Bins: {len(self.bins)} x {self.bin_size // 1_000_000}Mb\n"
            f"Total fragments: {self.total_fragments:,}\n"
            f"  Short (100-150bp): {self.total_short:,}\n"
            f"  Long (151-220bp): {self.total_long:,}\n"
            f"Genome-wide S/L ratio: {self.genome_wide_ratio:.3f}\n"
            f"Profile correlation: {self.profile_correlation:.3f}\n"
        )


class DELFIAnalyzer:
    """
    DELFI-style fragmentation analysis.

    Computes genome-wide short/long fragment ratios in 5Mb bins,
    following the approach from Cristiano et al. 2019.

    Parameters
    ----------
    bin_size : int, default 5_000_000
        Size of genomic bins (5 Mb recommended)
    short_range : tuple, default (100, 150)
        Size range for "short" fragments
    long_range : tuple, default (151, 220)
        Size range for "long" fragments
    min_mapq : int, default 30
        Minimum mapping quality
    exclude_chroms : list, optional
        Chromosomes to exclude (default: chrM, chrY, unplaced)

    Examples
    --------
    >>> analyzer = DELFIAnalyzer()
    >>> profile = analyzer.analyze("sample.bam")
    >>> features = profile.to_feature_vector()
    >>> print(f"Genome-wide ratio: {profile.genome_wide_ratio:.3f}")
    """

    def __init__(
        self,
        bin_size: int = DEFAULT_BIN_SIZE,
        short_range: tuple[int, int] = (SHORT_MIN, SHORT_MAX),
        long_range: tuple[int, int] = (LONG_MIN, LONG_MAX),
        min_mapq: int = 30,
        exclude_chroms: list[str] | None = None,
    ):
        self.bin_size = bin_size
        self.short_min, self.short_max = short_range
        self.long_min, self.long_max = long_range
        self.min_mapq = min_mapq
        self.exclude_chroms = exclude_chroms or [
            "chrM", "MT", "chrY", "chrUn", "random", "hap"
        ]

    def analyze(
        self,
        bam_path: str | Path,
        reference: str | Path | None = None,
        sample_name: str | None = None,
    ) -> DELFIProfile:
        """
        Compute DELFI fragmentation profile from BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file
        reference : str or Path, optional
            Reference FASTA for GC content
        sample_name : str, optional
            Sample identifier

        Returns
        -------
        DELFIProfile
            Computed fragmentation profile
        """
        bam_path = Path(bam_path)
        sample_name = sample_name or bam_path.stem

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Get chromosomes to process
            chroms = []
            for i in range(bam.nreferences):
                name = bam.get_reference_name(i)
                if not self._should_exclude(name):
                    chroms.append((name, bam.get_reference_length(name)))

            bins = []
            total_short = 0
            total_long = 0
            total_frags = 0

            for chrom_name, chrom_len in chroms:
                for bin_start in range(0, chrom_len, self.bin_size):
                    bin_end = min(bin_start + self.bin_size, chrom_len)

                    short_count = 0
                    long_count = 0
                    total_count = 0

                    for read in bam.fetch(chrom_name, bin_start, bin_end):
                        if not self._is_valid_read(read):
                            continue

                        frag_size = abs(read.template_length)

                        # Use fragment midpoint for binning
                        frag_start = read.reference_start
                        frag_end = frag_start + frag_size
                        midpoint = (frag_start + frag_end) // 2

                        if not (bin_start <= midpoint < bin_end):
                            continue

                        total_count += 1
                        total_frags += 1

                        if self.short_min <= frag_size <= self.short_max:
                            short_count += 1
                            total_short += 1
                        elif self.long_min <= frag_size <= self.long_max:
                            long_count += 1
                            total_long += 1

                    bins.append(DELFIBin(
                        chrom=chrom_name,
                        start=bin_start,
                        end=bin_end,
                        short_count=short_count,
                        long_count=long_count,
                        total_count=total_count,
                    ))

        # Compute genome-wide ratio
        genome_ratio = total_short / total_long if total_long > 0 else 0.0

        return DELFIProfile(
            bins=bins,
            bin_size=self.bin_size,
            total_fragments=total_frags,
            total_short=total_short,
            total_long=total_long,
            genome_wide_ratio=genome_ratio,
            sample_name=sample_name,
        )

    def _should_exclude(self, chrom: str) -> bool:
        """Check if chromosome should be excluded."""
        for pattern in self.exclude_chroms:
            if pattern.lower() in chrom.lower():
                return True
        return False

    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read passes filters."""
        if read.is_unmapped:
            return False
        if not read.is_proper_pair:
            return False
        if read.is_duplicate:
            return False
        if not read.is_read1:  # Avoid double-counting pairs
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        return True

    @staticmethod
    def compute_correlation(
        profile: DELFIProfile,
        reference_ratios: NDArray[np.float64],
    ) -> float:
        """
        Compute correlation with a reference profile.

        Parameters
        ----------
        profile : DELFIProfile
            Sample profile to compare
        reference_ratios : NDArray
            Reference ratio vector (e.g., median of healthy samples)

        Returns
        -------
        float
            Pearson correlation coefficient
        """
        sample_ratios = profile.to_ratio_vector()

        if len(sample_ratios) != len(reference_ratios):
            raise ValueError("Profile and reference must have same number of bins")

        # Handle edge cases
        if np.std(sample_ratios) == 0 or np.std(reference_ratios) == 0:
            return 0.0

        correlation = np.corrcoef(sample_ratios, reference_ratios)[0, 1]
        return float(correlation)


def compute_delfi_profile(
    bam_path: str | Path,
    bin_size: int = DEFAULT_BIN_SIZE,
) -> DELFIProfile:
    """
    Convenience function to compute DELFI profile.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    bin_size : int, default 5_000_000
        Bin size in bp

    Returns
    -------
    DELFIProfile
        Computed profile
    """
    analyzer = DELFIAnalyzer(bin_size=bin_size)
    return analyzer.analyze(bam_path)


def create_healthy_reference(profiles: list[DELFIProfile]) -> NDArray[np.float64]:
    """
    Create reference profile from healthy samples.

    Parameters
    ----------
    profiles : list[DELFIProfile]
        Profiles from healthy individuals

    Returns
    -------
    NDArray
        Median ratio at each genomic position
    """
    if not profiles:
        raise ValueError("Need at least one profile")

    # Stack all ratio vectors
    ratios = np.array([p.to_ratio_vector() for p in profiles])

    # Return median across samples
    return np.median(ratios, axis=0)
