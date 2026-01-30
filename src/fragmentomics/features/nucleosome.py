"""
Nucleosome positioning analysis using cfDNA fragmentation patterns.

Implements Windowed Protection Score (WPS) and related metrics for
inferring nucleosome positions from cfDNA fragment data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy import signal, ndimage
import pysam

logger = logging.getLogger(__name__)

# Nucleosome constants
NUCLEOSOME_SIZE = 147  # bp wrapped around histone
LINKER_SIZE = 20  # typical linker between nucleosomes
NUCLEOSOME_REPEAT = NUCLEOSOME_SIZE + LINKER_SIZE  # ~167bp


@dataclass
class WPSProfile:
    """
    Windowed Protection Score profile for a genomic region.
    
    WPS measures the difference between fragments spanning a position
    (protective) and fragments with endpoints at that position (cutting).
    High WPS indicates nucleosome protection; low WPS indicates accessible DNA.
    
    Attributes
    ----------
    positions : NDArray
        Genomic positions
    wps : NDArray
        WPS values at each position
    smoothed_wps : NDArray
        Smoothed WPS for peak calling
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    """
    
    positions: NDArray
    wps: NDArray
    smoothed_wps: NDArray
    chrom: str
    start: int
    end: int
    
    # Derived metrics
    mean_wps: float = 0.0
    std_wps: float = 0.0
    peak_positions: NDArray = field(default_factory=lambda: np.array([]))
    trough_positions: NDArray = field(default_factory=lambda: np.array([]))
    periodicity: float = 0.0
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "mean_wps": self.mean_wps,
            "std_wps": self.std_wps,
            "n_peaks": len(self.peak_positions),
            "n_troughs": len(self.trough_positions),
            "periodicity": self.periodicity,
        }
    
    def summary(self) -> str:
        """Return human-readable summary."""
        return (
            f"WPS Profile: {self.chrom}:{self.start:,}-{self.end:,}\n"
            f"===============================================\n"
            f"Mean WPS: {self.mean_wps:.2f}\n"
            f"Std WPS: {self.std_wps:.2f}\n"
            f"Nucleosome peaks: {len(self.peak_positions)}\n"
            f"Periodicity (bp): {self.periodicity:.1f}\n"
        )


class NucleosomeAnalyzer:
    """
    Analyzer for nucleosome positioning using cfDNA fragments.
    
    Computes Windowed Protection Score (WPS) and identifies nucleosome
    positions based on fragment coverage patterns.
    
    Parameters
    ----------
    window_size : int, default 120
        Window size for WPS calculation
    min_fragment_size : int, default 120
        Minimum fragment size to consider
    max_fragment_size : int, default 180
        Maximum fragment size (mononucleosome range)
    smoothing_bandwidth : int, default 15
        Bandwidth for Gaussian smoothing
    min_mapq : int, default 30
        Minimum mapping quality
        
    Examples
    --------
    >>> analyzer = NucleosomeAnalyzer()
    >>> profile = analyzer.compute_wps("sample.bam", "chr1", 1000000, 1001000)
    >>> print(f"Found {len(profile.peak_positions)} nucleosomes")
    """
    
    def __init__(
        self,
        window_size: int = 120,
        min_fragment_size: int = 120,
        max_fragment_size: int = 180,
        smoothing_bandwidth: int = 15,
        min_mapq: int = 30,
    ):
        self.window_size = window_size
        self.min_fragment_size = min_fragment_size
        self.max_fragment_size = max_fragment_size
        self.smoothing_bandwidth = smoothing_bandwidth
        self.min_mapq = min_mapq
    
    def compute_wps(
        self,
        bam_path: Union[str, Path],
        chrom: str,
        start: int,
        end: int,
    ) -> WPSProfile:
        """
        Compute Windowed Protection Score for a genomic region.
        
        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file
        chrom : str
            Chromosome name
        start : int
            Start position
        end : int
            End position
            
        Returns
        -------
        WPSProfile
            Computed WPS profile
        """
        bam_path = Path(bam_path)
        region_size = end - start
        
        # Initialize arrays for fragment tracking
        protection = np.zeros(region_size, dtype=np.int32)  # Spanning fragments
        cutting = np.zeros(region_size, dtype=np.int32)     # Endpoint fragments
        
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(chrom, start, end):
                if not self._is_valid_read(read):
                    continue
                
                frag_start = read.reference_start
                frag_end = frag_start + abs(read.template_length)
                frag_size = frag_end - frag_start
                
                # Filter by size (mononucleosome)
                if frag_size < self.min_fragment_size or frag_size > self.max_fragment_size:
                    continue
                
                # Convert to local coordinates
                local_start = max(0, frag_start - start)
                local_end = min(region_size, frag_end - start)
                
                # Protection: positions covered by fragment body
                if local_end > local_start:
                    protection[local_start:local_end] += 1
                
                # Cutting: fragment endpoints
                if 0 <= frag_start - start < region_size:
                    cutting[frag_start - start] += 1
                if 0 <= frag_end - start < region_size:
                    cutting[frag_end - start] += 1
        
        # Compute WPS: protection - cutting (windowed)
        half_window = self.window_size // 2
        wps = np.zeros(region_size, dtype=np.float64)
        
        for i in range(half_window, region_size - half_window):
            window_protection = protection[i - half_window:i + half_window].sum()
            window_cutting = cutting[i - half_window:i + half_window].sum()
            wps[i] = window_protection - window_cutting
        
        # Smooth WPS
        smoothed_wps = ndimage.gaussian_filter1d(wps, sigma=self.smoothing_bandwidth)
        
        # Find peaks (nucleosome centers) and troughs (linkers)
        peak_positions, trough_positions = self._find_peaks_troughs(
            smoothed_wps, start
        )
        
        # Compute periodicity
        periodicity = self._compute_periodicity(peak_positions)
        
        positions = np.arange(start, end)
        
        return WPSProfile(
            positions=positions,
            wps=wps,
            smoothed_wps=smoothed_wps,
            chrom=chrom,
            start=start,
            end=end,
            mean_wps=float(np.mean(wps)),
            std_wps=float(np.std(wps)),
            peak_positions=peak_positions,
            trough_positions=trough_positions,
            periodicity=periodicity,
        )
    
    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read is valid."""
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
        return True
    
    def _find_peaks_troughs(
        self,
        wps: NDArray,
        offset: int = 0,
    ) -> Tuple[NDArray, NDArray]:
        """Find nucleosome peaks and linker troughs."""
        # Find peaks
        peaks, _ = signal.find_peaks(
            wps,
            distance=100,  # Minimum distance between nucleosomes
            prominence=np.std(wps) * 0.5,
        )
        
        # Find troughs (inverse peaks)
        troughs, _ = signal.find_peaks(
            -wps,
            distance=50,
            prominence=np.std(wps) * 0.3,
        )
        
        return peaks + offset, troughs + offset
    
    def _compute_periodicity(self, peak_positions: NDArray) -> float:
        """Compute nucleosome periodicity from peak positions."""
        if len(peak_positions) < 3:
            return 0.0
        
        # Compute inter-peak distances
        distances = np.diff(peak_positions)
        
        if len(distances) == 0:
            return 0.0
        
        # Return median distance (robust to outliers)
        return float(np.median(distances))
    
    def compute_coverage_profile(
        self,
        bam_path: Union[str, Path],
        chrom: str,
        start: int,
        end: int,
    ) -> NDArray[np.int32]:
        """
        Compute simple fragment coverage profile.
        
        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file
        chrom : str
            Chromosome
        start : int
            Start position
        end : int
            End position
            
        Returns
        -------
        NDArray
            Coverage at each position
        """
        bam_path = Path(bam_path)
        region_size = end - start
        coverage = np.zeros(region_size, dtype=np.int32)
        
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(chrom, start, end):
                if not self._is_valid_read(read):
                    continue
                
                frag_start = max(0, read.reference_start - start)
                frag_end = min(region_size, read.reference_start + abs(read.template_length) - start)
                
                if frag_end > frag_start:
                    coverage[frag_start:frag_end] += 1
        
        return coverage


def compute_wps(
    bam_path: Union[str, Path],
    chrom: str,
    start: int,
    end: int,
    window_size: int = 120,
) -> WPSProfile:
    """
    Convenience function to compute WPS.
    
    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    chrom : str
        Chromosome
    start : int
        Start position
    end : int
        End position
    window_size : int, default 120
        WPS window size
        
    Returns
    -------
    WPSProfile
        Computed profile
    """
    analyzer = NucleosomeAnalyzer(window_size=window_size)
    return analyzer.compute_wps(bam_path, chrom, start, end)
