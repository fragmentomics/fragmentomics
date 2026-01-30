"""
BAM/CRAM file reading utilities.

Provides streaming access to fragment data from aligned sequencing reads.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Union

import pysam
import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


@dataclass
class Fragment:
    """Represents a single cfDNA fragment."""
    
    chrom: str
    start: int
    end: int
    size: int
    mapq: int
    gc_content: Optional[float] = None
    
    @property
    def midpoint(self) -> int:
        """Return the midpoint of the fragment."""
        return (self.start + self.end) // 2


@dataclass 
class ReadStats:
    """Statistics from BAM reading."""
    
    total_reads: int = 0
    proper_pairs: int = 0
    duplicates_skipped: int = 0
    low_mapq_skipped: int = 0
    fragments_extracted: int = 0
    
    def __str__(self) -> str:
        return (
            f"ReadStats:\n"
            f"  Total reads seen: {self.total_reads:,}\n"
            f"  Proper pairs: {self.proper_pairs:,}\n"
            f"  Duplicates skipped: {self.duplicates_skipped:,}\n"
            f"  Low MAPQ skipped: {self.low_mapq_skipped:,}\n"
            f"  Fragments extracted: {self.fragments_extracted:,}"
        )


class BamReader:
    """
    Streaming BAM/CRAM reader optimized for cfDNA fragment extraction.
    
    Extracts properly-paired reads and calculates fragment sizes from
    template lengths. Handles both BAM and CRAM formats transparently.
    
    Parameters
    ----------
    path : str or Path
        Path to BAM or CRAM file
    reference : str or Path, optional
        Path to reference FASTA (required for CRAM, optional for BAM)
    min_mapq : int, default 30
        Minimum mapping quality to include
    min_size : int, default 50
        Minimum fragment size to include
    max_size : int, default 1000
        Maximum fragment size to include
    skip_duplicates : bool, default True
        Whether to skip duplicate reads
        
    Examples
    --------
    >>> reader = BamReader("sample.bam", min_mapq=30)
    >>> sizes = reader.extract_sizes()
    >>> print(f"Median size: {np.median(sizes):.0f} bp")
    """
    
    def __init__(
        self,
        path: Union[str, Path],
        reference: Optional[Union[str, Path]] = None,
        min_mapq: int = 30,
        min_size: int = 50,
        max_size: int = 1000,
        skip_duplicates: bool = True,
    ):
        self.path = Path(path)
        self.reference = Path(reference) if reference else None
        self.min_mapq = min_mapq
        self.min_size = min_size
        self.max_size = max_size
        self.skip_duplicates = skip_duplicates
        self.stats = ReadStats()
        
        # Validate file exists
        if not self.path.exists():
            raise FileNotFoundError(f"BAM file not found: {self.path}")
        
        # Check for index
        index_extensions = ['.bai', '.csi']
        has_index = any(
            self.path.with_suffix(self.path.suffix + ext).exists() or
            Path(str(self.path) + ext).exists()
            for ext in index_extensions
        )
        if not has_index:
            logger.warning(f"No index found for {self.path}. Random access will be slow.")
    
    def _open_bam(self) -> pysam.AlignmentFile:
        """Open the BAM/CRAM file with appropriate settings."""
        kwargs = {"mode": "rb"}
        if self.reference:
            kwargs["reference_filename"] = str(self.reference)
        return pysam.AlignmentFile(str(self.path), **kwargs)
    
    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if a read passes all filters."""
        self.stats.total_reads += 1
        
        # Must be proper pair
        if not read.is_proper_pair:
            return False
        self.stats.proper_pairs += 1
        
        # Skip duplicates
        if self.skip_duplicates and read.is_duplicate:
            self.stats.duplicates_skipped += 1
            return False
        
        # Check mapping quality
        if read.mapping_quality < self.min_mapq:
            self.stats.low_mapq_skipped += 1
            return False
        
        # Only process read1 to avoid counting fragments twice
        if not read.is_read1:
            return False
        
        # Check fragment size is in range
        size = abs(read.template_length)
        if size < self.min_size or size > self.max_size:
            return False
        
        return True
    
    def iter_fragments(
        self,
        region: Optional[str] = None,
    ) -> Iterator[Fragment]:
        """
        Iterate over fragments in the BAM file.
        
        Parameters
        ----------
        region : str, optional
            Genomic region to fetch (e.g., "chr1:1000-2000")
            
        Yields
        ------
        Fragment
            Fragment objects with position and size information
        """
        with self._open_bam() as bam:
            iterator = bam.fetch(region=region) if region else bam.fetch()
            
            for read in iterator:
                if not self._is_valid_read(read):
                    continue
                
                size = abs(read.template_length)
                start = min(read.reference_start, read.next_reference_start)
                end = start + size
                
                self.stats.fragments_extracted += 1
                
                yield Fragment(
                    chrom=read.reference_name,
                    start=start,
                    end=end,
                    size=size,
                    mapq=read.mapping_quality,
                )
    
    def extract_sizes(
        self,
        region: Optional[str] = None,
        max_fragments: Optional[int] = None,
    ) -> NDArray[np.int32]:
        """
        Extract fragment sizes as a numpy array.
        
        Parameters
        ----------
        region : str, optional
            Genomic region to fetch
        max_fragments : int, optional
            Maximum number of fragments to extract (for testing/preview)
            
        Returns
        -------
        NDArray[np.int32]
            Array of fragment sizes in base pairs
        """
        sizes = []
        for i, fragment in enumerate(self.iter_fragments(region)):
            sizes.append(fragment.size)
            if max_fragments and i + 1 >= max_fragments:
                break
        
        return np.array(sizes, dtype=np.int32)
    
    def extract_positions(
        self,
        region: Optional[str] = None,
    ) -> tuple[NDArray, NDArray, NDArray]:
        """
        Extract fragment positions and sizes.
        
        Returns
        -------
        tuple of (chroms, midpoints, sizes)
        """
        chroms = []
        midpoints = []
        sizes = []
        
        for fragment in self.iter_fragments(region):
            chroms.append(fragment.chrom)
            midpoints.append(fragment.midpoint)
            sizes.append(fragment.size)
        
        return (
            np.array(chroms),
            np.array(midpoints, dtype=np.int64),
            np.array(sizes, dtype=np.int32),
        )


def read_fragments(
    bam_path: Union[str, Path],
    reference: Optional[Union[str, Path]] = None,
    region: Optional[str] = None,
    min_mapq: int = 30,
    min_size: int = 50,
    max_size: int = 1000,
) -> NDArray[np.int32]:
    """
    Convenience function to extract fragment sizes from a BAM file.
    
    Parameters
    ----------
    bam_path : str or Path
        Path to BAM/CRAM file
    reference : str or Path, optional
        Path to reference FASTA
    region : str, optional
        Genomic region to fetch
    min_mapq : int, default 30
        Minimum mapping quality
    min_size : int, default 50
        Minimum fragment size
    max_size : int, default 1000
        Maximum fragment size
        
    Returns
    -------
    NDArray[np.int32]
        Array of fragment sizes
        
    Examples
    --------
    >>> sizes = read_fragments("sample.bam")
    >>> print(f"Found {len(sizes):,} fragments")
    >>> print(f"Median size: {np.median(sizes):.0f} bp")
    """
    reader = BamReader(
        bam_path,
        reference=reference,
        min_mapq=min_mapq,
        min_size=min_size,
        max_size=max_size,
    )
    return reader.extract_sizes(region=region)
