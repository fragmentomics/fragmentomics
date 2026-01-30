"""
Genomic region handling (BED files, intervals).

Provides utilities for reading BED files and managing genomic intervals
for targeted analysis.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class GenomicRegion:
    """
    A genomic interval.

    Attributes
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position (0-based)
    end : int
        End position (exclusive)
    name : str, optional
        Region name/label
    score : float, optional
        Associated score
    strand : str, optional
        Strand (+, -, or .)
    """

    chrom: str
    start: int
    end: int
    name: str = ""
    score: float = 0.0
    strand: str = "."

    @property
    def length(self) -> int:
        """Region length in bp."""
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        """Midpoint of the region."""
        return (self.start + self.end) // 2

    def __str__(self) -> str:
        if self.name:
            return f"{self.chrom}:{self.start}-{self.end} ({self.name})"
        return f"{self.chrom}:{self.start}-{self.end}"

    def to_ucsc(self) -> str:
        """Return UCSC-style coordinate string."""
        return f"{self.chrom}:{self.start + 1}-{self.end}"

    def overlaps(self, other: GenomicRegion) -> bool:
        """Check if this region overlaps another."""
        if self.chrom != other.chrom:
            return False
        return self.start < other.end and other.start < self.end

    def contains(self, chrom: str, pos: int) -> bool:
        """Check if region contains a position."""
        return self.chrom == chrom and self.start <= pos < self.end

    def expand(self, upstream: int = 0, downstream: int = 0) -> GenomicRegion:
        """
        Expand region by specified amounts.

        Parameters
        ----------
        upstream : int
            Bases to add upstream (respects strand)
        downstream : int
            Bases to add downstream (respects strand)

        Returns
        -------
        GenomicRegion
            Expanded region
        """
        if self.strand == "-":
            new_start = max(0, self.start - downstream)
            new_end = self.end + upstream
        else:
            new_start = max(0, self.start - upstream)
            new_end = self.end + downstream

        return GenomicRegion(
            chrom=self.chrom,
            start=new_start,
            end=new_end,
            name=self.name,
            score=self.score,
            strand=self.strand,
        )


def read_bed(
    path: str | Path,
    skip_header: bool = False,
) -> Iterator[GenomicRegion]:
    """
    Read genomic regions from a BED file.

    Supports BED3, BED4, BED5, BED6, and BED12 formats.

    Parameters
    ----------
    path : str or Path
        Path to BED file (can be gzipped)
    skip_header : bool, default False
        Skip lines starting with # or track/browser

    Yields
    ------
    GenomicRegion
        Each region from the file

    Examples
    --------
    >>> for region in read_bed("promoters.bed"):
    ...     print(f"Processing {region}")
    """
    import gzip

    path = Path(path)

    # Handle gzipped files
    if path.suffix == ".gz":
        opener = gzip.open
        mode = "rt"
    else:
        opener = open
        mode = "r"

    with opener(path, mode) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Skip headers
            if skip_header or line.startswith(("#", "track", "browser")):
                if line.startswith(("#", "track", "browser")):
                    continue

            fields = line.split("\t")

            # Minimum BED3
            if len(fields) < 3:
                logger.warning(f"Skipping line {line_num}: fewer than 3 fields")
                continue

            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])

                # Optional fields
                name = fields[3] if len(fields) > 3 else ""
                score = float(fields[4]) if len(fields) > 4 else 0.0
                strand = fields[5] if len(fields) > 5 else "."

                yield GenomicRegion(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score,
                    strand=strand,
                )

            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping line {line_num}: {e}")
                continue


def read_bed_to_list(path: str | Path, skip_header: bool = True) -> list[GenomicRegion]:
    """
    Read all regions from BED file into a list.

    Parameters
    ----------
    path : str or Path
        Path to BED file
    skip_header : bool, default True
        Skip header lines

    Returns
    -------
    list[GenomicRegion]
        All regions from the file
    """
    return list(read_bed(path, skip_header=skip_header))


def write_bed(
    regions: list[GenomicRegion],
    path: str | Path,
    bed_format: int = 6,
) -> None:
    """
    Write regions to a BED file.

    Parameters
    ----------
    regions : list[GenomicRegion]
        Regions to write
    path : str or Path
        Output path
    bed_format : int, default 6
        Number of columns (3, 4, 5, or 6)
    """
    path = Path(path)

    with open(path, "w") as f:
        for region in regions:
            fields = [region.chrom, str(region.start), str(region.end)]

            if bed_format >= 4:
                fields.append(region.name)
            if bed_format >= 5:
                fields.append(str(int(region.score)))
            if bed_format >= 6:
                fields.append(region.strand)

            f.write("\t".join(fields) + "\n")


def parse_region_string(region: str) -> GenomicRegion:
    """
    Parse a region string like 'chr1:1000-2000'.

    Supports formats:
    - chr1:1000-2000
    - chr1:1,000-2,000
    - chr1:1000..2000

    Parameters
    ----------
    region : str
        Region string

    Returns
    -------
    GenomicRegion
        Parsed region

    Raises
    ------
    ValueError
        If region string is invalid
    """
    # Remove commas from numbers
    region = region.replace(",", "")

    # Handle different separators
    region = region.replace("..", "-")

    try:
        chrom, coords = region.split(":")
        start_str, end_str = coords.split("-")

        return GenomicRegion(
            chrom=chrom,
            start=int(start_str),
            end=int(end_str),
        )
    except (ValueError, IndexError) as e:
        raise ValueError(f"Invalid region format: {region}") from e


def merge_overlapping(
    regions: list[GenomicRegion],
    gap: int = 0,
) -> list[GenomicRegion]:
    """
    Merge overlapping or adjacent regions.

    Parameters
    ----------
    regions : list[GenomicRegion]
        Input regions (will be sorted)
    gap : int, default 0
        Merge regions within this distance

    Returns
    -------
    list[GenomicRegion]
        Merged regions
    """
    if not regions:
        return []

    # Sort by chromosome, then start position
    sorted_regions = sorted(regions, key=lambda r: (r.chrom, r.start))

    merged = [sorted_regions[0]]

    for region in sorted_regions[1:]:
        last = merged[-1]

        # Check if same chromosome and overlapping/adjacent
        if region.chrom == last.chrom and region.start <= last.end + gap:
            # Extend the last region
            merged[-1] = GenomicRegion(
                chrom=last.chrom,
                start=last.start,
                end=max(last.end, region.end),
                name=last.name if last.name else region.name,
                strand=last.strand if last.strand != "." else region.strand,
            )
        else:
            merged.append(region)

    return merged
