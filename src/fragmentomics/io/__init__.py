"""
I/O modules for reading BAM/CRAM files, BED files, and writing outputs.
"""

from fragmentomics.io.bam import BamReader, read_fragments
from fragmentomics.io.regions import (
    GenomicRegion,
    merge_overlapping,
    parse_region_string,
    read_bed,
    read_bed_to_list,
    write_bed,
)

__all__ = [
    "BamReader",
    "read_fragments",
    "GenomicRegion",
    "read_bed",
    "read_bed_to_list",
    "write_bed",
    "parse_region_string",
    "merge_overlapping",
]
