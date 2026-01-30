"""
I/O modules for reading BAM/CRAM files and writing outputs.
"""

from fragmentomics.io.bam import BamReader, read_fragments

__all__ = ["BamReader", "read_fragments"]
