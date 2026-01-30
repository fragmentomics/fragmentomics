"""
Feature extraction modules for cfDNA fragmentomics.
"""

from fragmentomics.features.gc import (
    GCBiasProfile,
    GCCorrector,
    compute_gc_bias,
)
from fragmentomics.features.motifs import (
    ALL_4MERS,
    EndMotifAnalyzer,
    EndMotifProfile,
    analyze_end_motifs,
)
from fragmentomics.features.nucleosome import (
    NucleosomeAnalyzer,
    WPSProfile,
    compute_wps,
)
from fragmentomics.features.sizes import (
    FragmentSizeAnalyzer,
    SizeDistribution,
    analyze_sizes,
)

__all__ = [
    # Size analysis
    "FragmentSizeAnalyzer",
    "SizeDistribution",
    "analyze_sizes",
    # Motif analysis
    "EndMotifAnalyzer",
    "EndMotifProfile",
    "analyze_end_motifs",
    "ALL_4MERS",
    # GC correction
    "GCCorrector",
    "GCBiasProfile",
    "compute_gc_bias",
    # Nucleosome positioning
    "NucleosomeAnalyzer",
    "WPSProfile",
    "compute_wps",
]
