"""
Feature extraction modules for cfDNA fragmentomics.
"""

from fragmentomics.features.sizes import (
    FragmentSizeAnalyzer,
    SizeDistribution,
    analyze_sizes,
)
from fragmentomics.features.motifs import (
    EndMotifAnalyzer,
    EndMotifProfile,
    analyze_end_motifs,
    ALL_4MERS,
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
]
