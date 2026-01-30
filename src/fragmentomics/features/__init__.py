"""
Feature extraction modules for cfDNA fragmentomics.
"""

from fragmentomics.features.coverage import (
    CoverageAnalyzer,
    CoverageBin,
    CoverageProfile,
    compute_coverage,
)
from fragmentomics.features.delfi import (
    DELFIAnalyzer,
    DELFIBin,
    DELFIProfile,
    compute_delfi_profile,
    create_healthy_reference,
)
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
    # Coverage/CNV
    "CoverageAnalyzer",
    "CoverageBin",
    "CoverageProfile",
    "compute_coverage",
    # DELFI (Cristiano et al.)
    "DELFIAnalyzer",
    "DELFIBin",
    "DELFIProfile",
    "compute_delfi_profile",
    "create_healthy_reference",
]
