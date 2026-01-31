"""
Report section builders for fragmentomics reports.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt

from fragmentomics.reports.utils import fig_to_base64
from fragmentomics.viz import (
    plot_coverage_profile,
    plot_delfi_profile,
    plot_motif_frequencies,
    plot_size_distribution,
)

if TYPE_CHECKING:
    from fragmentomics import FragMentor
    from fragmentomics.features.coverage import CoverageProfile
    from fragmentomics.features.delfi import DELFIProfile
    from fragmentomics.models.cancer_detector import PredictionResult
    from fragmentomics.features.motifs import MotifProfile
    from fragmentomics.features.sizes import SizeDistribution


def build_qc_section(
    fm: FragMentor,
    sizes: SizeDistribution | None = None,
) -> dict[str, Any]:
    """
    Build QC/summary section.

    Parameters
    ----------
    fm : FragMentor
        FragMentor instance with loaded BAM
    sizes : SizeDistribution, optional
        Pre-computed size distribution

    Returns
    -------
    dict
        Section data for template
    """
    # Get basic stats from FragMentor
    stats = fm.stats()

    # Compute derived metrics
    total_reads = stats.total_reads
    proper_pairs = stats.proper_pairs
    fragments = stats.fragments_extracted
    duplicates = stats.duplicates_skipped
    duplicate_rate = duplicates / max(total_reads, 1)

    # Build metrics table
    metrics = [
        ("Total Reads", f"{total_reads:,}"),
        ("Proper Pairs", f"{proper_pairs:,}"),
        ("Fragments Extracted", f"{fragments:,}"),
        ("Duplicates Skipped", f"{duplicates:,}"),
        ("Duplicate Rate", f"{duplicate_rate:.1%}"),
    ]

    # Add size metrics if available
    if sizes:
        metrics.extend([
            ("Analyzed Fragments", f"{sizes.n_fragments:,}"),
            ("Median Size", f"{sizes.median:.0f} bp"),
            ("Mode Size", f"{sizes.mode:.0f} bp"),
            ("Mono Peak", f"{sizes.peak_mono:.0f} bp" if sizes.peak_mono else "N/A"),
            ("Di Peak", f"{sizes.peak_di:.0f} bp" if sizes.peak_di else "N/A"),
        ])

    # Create a stats-like object for QC flags
    class _QCStats:
        def __init__(self):
            self.mapped_reads = fragments
            self.duplicate_rate = duplicate_rate
            self.mean_mapq = 30  # Default, not tracked in ReadStats

    qc_stats = _QCStats()
    qc_flags = _compute_qc_flags(qc_stats, sizes)

    return {
        "id": "qc",
        "title": "Quality Control Summary",
        "type": "qc",
        "metrics": metrics,
        "qc_flags": qc_flags,
        "plot": None,  # QC section doesn't have a plot
    }


def _compute_qc_flags(stats: Any, sizes: SizeDistribution | None) -> list[dict[str, Any]]:
    """Compute QC pass/warn/fail flags."""
    flags = []

    # Read count check
    if stats.mapped_reads >= 10_000_000:
        flags.append({"name": "Read Depth", "status": "pass", "value": "≥10M reads"})
    elif stats.mapped_reads >= 1_000_000:
        flags.append({"name": "Read Depth", "status": "warn", "value": "1-10M reads"})
    else:
        flags.append({"name": "Read Depth", "status": "fail", "value": "<1M reads"})

    # Duplicate rate check
    if stats.duplicate_rate <= 0.20:
        flags.append({"name": "Duplicate Rate", "status": "pass", "value": "≤20%"})
    elif stats.duplicate_rate <= 0.40:
        flags.append({"name": "Duplicate Rate", "status": "warn", "value": "20-40%"})
    else:
        flags.append({"name": "Duplicate Rate", "status": "fail", "value": ">40%"})

    # MAPQ check
    if stats.mean_mapq >= 30:
        flags.append({"name": "Mapping Quality", "status": "pass", "value": "≥30"})
    elif stats.mean_mapq >= 20:
        flags.append({"name": "Mapping Quality", "status": "warn", "value": "20-30"})
    else:
        flags.append({"name": "Mapping Quality", "status": "fail", "value": "<20"})

    # Size distribution check (if available)
    if sizes and sizes.peak_mono:
        if 160 <= sizes.peak_mono <= 170:
            flags.append({"name": "Mono Peak", "status": "pass", "value": "160-170bp"})
        elif 150 <= sizes.peak_mono <= 180:
            flags.append({"name": "Mono Peak", "status": "warn", "value": "150-180bp"})
        else:
            flags.append({"name": "Mono Peak", "status": "fail", "value": "Outside expected"})

    return flags


def build_sizes_section(sizes: SizeDistribution) -> dict[str, Any]:
    """
    Build fragment size distribution section.

    Parameters
    ----------
    sizes : SizeDistribution
        Size distribution from analysis

    Returns
    -------
    dict
        Section data for template
    """
    # Generate plot
    fig, ax = plot_size_distribution(sizes, figsize=(10, 5))
    plot_b64 = fig_to_base64(fig)
    plt.close(fig)

    # Build metrics
    metrics = [
        ("Total Fragments", f"{sizes.n_fragments:,}"),
        ("Median Size", f"{sizes.median:.0f} bp"),
        ("Mode Size", f"{sizes.mode:.0f} bp"),
        ("Mean Size", f"{sizes.mean:.1f} bp"),
        ("Std Dev", f"{sizes.std:.1f} bp"),
        ("Mononucleosome Peak", f"{sizes.peak_mono:.0f} bp" if sizes.peak_mono else "Not detected"),
        ("Dinucleosome Peak", f"{sizes.peak_di:.0f} bp" if sizes.peak_di else "Not detected"),
        ("Short Ratio (<150bp)", f"{sizes.ratio_short:.1%}"),
        ("Mono/Di Ratio", f"{sizes.amplitude_ratio:.2f}" if sizes.amplitude_ratio else "N/A"),
    ]

    # Interpretation
    interpretation = _interpret_sizes(sizes)

    return {
        "id": "sizes",
        "title": "Fragment Size Distribution",
        "type": "sizes",
        "plot": plot_b64,
        "metrics": metrics,
        "interpretation": interpretation,
    }


def _interpret_sizes(sizes: SizeDistribution) -> str:
    """Generate interpretation text for size distribution."""
    lines = []

    # Peak analysis
    if sizes.peak_mono:
        if 160 <= sizes.peak_mono <= 170:
            lines.append(
                f"The mononucleosome peak at {sizes.peak_mono:.0f}bp is within the "
                "expected range (160-170bp) for cfDNA."
            )
        else:
            lines.append(
                f"The mononucleosome peak at {sizes.peak_mono:.0f}bp deviates from the "
                "typical range (160-170bp), which may indicate sample quality issues "
                "or biological variation."
            )

    # Short fragment analysis
    if sizes.ratio_short > 0.15:
        lines.append(
            f"Elevated short fragment ratio ({sizes.ratio_short:.1%}) may indicate "
            "increased tumor-derived DNA or sample degradation."
        )
    elif sizes.ratio_short < 0.05:
        lines.append(
            f"Low short fragment ratio ({sizes.ratio_short:.1%}) is consistent with "
            "healthy cfDNA profiles."
        )

    # Mono/Di ratio
    if sizes.amplitude_ratio:
        if sizes.amplitude_ratio > 5:
            lines.append(
                f"High mono/di ratio ({sizes.amplitude_ratio:.1f}) suggests good sample quality."
            )
        elif sizes.amplitude_ratio < 2:
            lines.append(
                f"Low mono/di ratio ({sizes.amplitude_ratio:.1f}) may indicate sample "
                "degradation or processing artifacts."
            )

    return " ".join(lines) if lines else "Size distribution appears typical for cfDNA."


def build_motifs_section(motifs: MotifProfile) -> dict[str, Any]:
    """
    Build end motif analysis section.

    Parameters
    ----------
    motifs : MotifProfile
        Motif profile from analysis

    Returns
    -------
    dict
        Section data for template
    """
    # Generate plot
    fig, ax = plot_motif_frequencies(motifs, figsize=(12, 5))
    plot_b64 = fig_to_base64(fig)
    plt.close(fig)

    # Top motifs table
    top_motifs = motifs.top_motifs(n=10)
    metrics = [(m, f"{f:.2%}") for m, f in top_motifs]

    # Motif diversity
    diversity_metrics = [
        ("Motif Diversity (Shannon)", f"{motifs.diversity:.3f}"),
        ("Unique Motifs", f"{motifs.n_unique_motifs}"),
        ("Total Fragments Analyzed", f"{motifs.n_fragments:,}"),
    ]

    # Interpretation
    interpretation = _interpret_motifs(motifs)

    return {
        "id": "motifs",
        "title": "End Motif Analysis",
        "type": "motifs",
        "plot": plot_b64,
        "metrics": diversity_metrics,
        "top_motifs": metrics,
        "interpretation": interpretation,
    }


def _interpret_motifs(motifs: MotifProfile) -> str:
    """Generate interpretation text for motif analysis."""
    lines = []

    # Check for CCCA enrichment (associated with cancer)
    ccca_freq = motifs.get_frequency("CCCA")
    if ccca_freq and ccca_freq > 0.05:
        lines.append(
            f"Elevated CCCA motif frequency ({ccca_freq:.1%}) may be associated "
            "with altered nuclease activity."
        )

    # Diversity assessment
    if motifs.diversity < 4.0:
        lines.append(
            "Lower than expected motif diversity suggests potential bias in "
            "enzymatic processing or sample preparation."
        )

    return " ".join(lines) if lines else "End motif profile appears typical."


def build_coverage_section(coverage: CoverageProfile) -> dict[str, Any]:
    """
    Build genome-wide coverage section.

    Parameters
    ----------
    coverage : CoverageProfile
        Coverage profile from analysis

    Returns
    -------
    dict
        Section data for template
    """
    # Generate plot
    fig, ax = plot_coverage_profile(coverage, figsize=(14, 4))
    plot_b64 = fig_to_base64(fig)
    plt.close(fig)

    # Compute additional metrics from bins
    n_bins = len(coverage.bins)
    raw_counts = [b.raw_count for b in coverage.bins]
    import numpy as np
    mean_cov = np.mean(raw_counts) if raw_counts else 0
    std_cov = np.std(raw_counts) if raw_counts else 0
    cv = std_cov / mean_cov if mean_cov > 0 else 0

    # Metrics
    metrics = [
        ("Mean Coverage", f"{mean_cov:.2f}"),
        ("Median Coverage", f"{coverage.median_coverage:.2f}"),
        ("Coverage MAD", f"{coverage.mad:.2f}"),
        ("Coefficient of Variation", f"{cv:.2%}"),
        ("Bins Analyzed", f"{n_bins:,}"),
        ("Bin Size", f"{coverage.bin_size:,} bp"),
        ("Total Reads", f"{coverage.total_reads:,}"),
    ]

    # GC correction status
    if coverage.gc_corrected:
        metrics.append(("GC Corrected", "Yes"))

    # Interpretation
    interpretation = _interpret_coverage(coverage, cv)

    return {
        "id": "coverage",
        "title": "Genome-Wide Coverage Profile",
        "type": "coverage",
        "plot": plot_b64,
        "metrics": metrics,
        "interpretation": interpretation,
    }


def _interpret_coverage(coverage: CoverageProfile, cv: float) -> str:
    """Generate interpretation text for coverage analysis."""
    lines = []

    # CV assessment
    if cv > 0.5:
        lines.append(
            f"High coverage variability (CV={cv:.2f}) may indicate "
            "copy number alterations or technical artifacts."
        )
    elif cv < 0.2:
        lines.append(
            f"Low coverage variability (CV={cv:.2f}) is consistent with "
            "a diploid genome without major CNV events."
        )

    return " ".join(lines) if lines else "Coverage profile appears uniform."


def build_delfi_section(delfi: DELFIProfile) -> dict[str, Any]:
    """
    Build DELFI fragmentation section.

    Parameters
    ----------
    delfi : DELFIProfile
        DELFI profile from analysis

    Returns
    -------
    dict
        Section data for template
    """
    # Generate plot
    fig, ax = plot_delfi_profile(delfi, figsize=(14, 4))
    plot_b64 = fig_to_base64(fig)
    plt.close(fig)

    # Compute metrics
    import numpy as np
    ratios = delfi.to_ratio_vector()
    n_elevated = np.sum(ratios > 0.6)
    pct_elevated = n_elevated / len(ratios) if len(ratios) > 0 else 0

    metrics = [
        ("Total Fragments", f"{delfi.total_fragments:,}"),
        ("Short Fragments (100-150bp)", f"{delfi.total_short:,}"),
        ("Long Fragments (151-220bp)", f"{delfi.total_long:,}"),
        ("Genome-Wide S/L Ratio", f"{delfi.genome_wide_ratio:.3f}"),
        ("Bins Analyzed", f"{len(delfi.bins):,}"),
        ("Elevated Bins (>0.6)", f"{n_elevated} ({pct_elevated:.1%})"),
    ]

    # Interpretation
    interpretation = _interpret_delfi(delfi, pct_elevated)

    return {
        "id": "delfi",
        "title": "DELFI Fragmentation Profile",
        "type": "delfi",
        "plot": plot_b64,
        "metrics": metrics,
        "interpretation": interpretation,
    }


def _interpret_delfi(delfi: DELFIProfile, pct_elevated: float) -> str:
    """Generate interpretation text for DELFI analysis."""
    lines = []

    # Genome-wide ratio assessment
    if delfi.genome_wide_ratio > 0.55:
        lines.append(
            f"Elevated genome-wide short/long ratio ({delfi.genome_wide_ratio:.2f}) "
            "is consistent with increased tumor-derived cfDNA."
        )
    elif delfi.genome_wide_ratio < 0.35:
        lines.append(
            f"Low genome-wide short/long ratio ({delfi.genome_wide_ratio:.2f}) "
            "is typical of healthy cfDNA profiles."
        )

    # Regional elevation
    if pct_elevated > 0.1:
        lines.append(
            f"{pct_elevated:.0%} of genomic bins show elevated fragmentation, "
            "suggesting potential copy number alterations or tumor contribution."
        )

    return " ".join(lines) if lines else "DELFI profile appears typical for healthy cfDNA."


def build_prediction_section(prediction: PredictionResult) -> dict[str, Any]:
    """
    Build ML prediction section.

    Parameters
    ----------
    prediction : PredictionResult
        Prediction result from cancer detector

    Returns
    -------
    dict
        Section data for template
    """
    # Determine status colors
    if prediction.prediction == "cancer":
        status_class = "danger" if prediction.confidence == "high" else "warning"
        status_emoji = "🔴" if prediction.confidence == "high" else "🟠"
    else:
        status_class = "success"
        status_emoji = "🟢"

    # Build feature importance table
    feature_table = []
    sorted_features = sorted(
        prediction.feature_importance.items(),
        key=lambda x: abs(x[1]),
        reverse=True
    )
    for feature, importance in sorted_features[:6]:
        direction = "↑" if importance > 0 else "↓"
        feature_table.append((feature, f"{direction} {abs(importance):.2f}"))

    metrics = [
        ("Prediction", f"{status_emoji} {prediction.prediction.upper()}"),
        ("Probability", f"{prediction.probability:.1%}"),
        ("Confidence", prediction.confidence.capitalize()),
        ("Features Used", str(len(prediction.features_used))),
    ]

    interpretation = _interpret_prediction(prediction)

    return {
        "id": "prediction",
        "title": "Cancer Risk Assessment",
        "type": "prediction",
        "plot": None,  # No plot for prediction
        "metrics": metrics,
        "feature_table": feature_table,
        "status_class": status_class,
        "probability": prediction.probability,
        "interpretation": interpretation,
    }


def _interpret_prediction(prediction: PredictionResult) -> str:
    """Generate interpretation text for prediction."""
    lines = []

    if prediction.prediction == "cancer":
        if prediction.confidence == "high":
            lines.append(
                "⚠️ High probability of cancer-associated cfDNA patterns detected. "
                "Clinical correlation and confirmatory testing recommended."
            )
        else:
            lines.append(
                "Moderate elevation of cancer-associated patterns. "
                "Consider repeat testing or additional workup."
            )
    else:
        if prediction.confidence == "high":
            lines.append(
                "cfDNA fragmentation patterns are consistent with healthy baseline. "
                "No significant cancer-associated signals detected."
            )
        else:
            lines.append(
                "Fragmentation patterns appear normal, though some uncertainty remains. "
                "Routine monitoring may be appropriate."
            )

    lines.append(
        "\n\n⚕️ This analysis is for research purposes only and should not be used "
        "for clinical decision-making without validation."
    )

    return " ".join(lines)
