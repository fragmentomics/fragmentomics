"""
FragMentor — The definitive toolkit for cfDNA fragmentomics analysis.

From BAM to biological insight in minutes. See what others miss. 🧬

Package: fragmentomics
Brand: FragMentor

Quick Start
-----------
>>> from fragmentomics import FragMentor
>>> fm = FragMentor("sample.bam")
>>> sizes = fm.sizes()
>>> print(sizes.summary())

Or use the CLI:
    $ fragmentomics sizes sample.bam -o results/
"""

__version__ = "0.1.0"
__author__ = "Jordon"
__brand__ = "FragMentor"

# Core imports
from fragmentomics.cli import app
from fragmentomics.features.sizes import (
    FragmentSizeAnalyzer,
    SizeDistribution,
    analyze_sizes,
)
from fragmentomics.io.bam import BamReader, Fragment, read_fragments
from fragmentomics.viz.plots import (
    plot_size_comparison,
    plot_size_distribution,
)


class FragMentor:
    """
    Main interface for cfDNA fragmentomics analysis.

    Provides a convenient wrapper around the core analysis modules.

    Parameters
    ----------
    bam_path : str
        Path to BAM/CRAM file
    reference : str, optional
        Path to reference genome FASTA
    min_mapq : int, default 30
        Minimum mapping quality filter

    Examples
    --------
    >>> fm = FragMentor("sample.bam", reference="hg38.fa")
    >>>
    >>> # Analyze fragment sizes
    >>> sizes = fm.sizes()
    >>> print(f"Median: {sizes.median:.0f} bp")
    >>> print(f"Short fragment ratio: {sizes.ratio_short:.1%}")
    >>>
    >>> # Plot distribution
    >>> fig, ax = fm.plot_sizes()
    >>> fig.savefig("sizes.png")
    """

    def __init__(
        self,
        bam_path: str,
        reference: str = None,
        min_mapq: int = 30,
        min_size: int = 50,
        max_size: int = 1000,
    ):
        self.bam_path = bam_path
        self.reference = reference
        self.min_mapq = min_mapq
        self.min_size = min_size
        self.max_size = max_size

        self._reader = BamReader(
            bam_path,
            reference=reference,
            min_mapq=min_mapq,
            min_size=min_size,
            max_size=max_size,
        )
        self._size_dist = None
        self._motif_profile = None
        self._coverage_profile = None

    def sizes(self, region: str = None, force: bool = False) -> SizeDistribution:
        """
        Analyze fragment size distribution.

        Parameters
        ----------
        region : str, optional
            Genomic region to analyze (e.g., "chr1:1000-2000")
        force : bool, default False
            Force re-analysis even if cached

        Returns
        -------
        SizeDistribution
            Size distribution with computed features
        """
        if self._size_dist is None or force:
            raw_sizes = self._reader.extract_sizes(region=region)
            analyzer = FragmentSizeAnalyzer(
                min_size=self.min_size,
                max_size=min(self.max_size, 500),
            )
            self._size_dist = analyzer.analyze(raw_sizes)
        return self._size_dist

    def end_motifs(self, k: int = 4, force: bool = False):
        """
        Analyze fragment end motifs.

        Parameters
        ----------
        k : int, default 4
            K-mer size for motif analysis
        force : bool, default False
            Force re-analysis even if cached

        Returns
        -------
        EndMotifProfile
            End motif profile with frequencies and diversity metrics

        Raises
        ------
        ValueError
            If no reference genome was provided
        """
        if self.reference is None:
            raise ValueError("Reference genome required for motif analysis")

        if self._motif_profile is None or force:
            from fragmentomics.features.motifs import EndMotifAnalyzer

            analyzer = EndMotifAnalyzer(
                k=k,
                min_mapq=self.min_mapq,
            )
            self._motif_profile = analyzer.analyze_bam(self.bam_path, self.reference)
        return self._motif_profile

    def coverage(self, bin_size: int = 100_000, force: bool = False):
        """
        Analyze genome-wide coverage.

        Parameters
        ----------
        bin_size : int, default 100000
            Size of genomic bins in bp
        force : bool, default False
            Force re-analysis even if cached

        Returns
        -------
        CoverageProfile
            Coverage profile with per-bin counts and metrics
        """
        if self._coverage_profile is None or force:
            from fragmentomics.features.coverage import CoverageAnalyzer

            analyzer = CoverageAnalyzer(
                bin_size=bin_size,
                min_mapq=self.min_mapq,
            )
            self._coverage_profile = analyzer.compute_coverage(
                self.bam_path,
                reference=self.reference,
            )
        return self._coverage_profile

    def delfi(self, bin_size: int = 5_000_000, force: bool = False):
        """
        Compute DELFI fragmentation profile.

        Implements the approach from Cristiano et al. 2019,
        computing short/long fragment ratios in genome-wide bins.

        Parameters
        ----------
        bin_size : int, default 5_000_000
            Size of genomic bins (5 Mb recommended)
        force : bool, default False
            Force re-analysis even if cached

        Returns
        -------
        DELFIProfile
            DELFI profile with per-bin ratios and features
        """
        if not hasattr(self, "_delfi_profile") or self._delfi_profile is None or force:
            from fragmentomics.features.delfi import DELFIAnalyzer

            analyzer = DELFIAnalyzer(
                bin_size=bin_size,
                min_mapq=self.min_mapq,
            )
            self._delfi_profile = analyzer.analyze(self.bam_path)
        return self._delfi_profile

    def predict(self, force: bool = False):
        """
        Predict cancer probability from fragmentomics features.

        Uses fragment size distribution and DELFI-style features
        to compute a cancer probability score.

        Parameters
        ----------
        force : bool, default False
            Force re-prediction even if cached

        Returns
        -------
        PredictionResult
            Prediction with probability, confidence, and feature importance
        """
        if not hasattr(self, "_prediction") or self._prediction is None or force:
            from fragmentomics.models.cancer_detector import CancerDetector

            detector = CancerDetector()
            self._prediction = detector.predict_from_bam(self.bam_path)
        return self._prediction

    def plot_sizes(self, **kwargs):
        """
        Plot fragment size distribution.

        Parameters
        ----------
        **kwargs
            Arguments passed to plot_size_distribution

        Returns
        -------
        tuple[Figure, Axes]
            Matplotlib figure and axes
        """
        dist = self.sizes()
        return plot_size_distribution(dist, **kwargs)

    def stats(self):
        """
        Get BAM file statistics.

        Returns
        -------
        ReadStats
            Statistics including read counts, mapping rates, etc.
        """
        return self._reader.stats

    def __repr__(self):
        return f"FragMentor('{self.bam_path}')"


# Lazy import for reports (to avoid slow startup)
def __getattr__(name):
    if name == "ReportGenerator":
        from fragmentomics.reports import ReportGenerator
        return ReportGenerator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__brand__",
    # Main interface
    "FragMentor",
    # I/O
    "BamReader",
    "read_fragments",
    "Fragment",
    # Analysis
    "FragmentSizeAnalyzer",
    "SizeDistribution",
    "analyze_sizes",
    # Visualization
    "plot_size_distribution",
    "plot_size_comparison",
    # Reports
    "ReportGenerator",
    # CLI
    "app",
]
