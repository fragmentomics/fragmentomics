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

__version__ = "0.1.0-dev"
__author__ = "Jordon"
__brand__ = "FragMentor"

# Core imports
from fragmentomics.io.bam import BamReader, read_fragments, Fragment
from fragmentomics.features.sizes import (
    FragmentSizeAnalyzer,
    SizeDistribution,
    analyze_sizes,
)
from fragmentomics.viz.plots import (
    plot_size_distribution,
    plot_size_comparison,
)
from fragmentomics.cli import app


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
    
    @property
    def stats(self):
        """Return read statistics from BAM parsing."""
        return self._reader.stats
    
    def __repr__(self):
        return f"FragMentor('{self.bam_path}')"


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
    # CLI
    "app",
]
