"""
Fragment end motif analysis.

Analyzes the sequence composition at cfDNA fragment ends, which reflects
enzymatic cleavage patterns and can distinguish tissue origins.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


# Standard 4-mer motifs (256 total)
def generate_kmers(k: int = 4) -> list[str]:
    """Generate all possible k-mers."""
    bases = ["A", "C", "G", "T"]
    if k == 1:
        return bases
    return [b + kmer for b in bases for kmer in generate_kmers(k - 1)]


ALL_4MERS = generate_kmers(4)
KMER_TO_IDX = {kmer: i for i, kmer in enumerate(ALL_4MERS)}


@dataclass
class EndMotifProfile:
    """
    Fragment end motif profile.

    Attributes
    ----------
    motif_counts : Counter
        Raw counts of each motif
    motif_freqs : dict
        Normalized frequencies of each motif
    n_fragments : int
        Total fragments analyzed
    k : int
        K-mer size

    Derived Features
    ----------------
    entropy : float
        Shannon entropy of motif distribution
    top_motifs : list
        Most frequent motifs
    gc_content : float
        Average GC content at ends
    diversity : float
        Effective number of motifs (exp of entropy)
    """

    motif_counts: Counter = field(default_factory=Counter)
    motif_freqs: dict = field(default_factory=dict)
    n_fragments: int = 0
    n_ends: int = 0
    k: int = 4

    # Derived features
    entropy: float = 0.0
    diversity: float = 0.0
    gc_content: float = 0.0
    top_motifs: list = field(default_factory=list)

    def to_vector(self, normalize: bool = True) -> NDArray[np.float64]:
        """
        Convert motif profile to fixed-length feature vector.

        Parameters
        ----------
        normalize : bool, default True
            Return frequencies (True) or raw counts (False)

        Returns
        -------
        NDArray
            256-element vector for 4-mers
        """
        if self.k != 4:
            raise ValueError("to_vector only supports k=4")

        vec = np.zeros(256, dtype=np.float64)
        source = self.motif_freqs if normalize else self.motif_counts

        for motif, value in source.items():
            if motif in KMER_TO_IDX:
                vec[KMER_TO_IDX[motif]] = value

        return vec

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "n_fragments": self.n_fragments,
            "n_ends": self.n_ends,
            "k": self.k,
            "entropy": self.entropy,
            "diversity": self.diversity,
            "gc_content": self.gc_content,
            "top_motifs": self.top_motifs[:10],
            "motif_freqs": dict(self.motif_freqs),
        }

    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            "End Motif Profile Summary",
            "=" * 40,
            f"Fragments analyzed: {self.n_fragments:,}",
            f"Total ends: {self.n_ends:,}",
            f"K-mer size: {self.k}",
            "",
            "Features:",
            f"  Entropy: {self.entropy:.3f} bits",
            f"  Diversity: {self.diversity:.1f} effective motifs",
            f"  GC content: {self.gc_content:.1%}",
            "",
            "Top 10 motifs:",
        ]

        for motif, freq in self.top_motifs[:10]:
            lines.append(f"  {motif}: {freq:.3%}")

        return "\n".join(lines)


class EndMotifAnalyzer:
    """
    Analyzer for cfDNA fragment end motifs.

    Extracts k-mer sequences at fragment ends and computes
    frequency profiles and derived features.

    Parameters
    ----------
    k : int, default 4
        K-mer length
    min_mapq : int, default 30
        Minimum mapping quality
    both_ends : bool, default True
        Analyze both 5' and 3' ends

    Examples
    --------
    >>> analyzer = EndMotifAnalyzer(k=4)
    >>> profile = analyzer.analyze_bam("sample.bam", "reference.fa")
    >>> print(profile.summary())
    >>> vec = profile.to_vector()  # For ML
    """

    def __init__(
        self,
        k: int = 4,
        min_mapq: int = 30,
        both_ends: bool = True,
    ):
        self.k = k
        self.min_mapq = min_mapq
        self.both_ends = both_ends

    def analyze_bam(
        self,
        bam_path: str | Path,
        reference_path: str | Path,
        region: str | None = None,
        max_fragments: int | None = None,
    ) -> EndMotifProfile:
        """
        Analyze end motifs from a BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM/CRAM file
        reference_path : str or Path
            Path to reference FASTA
        region : str, optional
            Genomic region to analyze
        max_fragments : int, optional
            Maximum fragments to analyze

        Returns
        -------
        EndMotifProfile
            Motif profile with computed features
        """
        bam_path = Path(bam_path)
        reference_path = Path(reference_path)

        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        if not reference_path.exists():
            raise FileNotFoundError(f"Reference not found: {reference_path}")

        motif_counts = Counter()
        n_fragments = 0
        n_ends = 0

        ref = pysam.FastaFile(str(reference_path))

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                iterator = bam.fetch(region=region) if region else bam.fetch()

                for read in iterator:
                    if not self._is_valid_read(read):
                        continue

                    # Extract motifs at fragment ends
                    motifs = self._extract_motifs(read, ref)

                    for motif in motifs:
                        if motif and len(motif) == self.k:
                            motif_counts[motif] += 1
                            n_ends += 1

                    n_fragments += 1

                    if max_fragments and n_fragments >= max_fragments:
                        break
        finally:
            ref.close()

        if n_ends == 0:
            raise ValueError("No valid motifs extracted")

        # Compute frequencies
        total = sum(motif_counts.values())
        motif_freqs = {k: v / total for k, v in motif_counts.items()}

        # Compute derived features
        entropy = self._compute_entropy(motif_freqs)
        diversity = np.exp(entropy)
        gc_content = self._compute_gc_content(motif_freqs)
        top_motifs = motif_counts.most_common(20)
        top_motifs = [(m, c / total) for m, c in top_motifs]

        return EndMotifProfile(
            motif_counts=motif_counts,
            motif_freqs=motif_freqs,
            n_fragments=n_fragments,
            n_ends=n_ends,
            k=self.k,
            entropy=entropy,
            diversity=diversity,
            gc_content=gc_content,
            top_motifs=top_motifs,
        )

    def _is_valid_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read passes filters."""
        if read.is_unmapped:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if not read.is_proper_pair:
            return False
        if read.is_duplicate:
            return False
        if not read.is_read1:
            return False
        return True

    def _extract_motifs(
        self,
        read: pysam.AlignedSegment,
        ref: pysam.FastaFile,
    ) -> list[str]:
        """Extract k-mer motifs at fragment ends."""
        motifs = []
        chrom = read.reference_name

        try:
            chrom_len = ref.get_reference_length(chrom)
        except KeyError:
            return motifs

        # 5' end of fragment
        if read.is_reverse:
            # Read is on reverse strand, 5' end is at reference_end
            pos = read.reference_end
            if pos and pos + self.k <= chrom_len:
                seq = ref.fetch(chrom, pos, pos + self.k)
                # Reverse complement for 5' end on reverse strand
                seq = self._reverse_complement(seq)
                motifs.append(seq.upper())
        else:
            # Read is on forward strand, 5' end is at reference_start
            pos = read.reference_start
            if pos >= self.k:
                seq = ref.fetch(chrom, pos - self.k, pos)
                motifs.append(seq.upper())

        if self.both_ends:
            # 3' end of fragment (opposite end)
            if read.is_reverse:
                mate_pos = read.next_reference_start
                if mate_pos >= self.k:
                    seq = ref.fetch(chrom, mate_pos - self.k, mate_pos)
                    motifs.append(seq.upper())
            else:
                mate_end = read.reference_start + abs(read.template_length)
                if mate_end + self.k <= chrom_len:
                    seq = ref.fetch(chrom, mate_end, mate_end + self.k)
                    seq = self._reverse_complement(seq)
                    motifs.append(seq.upper())

        return motifs

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement of sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement.get(b, "N") for b in reversed(seq))

    @staticmethod
    def _compute_entropy(freqs: dict) -> float:
        """Compute Shannon entropy of distribution."""
        entropy = 0.0
        for freq in freqs.values():
            if freq > 0:
                entropy -= freq * np.log2(freq)
        return entropy

    @staticmethod
    def _compute_gc_content(freqs: dict) -> float:
        """Compute average GC content from motif frequencies."""
        total_gc = 0.0
        total_weight = 0.0

        for motif, freq in freqs.items():
            gc = sum(1 for b in motif if b in "GC")
            total_gc += gc * freq
            total_weight += len(motif) * freq

        return total_gc / total_weight if total_weight > 0 else 0.0


def analyze_end_motifs(
    bam_path: str | Path,
    reference_path: str | Path,
    k: int = 4,
    region: str | None = None,
) -> EndMotifProfile:
    """
    Convenience function to analyze end motifs.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    reference_path : str or Path
        Path to reference FASTA
    k : int, default 4
        K-mer length
    region : str, optional
        Genomic region to analyze

    Returns
    -------
    EndMotifProfile
        Motif profile with features

    Examples
    --------
    >>> profile = analyze_end_motifs("sample.bam", "hg38.fa")
    >>> print(f"Entropy: {profile.entropy:.2f}")
    >>> print(f"Top motif: {profile.top_motifs[0]}")
    """
    analyzer = EndMotifAnalyzer(k=k)
    return analyzer.analyze_bam(bam_path, reference_path, region=region)
