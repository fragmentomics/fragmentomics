"""
Tests for end motif analysis module.
"""

import numpy as np
import pytest
from collections import Counter

from fragmentomics.features.motifs import (
    EndMotifAnalyzer,
    EndMotifProfile,
    generate_kmers,
    ALL_4MERS,
    KMER_TO_IDX,
)


class TestGenerateKmers:
    """Tests for k-mer generation."""
    
    def test_1mers(self):
        """Test 1-mer generation."""
        kmers = generate_kmers(1)
        assert kmers == ['A', 'C', 'G', 'T']
    
    def test_2mers(self):
        """Test 2-mer generation."""
        kmers = generate_kmers(2)
        assert len(kmers) == 16
        assert 'AA' in kmers
        assert 'TT' in kmers
    
    def test_4mers(self):
        """Test 4-mer generation."""
        assert len(ALL_4MERS) == 256
        assert 'AAAA' in ALL_4MERS
        assert 'TTTT' in ALL_4MERS
        assert 'ACGT' in ALL_4MERS
    
    def test_kmer_index(self):
        """Test k-mer to index mapping."""
        assert len(KMER_TO_IDX) == 256
        assert KMER_TO_IDX['AAAA'] == 0
        # Each k-mer should have unique index
        assert len(set(KMER_TO_IDX.values())) == 256


class TestEndMotifProfile:
    """Tests for EndMotifProfile dataclass."""
    
    @pytest.fixture
    def sample_profile(self):
        """Create a sample profile for testing."""
        counts = Counter({
            'CCCA': 100,
            'CCCC': 80,
            'CCCT': 60,
            'CCCG': 40,
            'AAAA': 20,
        })
        total = sum(counts.values())
        freqs = {k: v / total for k, v in counts.items()}
        
        return EndMotifProfile(
            motif_counts=counts,
            motif_freqs=freqs,
            n_fragments=150,
            n_ends=300,
            k=4,
            entropy=2.1,
            diversity=4.3,
            gc_content=0.65,
            top_motifs=[('CCCA', 100/total), ('CCCC', 80/total)],
        )
    
    def test_to_vector(self, sample_profile):
        """Test conversion to feature vector."""
        vec = sample_profile.to_vector(normalize=True)
        
        assert isinstance(vec, np.ndarray)
        assert len(vec) == 256
        assert vec.sum() > 0
        assert vec.sum() <= 1.0 + 1e-6  # Frequencies sum to ~1
    
    def test_to_vector_counts(self, sample_profile):
        """Test conversion to raw count vector."""
        vec = sample_profile.to_vector(normalize=False)
        
        assert vec.sum() == 300  # Total counts
    
    def test_to_dict(self, sample_profile):
        """Test serialization to dict."""
        d = sample_profile.to_dict()
        
        assert isinstance(d, dict)
        assert d['n_fragments'] == 150
        assert d['k'] == 4
        assert 'entropy' in d
        assert 'motif_freqs' in d
    
    def test_summary(self, sample_profile):
        """Test summary string generation."""
        summary = sample_profile.summary()
        
        assert isinstance(summary, str)
        assert 'End Motif Profile' in summary
        assert 'Entropy' in summary
        assert 'CCCA' in summary


class TestEndMotifAnalyzer:
    """Tests for EndMotifAnalyzer class."""
    
    def test_init_defaults(self):
        """Test default initialization."""
        analyzer = EndMotifAnalyzer()
        
        assert analyzer.k == 4
        assert analyzer.min_mapq == 30
        assert analyzer.both_ends is True
    
    def test_init_custom(self):
        """Test custom initialization."""
        analyzer = EndMotifAnalyzer(k=6, min_mapq=20, both_ends=False)
        
        assert analyzer.k == 6
        assert analyzer.min_mapq == 20
        assert analyzer.both_ends is False
    
    def test_reverse_complement(self):
        """Test reverse complement function."""
        rc = EndMotifAnalyzer._reverse_complement
        
        assert rc('ACGT') == 'ACGT'  # Palindrome
        assert rc('AAAA') == 'TTTT'
        assert rc('CCCC') == 'GGGG'
        assert rc('AACG') == 'CGTT'
    
    def test_compute_entropy(self):
        """Test entropy computation."""
        entropy = EndMotifAnalyzer._compute_entropy
        
        # Uniform distribution has max entropy
        uniform = {f'M{i}': 0.25 for i in range(4)}
        assert entropy(uniform) == pytest.approx(2.0, abs=0.01)
        
        # Single motif has zero entropy
        single = {'AAAA': 1.0}
        assert entropy(single) == 0.0
    
    def test_compute_gc_content(self):
        """Test GC content computation."""
        gc = EndMotifAnalyzer._compute_gc_content
        
        # All GC
        all_gc = {'GGGG': 0.5, 'CCCC': 0.5}
        assert gc(all_gc) == 1.0
        
        # All AT
        all_at = {'AAAA': 0.5, 'TTTT': 0.5}
        assert gc(all_at) == 0.0
        
        # Mixed
        mixed = {'ACGT': 1.0}
        assert gc(mixed) == 0.5
    
    def test_file_not_found(self):
        """Test error on missing file."""
        analyzer = EndMotifAnalyzer()
        
        with pytest.raises(FileNotFoundError):
            analyzer.analyze_bam("nonexistent.bam", "ref.fa")
