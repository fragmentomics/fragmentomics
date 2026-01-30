"""
Tests for GC bias correction module.
"""

import numpy as np
import pytest

from fragmentomics.features.gc import (
    GCCorrector,
    GCBiasProfile,
)


class TestGCCorrector:
    """Tests for GCCorrector class."""
    
    @pytest.fixture
    def corrector(self):
        """Create a corrector instance."""
        return GCCorrector(n_bins=101)
    
    def test_init_defaults(self):
        """Test default initialization."""
        corrector = GCCorrector()
        assert corrector.n_bins == 101
        assert corrector.smoothing_span == 0.3
        assert corrector.min_count == 10
    
    def test_compute_gc_content(self, corrector):
        """Test GC content calculation."""
        assert corrector.compute_gc_content("GGCC") == 1.0
        assert corrector.compute_gc_content("AATT") == 0.0
        assert corrector.compute_gc_content("ACGT") == 0.5
        assert corrector.compute_gc_content("GCGCATAT") == 0.5
        assert corrector.compute_gc_content("") == 0.0
    
    def test_compute_bias_uniform(self, corrector):
        """Test bias computation with uniform GC distribution."""
        # Uniform GC values should have minimal bias
        np.random.seed(42)
        gc_values = np.random.uniform(0.3, 0.7, size=10000)
        
        profile = corrector.compute_bias_from_fragments(gc_values)
        
        assert isinstance(profile, GCBiasProfile)
        assert 0.3 < profile.mean_gc < 0.7
        assert profile.std_gc > 0
    
    def test_compute_bias_skewed(self, corrector):
        """Test bias computation with skewed GC distribution."""
        np.random.seed(42)
        # Heavily skewed toward high GC
        gc_values = np.random.beta(5, 2, size=10000)
        
        profile = corrector.compute_bias_from_fragments(gc_values)
        
        # Mean should be high
        assert profile.mean_gc > 0.5
    
    def test_correction_factors_shape(self, corrector):
        """Test that correction factors have correct shape."""
        np.random.seed(42)
        gc_values = np.random.uniform(0.2, 0.8, size=5000)
        
        profile = corrector.compute_bias_from_fragments(gc_values)
        
        assert len(profile.correction_factors) == corrector.n_bins
        assert len(profile.gc_bins) == corrector.n_bins
    
    def test_correct_values(self, corrector):
        """Test value correction."""
        np.random.seed(42)
        gc_values = np.random.uniform(0.2, 0.8, size=5000)
        
        profile = corrector.compute_bias_from_fragments(gc_values)
        
        # Correct some values
        values = np.ones(100)
        gc_for_correction = np.random.uniform(0.2, 0.8, size=100)
        
        corrected = corrector.correct_values(values, gc_for_correction, profile)
        
        assert len(corrected) == len(values)
        assert corrected.dtype == values.dtype


class TestGCBiasProfile:
    """Tests for GCBiasProfile dataclass."""
    
    @pytest.fixture
    def sample_profile(self):
        """Create sample profile."""
        n_bins = 101
        return GCBiasProfile(
            gc_bins=np.linspace(0, 1, n_bins),
            observed_counts=np.random.poisson(100, n_bins),
            expected_counts=np.full(n_bins, 100),
            correction_factors=np.ones(n_bins),
            r_squared=0.15,
            gc_dropout=0.05,
            mean_gc=0.42,
            std_gc=0.08,
        )
    
    def test_to_dict(self, sample_profile):
        """Test serialization."""
        d = sample_profile.to_dict()
        
        assert isinstance(d, dict)
        assert d['r_squared'] == 0.15
        assert d['mean_gc'] == 0.42
        assert d['n_bins'] == 101
    
    def test_summary(self, sample_profile):
        """Test summary generation."""
        summary = sample_profile.summary()
        
        assert isinstance(summary, str)
        assert 'GC Bias Profile' in summary
        assert '42.0%' in summary  # mean_gc
        assert '0.150' in summary  # r_squared
