"""
Tests for fragment size analysis module.
"""

import numpy as np
import pytest

from fragmentomics.features.sizes import (
    FragmentSizeAnalyzer,
    SizeDistribution,
    analyze_sizes,
)


class TestFragmentSizeAnalyzer:
    """Tests for FragmentSizeAnalyzer class."""
    
    def test_basic_analysis(self, synthetic_sizes):
        """Test basic size analysis on synthetic data."""
        analyzer = FragmentSizeAnalyzer()
        dist = analyzer.analyze(synthetic_sizes)
        
        assert isinstance(dist, SizeDistribution)
        assert dist.n_fragments > 0
        assert dist.mean > 0
        assert dist.median > 0
        assert dist.mode > 0
    
    def test_healthy_characteristics(self, synthetic_sizes):
        """Test that healthy sample has expected characteristics."""
        dist = analyze_sizes(synthetic_sizes)
        
        # Mononucleosome peak should be around 160-175bp
        assert 160 <= dist.peak_mono <= 175
        
        # Median should be in mononucleosome range
        assert 150 <= dist.median <= 180
        
        # Short fragment ratio should be low in healthy
        assert dist.ratio_short < 0.20
        
        # Mononucleosome ratio should be high
        assert dist.ratio_mono > 0.30
    
    def test_cancer_characteristics(self, synthetic_cancer_sizes):
        """Test that cancer sample shows elevated short fragments."""
        dist = analyze_sizes(synthetic_cancer_sizes)
        
        # Cancer samples have more short fragments
        assert dist.ratio_short > 0.15
        
        # Mode may be shifted
        assert dist.mode > 0
    
    def test_short_ratio_difference(self, synthetic_sizes, synthetic_cancer_sizes):
        """Test that cancer has higher short fragment ratio than healthy."""
        healthy_dist = analyze_sizes(synthetic_sizes)
        cancer_dist = analyze_sizes(synthetic_cancer_sizes)
        
        assert cancer_dist.ratio_short > healthy_dist.ratio_short
    
    def test_empty_input_raises(self):
        """Test that empty input raises ValueError."""
        analyzer = FragmentSizeAnalyzer()
        with pytest.raises(ValueError, match="No fragments"):
            analyzer.analyze(np.array([], dtype=np.int32))
    
    def test_out_of_range_input(self):
        """Test handling of all sizes outside analysis range."""
        analyzer = FragmentSizeAnalyzer(min_size=100, max_size=200)
        sizes = np.array([50, 60, 70, 250, 300], dtype=np.int32)
        
        with pytest.raises(ValueError, match="No fragments in range"):
            analyzer.analyze(sizes)
    
    def test_custom_parameters(self, synthetic_sizes):
        """Test analyzer with custom parameters."""
        analyzer = FragmentSizeAnalyzer(
            bin_size=5,
            min_size=100,
            max_size=400,
            smooth_window=3,
        )
        dist = analyzer.analyze(synthetic_sizes)
        
        assert dist.min_size >= 100
        assert dist.max_size <= 400
    
    def test_to_dict(self, synthetic_sizes):
        """Test serialization to dictionary."""
        dist = analyze_sizes(synthetic_sizes)
        d = dist.to_dict()
        
        assert isinstance(d, dict)
        assert "n_fragments" in d
        assert "mean" in d
        assert "median" in d
        assert "ratio_short" in d
        assert "peak_mono" in d
    
    def test_summary(self, synthetic_sizes):
        """Test summary string generation."""
        dist = analyze_sizes(synthetic_sizes)
        summary = dist.summary()
        
        assert isinstance(summary, str)
        assert "Fragment Size Distribution Summary" in summary
        assert "Median" in summary
        assert "Short fragments" in summary


class TestAnalyzeSizes:
    """Tests for the analyze_sizes convenience function."""
    
    def test_basic_usage(self, synthetic_sizes):
        """Test basic function usage."""
        dist = analyze_sizes(synthetic_sizes)
        assert dist.n_fragments > 0
    
    def test_with_parameters(self, synthetic_sizes):
        """Test with custom parameters."""
        dist = analyze_sizes(
            synthetic_sizes,
            bin_size=2,
            min_size=100,
            max_size=400,
        )
        assert dist.min_size >= 100


class TestSizeDistribution:
    """Tests for SizeDistribution dataclass."""
    
    def test_attributes(self, synthetic_sizes):
        """Test all expected attributes are present."""
        dist = analyze_sizes(synthetic_sizes)
        
        # Raw data
        assert hasattr(dist, "sizes")
        assert hasattr(dist, "counts")
        assert hasattr(dist, "bin_centers")
        
        # Statistics
        assert hasattr(dist, "mean")
        assert hasattr(dist, "median")
        assert hasattr(dist, "mode")
        assert hasattr(dist, "std")
        
        # Features
        assert hasattr(dist, "ratio_short")
        assert hasattr(dist, "ratio_mono")
        assert hasattr(dist, "ratio_di")
        assert hasattr(dist, "peak_mono")
        assert hasattr(dist, "periodicity_10bp")
    
    def test_ratios_sum_reasonable(self, synthetic_sizes):
        """Test that size ratios are reasonable (not summing to > 1)."""
        dist = analyze_sizes(synthetic_sizes)
        
        # Each ratio should be between 0 and 1
        assert 0 <= dist.ratio_short <= 1
        assert 0 <= dist.ratio_mono <= 1
        assert 0 <= dist.ratio_di <= 1
        
        # Note: these ranges overlap, so sum can exceed 1
    
    def test_periodicity_bounded(self, synthetic_sizes):
        """Test that periodicity is bounded."""
        dist = analyze_sizes(synthetic_sizes)
        
        # Periodicity from autocorrelation should be bounded
        assert -1 <= dist.periodicity_10bp <= 1
