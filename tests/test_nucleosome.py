"""
Tests for nucleosome positioning / WPS analysis module.
"""

import numpy as np
import pytest
from pathlib import Path

from fragmentomics.features.nucleosome import (
    NucleosomeAnalyzer,
    WPSProfile,
    compute_wps,
    NUCLEOSOME_SIZE,
    NUCLEOSOME_REPEAT,
)


# Test data directory
DATA_DIR = Path(__file__).parent / "data"


class TestNucleosomeConstants:
    """Test nucleosome-related constants."""

    def test_nucleosome_size(self):
        """Nucleosome wraps 147bp of DNA."""
        assert NUCLEOSOME_SIZE == 147

    def test_nucleosome_repeat(self):
        """Nucleosome repeat length is ~167bp."""
        assert NUCLEOSOME_REPEAT == 167


class TestWPSProfile:
    """Tests for WPSProfile dataclass."""

    @pytest.fixture
    def sample_profile(self):
        """Create sample WPS profile."""
        n = 1000
        return WPSProfile(
            positions=np.arange(1000000, 1000000 + n),
            wps=np.random.randn(n) * 10,
            smoothed_wps=np.random.randn(n) * 8,
            chrom="chr1",
            start=1000000,
            end=1000000 + n,
            mean_wps=0.5,
            std_wps=10.2,
            peak_positions=np.array([1000100, 1000267, 1000434, 1000601, 1000768]),
            trough_positions=np.array([1000183, 1000350, 1000517, 1000684]),
            periodicity=167.0,
        )

    def test_to_dict(self, sample_profile):
        """Test serialization."""
        d = sample_profile.to_dict()

        assert isinstance(d, dict)
        assert d["chrom"] == "chr1"
        assert d["start"] == 1000000
        assert d["end"] == 1001000
        assert d["mean_wps"] == 0.5
        assert d["n_peaks"] == 5
        assert d["n_troughs"] == 4
        assert d["periodicity"] == 167.0

    def test_summary(self, sample_profile):
        """Test summary generation."""
        summary = sample_profile.summary()

        assert isinstance(summary, str)
        assert "WPS Profile" in summary
        assert "chr1" in summary
        assert "1,000,000" in summary
        assert "Nucleosome peaks: 5" in summary
        assert "167.0" in summary


class TestNucleosomeAnalyzer:
    """Tests for NucleosomeAnalyzer class."""

    def test_init_defaults(self):
        """Test default initialization."""
        analyzer = NucleosomeAnalyzer()

        assert analyzer.window_size == 120
        assert analyzer.min_fragment_size == 120
        assert analyzer.max_fragment_size == 180
        assert analyzer.smoothing_bandwidth == 15
        assert analyzer.min_mapq == 30

    def test_init_custom(self):
        """Test custom initialization."""
        analyzer = NucleosomeAnalyzer(
            window_size=100,
            min_fragment_size=100,
            max_fragment_size=200,
            smoothing_bandwidth=20,
            min_mapq=20,
        )

        assert analyzer.window_size == 100
        assert analyzer.min_fragment_size == 100
        assert analyzer.max_fragment_size == 200
        assert analyzer.smoothing_bandwidth == 20
        assert analyzer.min_mapq == 20

    def test_compute_periodicity_empty(self):
        """Test periodicity with empty peaks."""
        analyzer = NucleosomeAnalyzer()

        # Empty array
        result = analyzer._compute_periodicity(np.array([]))
        assert result == 0.0

        # Single peak
        result = analyzer._compute_periodicity(np.array([100]))
        assert result == 0.0

        # Two peaks (not enough for reliable periodicity)
        result = analyzer._compute_periodicity(np.array([100, 267]))
        assert result == 0.0

    def test_compute_periodicity_regular(self):
        """Test periodicity with regular spacing."""
        analyzer = NucleosomeAnalyzer()

        # Regular ~167bp spacing
        peaks = np.array([1000, 1167, 1334, 1501, 1668])
        result = analyzer._compute_periodicity(peaks)

        assert 160 < result < 175  # Should be around 167

    def test_compute_periodicity_irregular(self):
        """Test periodicity with irregular spacing uses median."""
        analyzer = NucleosomeAnalyzer()

        # Irregular spacing: 167, 180, 165, 300 (outlier)
        peaks = np.array([1000, 1167, 1347, 1512, 1812])
        result = analyzer._compute_periodicity(peaks)

        # Median should be robust to outlier
        assert 160 < result < 190

    def test_find_peaks_troughs_synthetic(self):
        """Test peak/trough finding on synthetic nucleosome pattern."""
        analyzer = NucleosomeAnalyzer()

        # Create synthetic nucleosomal WPS signal
        # Peaks at nucleosome centers, troughs at linkers
        x = np.arange(2000)
        # Nucleosome repeat ~167bp
        wps = np.sin(2 * np.pi * x / 167) * 50 + np.random.randn(2000) * 5

        peaks, troughs = analyzer._find_peaks_troughs(wps, offset=0)

        # Should find approximately 12 nucleosomes in 2000bp
        assert 8 < len(peaks) < 16

        # Peaks and troughs should alternate
        if len(peaks) > 1 and len(troughs) > 0:
            # Check spacing is approximately correct
            peak_spacing = np.median(np.diff(peaks))
            assert 130 < peak_spacing < 200


class TestComputeWPSWithData:
    """Integration tests with actual BAM files."""

    @pytest.fixture
    def healthy_bam(self):
        """Path to healthy sample BAM."""
        return DATA_DIR / "healthy_sample.bam"

    def test_compute_wps_basic(self, healthy_bam, ensure_test_data):
        """Test basic WPS computation."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        analyzer = NucleosomeAnalyzer()
        profile = analyzer.compute_wps(
            healthy_bam,
            chrom="chr1",
            start=10000,
            end=50000,
        )

        assert isinstance(profile, WPSProfile)
        assert profile.chrom == "chr1"
        assert profile.start == 10000
        assert profile.end == 50000
        assert len(profile.positions) == 40000
        assert len(profile.wps) == 40000
        assert len(profile.smoothed_wps) == 40000

    def test_compute_wps_convenience_function(self, healthy_bam, ensure_test_data):
        """Test convenience function."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        profile = compute_wps(
            healthy_bam,
            chrom="chr1",
            start=10000,
            end=30000,
            window_size=100,
        )

        assert isinstance(profile, WPSProfile)
        assert len(profile.wps) == 20000

    def test_compute_coverage_profile(self, healthy_bam, ensure_test_data):
        """Test coverage profile computation."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        analyzer = NucleosomeAnalyzer()
        coverage = analyzer.compute_coverage_profile(
            healthy_bam,
            chrom="chr1",
            start=10000,
            end=50000,
        )

        assert isinstance(coverage, np.ndarray)
        assert len(coverage) == 40000
        assert coverage.dtype == np.int32
        # Should have some coverage
        assert coverage.sum() > 0


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_region(self):
        """Test with region that has no reads."""
        analyzer = NucleosomeAnalyzer()

        # This would need a real BAM with an empty region
        # For now, test the periodicity edge cases
        assert analyzer._compute_periodicity(np.array([])) == 0.0

    def test_very_small_region(self):
        """Test with very small region."""
        # Create minimal WPS array
        wps = np.zeros(100)
        analyzer = NucleosomeAnalyzer()

        peaks, troughs = analyzer._find_peaks_troughs(wps)

        # Should not crash, may return empty arrays
        assert isinstance(peaks, np.ndarray)
        assert isinstance(troughs, np.ndarray)
