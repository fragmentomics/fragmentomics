"""
Tests for DELFI-style fragmentation analysis.
"""

import numpy as np
import pytest
from pathlib import Path

from fragmentomics.features.delfi import (
    DELFIAnalyzer,
    DELFIBin,
    DELFIProfile,
    compute_delfi_profile,
    create_healthy_reference,
    DEFAULT_BIN_SIZE,
    SHORT_MIN,
    SHORT_MAX,
    LONG_MIN,
    LONG_MAX,
)


DATA_DIR = Path(__file__).parent / "data"


class TestDELFIConstants:
    """Test DELFI constants match the paper."""

    def test_bin_size(self):
        """Default bin size is 5Mb per Cristiano et al."""
        assert DEFAULT_BIN_SIZE == 5_000_000

    def test_short_range(self):
        """Short fragments: 100-150bp."""
        assert SHORT_MIN == 100
        assert SHORT_MAX == 150

    def test_long_range(self):
        """Long fragments: 151-220bp."""
        assert LONG_MIN == 151
        assert LONG_MAX == 220


class TestDELFIBin:
    """Tests for DELFIBin dataclass."""

    def test_ratio_calculation(self):
        """Test short/long ratio."""
        bin_ = DELFIBin(
            chrom="chr1",
            start=0,
            end=5_000_000,
            short_count=100,
            long_count=200,
            total_count=300,
        )
        assert bin_.ratio == 0.5

    def test_ratio_zero_long(self):
        """Test ratio when no long fragments."""
        bin_ = DELFIBin(
            chrom="chr1",
            start=0,
            end=5_000_000,
            short_count=100,
            long_count=0,
            total_count=100,
        )
        assert bin_.ratio == 0.0

    def test_coverage(self):
        """Test coverage calculation."""
        bin_ = DELFIBin(
            chrom="chr1",
            start=0,
            end=5_000_000,
            short_count=0,
            long_count=0,
            total_count=5000,
        )
        # 5000 fragments / 5 Mb = 1000 per Mb
        assert bin_.coverage == 1000.0


class TestDELFIProfile:
    """Tests for DELFIProfile dataclass."""

    @pytest.fixture
    def sample_profile(self):
        """Create sample profile."""
        bins = []
        for i in range(10):
            bins.append(DELFIBin(
                chrom="chr1",
                start=i * 5_000_000,
                end=(i + 1) * 5_000_000,
                short_count=100 + i * 5,
                long_count=200 - i * 3,
                total_count=300,
            ))
        return DELFIProfile(
            bins=bins,
            total_fragments=3000,
            total_short=1000,
            total_long=1700,
            genome_wide_ratio=1000 / 1700,
            sample_name="test_sample",
        )

    def test_to_feature_vector(self, sample_profile):
        """Test feature vector generation."""
        features = sample_profile.to_feature_vector()

        # Should have ratios + coverages
        assert len(features) == 20  # 10 ratios + 10 coverages

    def test_to_ratio_vector(self, sample_profile):
        """Test ratio vector."""
        ratios = sample_profile.to_ratio_vector()

        assert len(ratios) == 10
        # First ratio: 100/200 = 0.5
        assert abs(ratios[0] - 0.5) < 0.01

    def test_to_dict(self, sample_profile):
        """Test serialization."""
        d = sample_profile.to_dict()

        assert d["sample_name"] == "test_sample"
        assert d["n_bins"] == 10
        assert d["total_fragments"] == 3000

    def test_summary(self, sample_profile):
        """Test summary generation."""
        summary = sample_profile.summary()

        assert "DELFI Profile" in summary
        assert "test_sample" in summary
        assert "3,000" in summary


class TestDELFIAnalyzer:
    """Tests for DELFIAnalyzer class."""

    def test_init_defaults(self):
        """Test default initialization."""
        analyzer = DELFIAnalyzer()

        assert analyzer.bin_size == 5_000_000
        assert analyzer.short_min == 100
        assert analyzer.short_max == 150
        assert analyzer.long_min == 151
        assert analyzer.long_max == 220
        assert analyzer.min_mapq == 30

    def test_init_custom(self):
        """Test custom initialization."""
        analyzer = DELFIAnalyzer(
            bin_size=1_000_000,
            short_range=(80, 140),
            long_range=(141, 200),
            min_mapq=20,
        )

        assert analyzer.bin_size == 1_000_000
        assert analyzer.short_min == 80
        assert analyzer.short_max == 140

    def test_exclude_chroms(self):
        """Test chromosome exclusion."""
        analyzer = DELFIAnalyzer()

        assert analyzer._should_exclude("chrM")
        assert analyzer._should_exclude("chrY")
        assert analyzer._should_exclude("chr1_random")
        assert not analyzer._should_exclude("chr1")
        assert not analyzer._should_exclude("chr22")


class TestDELFIWithData:
    """Integration tests with BAM files."""

    @pytest.fixture
    def healthy_bam(self):
        """Path to healthy sample BAM."""
        return DATA_DIR / "healthy_sample.bam"

    def test_analyze_basic(self, healthy_bam, ensure_test_data):
        """Test basic DELFI analysis."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        # Use smaller bin size for test data
        analyzer = DELFIAnalyzer(bin_size=100_000)
        profile = analyzer.analyze(healthy_bam)

        assert isinstance(profile, DELFIProfile)
        assert len(profile.bins) > 0
        assert profile.total_fragments > 0

    def test_convenience_function(self, healthy_bam, ensure_test_data):
        """Test compute_delfi_profile function."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        profile = compute_delfi_profile(healthy_bam, bin_size=100_000)

        assert isinstance(profile, DELFIProfile)


class TestCorrelation:
    """Test correlation computation."""

    def test_compute_correlation_identical(self):
        """Identical profiles should have correlation 1.0."""
        # Use varying ratios so correlation is meaningful
        bins = [
            DELFIBin("chr1", i * 5_000_000, (i + 1) * 5_000_000, 100 + i * 10, 200, 300)
            for i in range(10)
        ]
        profile = DELFIProfile(bins=bins, total_fragments=3000)
        reference = profile.to_ratio_vector()

        corr = DELFIAnalyzer.compute_correlation(profile, reference)
        assert abs(corr - 1.0) < 0.001

    def test_compute_correlation_different(self):
        """Different profiles should have lower correlation."""
        # Profile 1: increasing ratios
        bins1 = [
            DELFIBin("chr1", i * 5_000_000, (i + 1) * 5_000_000, 100 + i * 10, 200, 300)
            for i in range(10)
        ]
        profile1 = DELFIProfile(bins=bins1, total_fragments=3000)

        # Profile 2: decreasing ratios
        bins2 = [
            DELFIBin("chr1", i * 5_000_000, (i + 1) * 5_000_000, 200 - i * 10, 200, 300)
            for i in range(10)
        ]
        profile2 = DELFIProfile(bins=bins2, total_fragments=3000)

        corr = DELFIAnalyzer.compute_correlation(profile1, profile2.to_ratio_vector())
        # Should be negative correlation
        assert corr < 0


class TestHealthyReference:
    """Test healthy reference creation."""

    def test_create_reference(self):
        """Test creating healthy reference."""
        # Create multiple profiles
        profiles = []
        for seed in range(5):
            np.random.seed(seed)
            bins = [
                DELFIBin(
                    "chr1",
                    i * 5_000_000,
                    (i + 1) * 5_000_000,
                    int(100 + np.random.normal(0, 5)),
                    200,
                    300,
                )
                for i in range(10)
            ]
            profiles.append(DELFIProfile(bins=bins, total_fragments=3000))

        reference = create_healthy_reference(profiles)

        assert len(reference) == 10
        # Reference should be median of individual ratios
        assert all(0 < r < 1 for r in reference)

    def test_create_reference_empty_raises(self):
        """Empty profiles list should raise."""
        with pytest.raises(ValueError):
            create_healthy_reference([])
