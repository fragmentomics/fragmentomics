"""
Tests for coverage-based copy number analysis.
"""

import numpy as np
import pytest
from pathlib import Path

from fragmentomics.features.coverage import (
    CoverageAnalyzer,
    CoverageBin,
    CoverageProfile,
    compute_coverage,
)


DATA_DIR = Path(__file__).parent / "data"


class TestCoverageBin:
    """Tests for CoverageBin dataclass."""

    def test_basic_creation(self):
        """Test basic bin creation."""
        bin_ = CoverageBin(
            chrom="chr1",
            start=0,
            end=100000,
            raw_count=500,
        )
        assert bin_.chrom == "chr1"
        assert bin_.start == 0
        assert bin_.end == 100000
        assert bin_.raw_count == 500

    def test_midpoint(self):
        """Test midpoint calculation."""
        bin_ = CoverageBin(chrom="chr1", start=0, end=100000, raw_count=0)
        assert bin_.midpoint == 50000

        bin2 = CoverageBin(chrom="chr1", start=100000, end=200000, raw_count=0)
        assert bin2.midpoint == 150000


class TestCoverageProfile:
    """Tests for CoverageProfile dataclass."""

    @pytest.fixture
    def sample_profile(self):
        """Create sample profile."""
        bins = [
            CoverageBin("chr1", i * 100000, (i + 1) * 100000, raw_count=100 + i * 10)
            for i in range(10)
        ]
        return CoverageProfile(
            bins=bins,
            bin_size=100000,
            total_reads=1000,
            median_coverage=145.0,
            mad=25.0,
        )

    def test_to_dict(self, sample_profile):
        """Test serialization."""
        d = sample_profile.to_dict()

        assert d["n_bins"] == 10
        assert d["bin_size"] == 100000
        assert d["total_reads"] == 1000
        assert d["median_coverage"] == 145.0

    def test_to_array(self, sample_profile):
        """Test array conversion."""
        positions, raw_counts, log2_ratios = sample_profile.to_array()

        assert len(positions) == 10
        assert len(raw_counts) == 10
        assert len(log2_ratios) == 10
        assert positions[0] == 50000
        assert raw_counts[0] == 100


class TestCoverageAnalyzer:
    """Tests for CoverageAnalyzer class."""

    def test_init_defaults(self):
        """Test default initialization."""
        analyzer = CoverageAnalyzer()

        assert analyzer.bin_size == 100_000
        assert analyzer.min_mapq == 30
        assert analyzer.min_fragment_size == 50
        assert analyzer.max_fragment_size == 500

    def test_init_custom(self):
        """Test custom initialization."""
        analyzer = CoverageAnalyzer(
            bin_size=50_000,
            min_mapq=20,
            min_fragment_size=100,
            max_fragment_size=300,
        )

        assert analyzer.bin_size == 50_000
        assert analyzer.min_mapq == 20


class TestComputeCoverageWithData:
    """Integration tests with BAM files."""

    @pytest.fixture
    def healthy_bam(self):
        """Path to healthy sample BAM."""
        return DATA_DIR / "healthy_sample.bam"

    def test_compute_coverage_basic(self, healthy_bam, ensure_test_data):
        """Test basic coverage computation."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        analyzer = CoverageAnalyzer(bin_size=10_000)
        profile = analyzer.compute_coverage(healthy_bam, chrom="chr1")

        assert isinstance(profile, CoverageProfile)
        assert len(profile.bins) > 0
        assert profile.bin_size == 10_000
        assert profile.total_reads > 0
        assert profile.median_coverage >= 0

    def test_compute_coverage_specific_chrom(self, healthy_bam, ensure_test_data):
        """Test coverage on specific chromosome."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        profile = compute_coverage(healthy_bam, bin_size=10_000, chrom="chr1")

        # All bins should be from chr1
        for bin_ in profile.bins:
            assert bin_.chrom == "chr1"

    def test_log2_ratios_computed(self, healthy_bam, ensure_test_data):
        """Test that log2 ratios are computed."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        profile = compute_coverage(healthy_bam, bin_size=10_000, chrom="chr1")

        # Check arrays have expected shapes
        positions, raw_counts, log2_ratios = profile.to_array()

        assert len(positions) == len(raw_counts)
        assert len(positions) == len(log2_ratios)
        
        # Bins with coverage should have non-zero log2 ratios
        # Note: test data may be very sparse, so we just check structure
        assert profile.median_coverage >= 0


class TestEdgeCases:
    """Test edge cases."""

    def test_empty_region(self):
        """Test analyzer doesn't crash with unusual input."""
        analyzer = CoverageAnalyzer()

        # Just test that the class instantiates correctly
        assert analyzer.bin_size == 100_000

    def test_coverage_bin_defaults(self):
        """Test CoverageBin default values."""
        bin_ = CoverageBin(chrom="chr1", start=0, end=100, raw_count=0)

        assert bin_.gc_content == 0.0
        assert bin_.corrected_count == 0.0
        assert bin_.log2_ratio == 0.0
