"""
Tests for genomic region handling (BED files).
"""

import tempfile
from pathlib import Path

import pytest

from fragmentomics.io.regions import (
    GenomicRegion,
    read_bed,
    read_bed_to_list,
    write_bed,
    parse_region_string,
    merge_overlapping,
)


class TestGenomicRegion:
    """Tests for GenomicRegion dataclass."""

    def test_basic_creation(self):
        """Test basic region creation."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)

        assert region.chrom == "chr1"
        assert region.start == 1000
        assert region.end == 2000
        assert region.name == ""
        assert region.strand == "."

    def test_with_all_fields(self):
        """Test creation with all fields."""
        region = GenomicRegion(
            chrom="chr1",
            start=1000,
            end=2000,
            name="promoter1",
            score=100.5,
            strand="+",
        )

        assert region.name == "promoter1"
        assert region.score == 100.5
        assert region.strand == "+"

    def test_length(self):
        """Test length property."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)
        assert region.length == 1000

    def test_midpoint(self):
        """Test midpoint property."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)
        assert region.midpoint == 1500

    def test_str(self):
        """Test string representation."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)
        assert str(region) == "chr1:1000-2000"

        region_named = GenomicRegion(chrom="chr1", start=1000, end=2000, name="test")
        assert "test" in str(region_named)

    def test_to_ucsc(self):
        """Test UCSC coordinate format (1-based)."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)
        assert region.to_ucsc() == "chr1:1001-2000"

    def test_overlaps_true(self):
        """Test overlapping regions."""
        r1 = GenomicRegion(chrom="chr1", start=1000, end=2000)
        r2 = GenomicRegion(chrom="chr1", start=1500, end=2500)

        assert r1.overlaps(r2)
        assert r2.overlaps(r1)

    def test_overlaps_false_different_chrom(self):
        """Test non-overlapping regions on different chromosomes."""
        r1 = GenomicRegion(chrom="chr1", start=1000, end=2000)
        r2 = GenomicRegion(chrom="chr2", start=1000, end=2000)

        assert not r1.overlaps(r2)

    def test_overlaps_false_adjacent(self):
        """Test adjacent but non-overlapping regions."""
        r1 = GenomicRegion(chrom="chr1", start=1000, end=2000)
        r2 = GenomicRegion(chrom="chr1", start=2000, end=3000)

        assert not r1.overlaps(r2)

    def test_contains(self):
        """Test position containment."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000)

        assert region.contains("chr1", 1000)
        assert region.contains("chr1", 1500)
        assert region.contains("chr1", 1999)
        assert not region.contains("chr1", 2000)  # End is exclusive
        assert not region.contains("chr1", 999)
        assert not region.contains("chr2", 1500)

    def test_expand_plus_strand(self):
        """Test expansion on plus strand."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000, strand="+")
        expanded = region.expand(upstream=100, downstream=200)

        assert expanded.start == 900
        assert expanded.end == 2200

    def test_expand_minus_strand(self):
        """Test expansion on minus strand (swapped directions)."""
        region = GenomicRegion(chrom="chr1", start=1000, end=2000, strand="-")
        expanded = region.expand(upstream=100, downstream=200)

        assert expanded.start == 800  # downstream is towards lower coords
        assert expanded.end == 2100  # upstream is towards higher coords

    def test_expand_no_negative_start(self):
        """Test that expansion doesn't create negative start."""
        region = GenomicRegion(chrom="chr1", start=50, end=100)
        expanded = region.expand(upstream=100)

        assert expanded.start == 0


class TestParseRegionString:
    """Tests for parse_region_string function."""

    def test_basic_format(self):
        """Test basic chr:start-end format."""
        region = parse_region_string("chr1:1000-2000")

        assert region.chrom == "chr1"
        assert region.start == 1000
        assert region.end == 2000

    def test_with_commas(self):
        """Test format with commas in numbers."""
        region = parse_region_string("chr1:1,000,000-2,000,000")

        assert region.start == 1000000
        assert region.end == 2000000

    def test_with_dots(self):
        """Test format with .. separator."""
        region = parse_region_string("chr1:1000..2000")

        assert region.start == 1000
        assert region.end == 2000

    def test_invalid_format(self):
        """Test invalid format raises error."""
        with pytest.raises(ValueError):
            parse_region_string("invalid")

        with pytest.raises(ValueError):
            parse_region_string("chr1:1000")


class TestReadWriteBed:
    """Tests for BED file I/O."""

    @pytest.fixture
    def sample_bed_content(self):
        """Sample BED content."""
        return """chr1\t1000\t2000\tgene1\t100\t+
chr1\t3000\t4000\tgene2\t200\t-
chr2\t5000\t6000\tgene3\t300\t.
"""

    def test_read_bed_basic(self, sample_bed_content):
        """Test reading BED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(sample_bed_content)
            path = f.name

        try:
            regions = list(read_bed(path))

            assert len(regions) == 3
            assert regions[0].chrom == "chr1"
            assert regions[0].start == 1000
            assert regions[0].end == 2000
            assert regions[0].name == "gene1"
            assert regions[0].score == 100.0
            assert regions[0].strand == "+"
        finally:
            Path(path).unlink()

    def test_read_bed_with_header(self):
        """Test reading BED with header lines."""
        content = """# Header comment
track name=test
chr1\t1000\t2000
chr1\t3000\t4000
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            path = f.name

        try:
            regions = list(read_bed(path, skip_header=True))
            assert len(regions) == 2
        finally:
            Path(path).unlink()

    def test_read_bed3(self):
        """Test reading BED3 format."""
        content = "chr1\t1000\t2000\nchr1\t3000\t4000\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            path = f.name

        try:
            regions = list(read_bed(path))

            assert len(regions) == 2
            assert regions[0].name == ""
            assert regions[0].strand == "."
        finally:
            Path(path).unlink()

    def test_read_bed_to_list(self, sample_bed_content):
        """Test read_bed_to_list convenience function."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(sample_bed_content)
            path = f.name

        try:
            regions = read_bed_to_list(path)
            assert isinstance(regions, list)
            assert len(regions) == 3
        finally:
            Path(path).unlink()

    def test_write_bed6(self):
        """Test writing BED6 format."""
        regions = [
            GenomicRegion("chr1", 1000, 2000, "gene1", 100, "+"),
            GenomicRegion("chr1", 3000, 4000, "gene2", 200, "-"),
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            path = f.name

        try:
            write_bed(regions, path, bed_format=6)

            # Read back and verify
            read_back = read_bed_to_list(path)
            assert len(read_back) == 2
            assert read_back[0].name == "gene1"
            assert read_back[0].strand == "+"
        finally:
            Path(path).unlink()

    def test_write_bed3(self):
        """Test writing BED3 format."""
        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 3000, 4000),
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            path = f.name

        try:
            write_bed(regions, path, bed_format=3)

            with open(path) as f:
                lines = f.readlines()

            assert len(lines) == 2
            assert lines[0].strip() == "chr1\t1000\t2000"
        finally:
            Path(path).unlink()


class TestMergeOverlapping:
    """Tests for merge_overlapping function."""

    def test_merge_overlapping_regions(self):
        """Test merging overlapping regions."""
        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 1500, 2500),
            GenomicRegion("chr1", 2400, 3000),
        ]

        merged = merge_overlapping(regions)

        assert len(merged) == 1
        assert merged[0].start == 1000
        assert merged[0].end == 3000

    def test_merge_separate_regions(self):
        """Test that non-overlapping regions stay separate."""
        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 3000, 4000),
        ]

        merged = merge_overlapping(regions)

        assert len(merged) == 2

    def test_merge_different_chroms(self):
        """Test that different chromosomes stay separate."""
        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr2", 1000, 2000),
        ]

        merged = merge_overlapping(regions)

        assert len(merged) == 2

    def test_merge_with_gap(self):
        """Test merging with gap tolerance."""
        regions = [
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 2010, 3000),  # 10bp gap
        ]

        # Without gap
        merged = merge_overlapping(regions, gap=0)
        assert len(merged) == 2

        # With gap
        merged = merge_overlapping(regions, gap=20)
        assert len(merged) == 1

    def test_merge_empty_list(self):
        """Test merging empty list."""
        assert merge_overlapping([]) == []

    def test_merge_unsorted(self):
        """Test that unsorted input is handled."""
        regions = [
            GenomicRegion("chr1", 3000, 4000),
            GenomicRegion("chr1", 1000, 2000),
            GenomicRegion("chr1", 1500, 2500),
        ]

        merged = merge_overlapping(regions)

        # Should be sorted and merged properly
        assert merged[0].start == 1000
