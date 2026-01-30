"""Integration tests for the fragmentomics CLI."""

import json
import subprocess
import tempfile
from pathlib import Path

import pytest


# Path to test data
DATA_DIR = Path(__file__).parent / "data"
HEALTHY_BAM = DATA_DIR / "healthy_sample.bam"
CANCER_BAM = DATA_DIR / "cancer_sample.bam"


@pytest.fixture(scope="module")
def ensure_test_data():
    """Ensure test BAM files exist."""
    if not HEALTHY_BAM.exists() or not CANCER_BAM.exists():
        # Generate test data
        subprocess.run(
            ["python", str(DATA_DIR / "generate_test_bam.py")],
            check=True,
        )
    assert HEALTHY_BAM.exists(), "healthy_sample.bam not found"
    assert CANCER_BAM.exists(), "cancer_sample.bam not found"


class TestCLIHelp:
    """Test CLI help and version commands."""

    def test_help(self):
        """Test main help command."""
        result = subprocess.run(
            ["fragmentomics", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "FragMentor" in result.stdout
        assert "sizes" in result.stdout
        assert "motifs" in result.stdout

    def test_version(self):
        """Test version flag."""
        result = subprocess.run(
            ["fragmentomics", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "0." in result.stdout  # Version starts with 0.

    def test_sizes_help(self):
        """Test sizes subcommand help."""
        result = subprocess.run(
            ["fragmentomics", "sizes", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "BAM" in result.stdout or "bam" in result.stdout.lower()


class TestSizesCommand:
    """Test the sizes analysis command."""

    def test_sizes_basic(self, ensure_test_data):
        """Test basic size analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                ["fragmentomics", "sizes", str(HEALTHY_BAM), "-o", tmpdir],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            
            # Check outputs exist
            outdir = Path(tmpdir)
            assert (outdir / "size_stats.json").exists()
            assert (outdir / "size_distribution.png").exists()
            
            # Verify JSON content
            with open(outdir / "size_stats.json") as f:
                stats = json.load(f)
            
            assert "n_fragments" in stats
            assert stats["n_fragments"] > 0
            assert "mean" in stats
            assert "ratio_short" in stats

    def test_sizes_healthy_vs_cancer(self, ensure_test_data):
        """Test that healthy and cancer samples show expected differences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Analyze healthy
            healthy_out = Path(tmpdir) / "healthy"
            healthy_out.mkdir()
            subprocess.run(
                ["fragmentomics", "sizes", str(HEALTHY_BAM), "-o", str(healthy_out)],
                check=True,
                capture_output=True,
            )
            
            # Analyze cancer
            cancer_out = Path(tmpdir) / "cancer"
            cancer_out.mkdir()
            subprocess.run(
                ["fragmentomics", "sizes", str(CANCER_BAM), "-o", str(cancer_out)],
                check=True,
                capture_output=True,
            )
            
            # Load results
            with open(healthy_out / "size_stats.json") as f:
                healthy_stats = json.load(f)
            with open(cancer_out / "size_stats.json") as f:
                cancer_stats = json.load(f)
            
            # Cancer should have more short fragments
            assert cancer_stats["ratio_short"] > healthy_stats["ratio_short"]
            # Cancer should have lower mean fragment size
            assert cancer_stats["mean"] < healthy_stats["mean"]

    def test_sizes_region_filter(self, ensure_test_data):
        """Test region filtering."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [
                    "fragmentomics", "sizes", str(HEALTHY_BAM),
                    "-o", tmpdir,
                    "--region", "chr1:10000-100000",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

    def test_sizes_json_output(self, ensure_test_data):
        """Test JSON-only output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [
                    "fragmentomics", "sizes", str(HEALTHY_BAM),
                    "-o", tmpdir,
                    "--no-plot",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            
            outdir = Path(tmpdir)
            assert (outdir / "size_stats.json").exists()
            # Plot should not exist with --no-plot
            assert not (outdir / "size_distribution.png").exists()

    def test_sizes_invalid_bam(self):
        """Test error handling for invalid BAM file."""
        result = subprocess.run(
            ["fragmentomics", "sizes", "/nonexistent/file.bam", "-o", "/tmp"],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0


class TestInfoCommand:
    """Test the info command."""

    def test_info_basic(self, ensure_test_data):
        """Test basic BAM info."""
        result = subprocess.run(
            ["fragmentomics", "info", str(HEALTHY_BAM)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "chr1" in result.stdout or "reads" in result.stdout.lower()


class TestMotifsCommand:
    """Test the motifs analysis command."""

    def test_motifs_help(self):
        """Test motifs help."""
        result = subprocess.run(
            ["fragmentomics", "motifs", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "motif" in result.stdout.lower()


class TestExtractCommand:
    """Test the extract command."""

    def test_extract_help(self):
        """Test extract help."""
        result = subprocess.run(
            ["fragmentomics", "extract", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
