"""Pytest configuration and shared fixtures."""

import random
import subprocess
from pathlib import Path

import numpy as np
import pytest


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session", autouse=True)
def ensure_test_data():
    """Ensure test BAM files exist before running tests."""
    healthy_bam = DATA_DIR / "healthy_sample.bam"
    cancer_bam = DATA_DIR / "cancer_sample.bam"
    
    if not healthy_bam.exists() or not cancer_bam.exists():
        generate_script = DATA_DIR / "generate_test_bam.py"
        if generate_script.exists():
            subprocess.run(
                ["python", str(generate_script)],
                check=True,
                capture_output=True,
            )


@pytest.fixture
def healthy_bam():
    """Path to healthy sample BAM file."""
    return DATA_DIR / "healthy_sample.bam"


@pytest.fixture
def cancer_bam():
    """Path to cancer sample BAM file."""
    return DATA_DIR / "cancer_sample.bam"


@pytest.fixture
def data_dir():
    """Path to test data directory."""
    return DATA_DIR


@pytest.fixture
def synthetic_sizes():
    """Generate synthetic healthy cfDNA fragment sizes.
    
    Simulates a healthy sample with strong nucleosomal pattern:
    - Primary peak at ~167bp (mononucleosome)
    - Secondary peak at ~334bp (dinucleosome)
    - Few short fragments
    """
    np.random.seed(42)
    
    # Mononucleosome peak (~80% of fragments)
    mono = np.random.normal(167, 15, 800)
    
    # Dinucleosome peak (~15% of fragments)
    di = np.random.normal(334, 20, 150)
    
    # Short fragments (~5% of fragments)
    short = np.random.normal(100, 20, 50)
    
    sizes = np.concatenate([mono, di, short])
    sizes = np.clip(sizes, 30, 600).astype(int)
    
    return sizes


@pytest.fixture
def synthetic_cancer_sizes():
    """Generate synthetic cancer cfDNA fragment sizes.
    
    Simulates a cancer sample with:
    - Weaker nucleosomal pattern
    - More short fragments (tumor-derived)
    - Shifted size distribution
    """
    np.random.seed(123)
    
    # Mononucleosome peak (weaker, ~50% of fragments)
    mono = np.random.normal(167, 25, 500)
    
    # Short fragments (elevated, ~35% of fragments)
    short = np.random.normal(100, 25, 350)
    
    # Very short tumor fragments (~15% of fragments)
    tumor = np.random.normal(70, 15, 150)
    
    sizes = np.concatenate([mono, short, tumor])
    sizes = np.clip(sizes, 30, 600).astype(int)
    
    return sizes
