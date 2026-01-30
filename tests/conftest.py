"""
Pytest configuration and fixtures for fragmentomics tests.
"""

import numpy as np
import pytest
from pathlib import Path


@pytest.fixture
def synthetic_sizes():
    """
    Generate synthetic fragment sizes mimicking healthy cfDNA.
    
    Returns array with:
    - Mononucleosome peak around 167bp
    - Dinucleosome shoulder around 334bp
    - Typical cfDNA distribution characteristics
    """
    np.random.seed(42)
    
    # Mononucleosome peak (majority of fragments)
    mono = np.random.normal(167, 20, size=8000)
    
    # Dinucleosome shoulder
    di = np.random.normal(334, 30, size=1500)
    
    # Short fragments (small proportion in healthy)
    short = np.random.normal(120, 15, size=300)
    
    # Some noise
    noise = np.random.uniform(50, 500, size=200)
    
    # Combine and filter
    sizes = np.concatenate([mono, di, short, noise])
    sizes = sizes[(sizes >= 50) & (sizes <= 500)]
    sizes = sizes.astype(np.int32)
    
    return sizes


@pytest.fixture
def synthetic_cancer_sizes():
    """
    Generate synthetic fragment sizes mimicking cancer cfDNA.
    
    Cancer samples typically have:
    - Higher proportion of short fragments (<150bp)
    - More heterogeneous distribution
    - Shifted mononucleosome peak
    """
    np.random.seed(43)
    
    # Mononucleosome peak (slightly shifted/broader)
    mono = np.random.normal(165, 25, size=6000)
    
    # Dinucleosome
    di = np.random.normal(330, 35, size=1000)
    
    # Elevated short fragments (tumor signature)
    short = np.random.normal(115, 20, size=2500)
    
    # More noise
    noise = np.random.uniform(50, 500, size=500)
    
    sizes = np.concatenate([mono, di, short, noise])
    sizes = sizes[(sizes >= 50) & (sizes <= 500)]
    sizes = sizes.astype(np.int32)
    
    return sizes


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Create a temporary output directory."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir
