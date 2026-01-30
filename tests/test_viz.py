"""
Tests for visualization module.
"""

import numpy as np
import pytest
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing
import matplotlib.pyplot as plt

from fragmentomics.features.sizes import analyze_sizes
from fragmentomics.viz.plots import (
    plot_size_distribution,
    plot_size_comparison,
)


class TestPlotSizeDistribution:
    """Tests for plot_size_distribution function."""
    
    def test_basic_plot(self, synthetic_sizes):
        """Test basic plotting works."""
        dist = analyze_sizes(synthetic_sizes)
        fig, ax = plot_size_distribution(dist)
        
        assert fig is not None
        assert ax is not None
        plt.close(fig)
    
    def test_with_title(self, synthetic_sizes):
        """Test plotting with custom title."""
        dist = analyze_sizes(synthetic_sizes)
        fig, ax = plot_size_distribution(dist, title="Test Sample")
        
        assert ax.get_title() == "Test Sample"
        plt.close(fig)
    
    def test_without_features(self, synthetic_sizes):
        """Test plotting without feature annotations."""
        dist = analyze_sizes(synthetic_sizes)
        fig, ax = plot_size_distribution(dist, show_features=False)
        
        assert fig is not None
        plt.close(fig)
    
    def test_custom_color(self, synthetic_sizes):
        """Test plotting with custom color."""
        dist = analyze_sizes(synthetic_sizes)
        fig, ax = plot_size_distribution(dist, color="red")
        
        assert fig is not None
        plt.close(fig)
    
    def test_custom_figsize(self, synthetic_sizes):
        """Test plotting with custom figure size."""
        dist = analyze_sizes(synthetic_sizes)
        fig, ax = plot_size_distribution(dist, figsize=(8, 4))
        
        assert fig.get_figwidth() == 8
        assert fig.get_figheight() == 4
        plt.close(fig)
    
    def test_existing_axes(self, synthetic_sizes):
        """Test plotting on existing axes."""
        dist = analyze_sizes(synthetic_sizes)
        
        fig, ax = plt.subplots()
        fig2, ax2 = plot_size_distribution(dist, ax=ax)
        
        assert ax2 is ax
        plt.close(fig)


class TestPlotSizeComparison:
    """Tests for plot_size_comparison function."""
    
    def test_compare_two(self, synthetic_sizes, synthetic_cancer_sizes):
        """Test comparing two distributions."""
        dist1 = analyze_sizes(synthetic_sizes)
        dist2 = analyze_sizes(synthetic_cancer_sizes)
        
        fig, ax = plot_size_comparison(
            [dist1, dist2],
            ["Healthy", "Cancer"],
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_mismatched_labels_raises(self, synthetic_sizes):
        """Test that mismatched labels raises error."""
        dist = analyze_sizes(synthetic_sizes)
        
        with pytest.raises(ValueError, match="must match"):
            plot_size_comparison([dist], ["A", "B"])
    
    def test_custom_colors(self, synthetic_sizes, synthetic_cancer_sizes):
        """Test with custom colors."""
        dist1 = analyze_sizes(synthetic_sizes)
        dist2 = analyze_sizes(synthetic_cancer_sizes)
        
        fig, ax = plot_size_comparison(
            [dist1, dist2],
            ["A", "B"],
            colors=["blue", "red"],
        )
        
        assert fig is not None
        plt.close(fig)
