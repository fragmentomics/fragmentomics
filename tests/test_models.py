"""
Tests for cancer detection model.
"""

import pytest
from pathlib import Path

from fragmentomics.models.cancer_detector import (
    CancerDetector,
    PredictionResult,
    predict_cancer,
    load_model,
)


DATA_DIR = Path(__file__).parent / "data"


class TestPredictionResult:
    """Tests for PredictionResult dataclass."""

    def test_to_dict(self):
        """Test serialization."""
        result = PredictionResult(
            probability=0.75,
            prediction="cancer",
            confidence="high",
            features_used=["ratio_short", "mean_size"],
            feature_importance={"ratio_short": 0.6, "mean_size": 0.4},
        )
        d = result.to_dict()

        assert d["probability"] == 0.75
        assert d["prediction"] == "cancer"
        assert d["confidence"] == "high"

    def test_summary(self):
        """Test summary generation."""
        result = PredictionResult(
            probability=0.75,
            prediction="cancer",
            confidence="high",
            features_used=[],
            feature_importance={},
        )
        summary = result.summary()

        assert "CANCER" in summary
        assert "75.0%" in summary


class TestCancerDetector:
    """Tests for CancerDetector class."""

    def test_init_default(self):
        """Test default initialization uses heuristic model."""
        detector = CancerDetector()

        assert detector.model is not None
        assert detector.model["type"] == "heuristic"
        assert detector.threshold == 0.5

    def test_init_custom_threshold(self):
        """Test custom threshold."""
        detector = CancerDetector(threshold=0.7)
        assert detector.threshold == 0.7

    def test_predict_healthy_features(self):
        """Test prediction with healthy-like features."""
        detector = CancerDetector()

        # Healthy-like features (high mono ratio, low short ratio)
        features = {
            "ratio_short": 0.10,
            "ratio_mono": 0.70,
            "mean_size": 168,
            "std_size": 35,
            "periodicity_10bp": 0.85,
            "genome_ratio": 0.35,
        }

        result = detector.predict(features)

        assert result.probability < 0.5
        assert result.prediction == "healthy"

    def test_predict_cancer_features(self):
        """Test prediction with cancer-like features."""
        detector = CancerDetector()

        # Cancer-like features (high short ratio, low mean size)
        features = {
            "ratio_short": 0.35,
            "ratio_mono": 0.45,
            "mean_size": 140,
            "std_size": 50,
            "periodicity_10bp": 0.5,
            "genome_ratio": 0.7,
        }

        result = detector.predict(features)

        assert result.probability > 0.5
        assert result.prediction == "cancer"

    def test_feature_importance_sums(self):
        """Test that feature importances sum to ~1."""
        detector = CancerDetector()

        features = {
            "ratio_short": 0.20,
            "ratio_mono": 0.60,
            "mean_size": 160,
            "std_size": 40,
            "periodicity_10bp": 0.7,
            "genome_ratio": 0.5,
        }

        result = detector.predict(features)
        total_importance = sum(result.feature_importance.values())

        assert 0.99 < total_importance < 1.01


class TestPredictFromBAM:
    """Integration tests with BAM files."""

    @pytest.fixture
    def healthy_bam(self):
        """Path to healthy sample BAM."""
        return DATA_DIR / "healthy_sample.bam"

    @pytest.fixture
    def cancer_bam(self):
        """Path to cancer sample BAM."""
        return DATA_DIR / "cancer_sample.bam"

    def test_predict_healthy_bam(self, healthy_bam, ensure_test_data):
        """Test prediction on healthy sample."""
        if not healthy_bam.exists():
            pytest.skip("Test BAM not available")

        detector = CancerDetector()
        result = detector.predict_from_bam(healthy_bam)

        assert isinstance(result, PredictionResult)
        assert 0 <= result.probability <= 1
        # Healthy sample should have lower probability
        assert result.probability < 0.6

    def test_predict_cancer_bam(self, cancer_bam, ensure_test_data):
        """Test prediction on cancer sample."""
        if not cancer_bam.exists():
            pytest.skip("Test BAM not available")

        detector = CancerDetector()
        result = detector.predict_from_bam(cancer_bam)

        assert isinstance(result, PredictionResult)
        # Cancer sample should have higher probability
        assert result.probability > 0.5

    def test_healthy_vs_cancer_difference(self, healthy_bam, cancer_bam, ensure_test_data):
        """Test that cancer sample has higher probability than healthy."""
        if not healthy_bam.exists() or not cancer_bam.exists():
            pytest.skip("Test BAMs not available")

        detector = CancerDetector()
        healthy_result = detector.predict_from_bam(healthy_bam)
        cancer_result = detector.predict_from_bam(cancer_bam)

        # Cancer should have higher probability
        assert cancer_result.probability > healthy_result.probability


class TestConvenienceFunctions:
    """Test module-level convenience functions."""

    def test_load_model_default(self):
        """Test loading default model."""
        detector = load_model()

        assert isinstance(detector, CancerDetector)
        assert detector.model is not None

    def test_predict_cancer_function(self, ensure_test_data):
        """Test predict_cancer convenience function."""
        bam_path = DATA_DIR / "healthy_sample.bam"
        if not bam_path.exists():
            pytest.skip("Test BAM not available")

        result = predict_cancer(bam_path)

        assert isinstance(result, PredictionResult)
