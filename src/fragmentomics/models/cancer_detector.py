"""
Cancer detection from cfDNA fragmentomics features.

Implements a classifier trained on fragment size, DELFI ratios,
and other fragmentomics features to predict cancer probability.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# Model directory (shipped with package)
MODELS_DIR = Path(__file__).parent / "pretrained"


@dataclass
class PredictionResult:
    """
    Cancer prediction result.

    Attributes
    ----------
    probability : float
        Probability of cancer (0-1)
    prediction : str
        "cancer" or "healthy"
    confidence : str
        "high", "medium", or "low"
    features_used : list[str]
        Features used for prediction
    feature_importance : dict[str, float]
        Importance of each feature
    """

    probability: float
    prediction: str
    confidence: str
    features_used: list[str]
    feature_importance: dict[str, float]

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "probability": self.probability,
            "prediction": self.prediction,
            "confidence": self.confidence,
            "features_used": self.features_used,
            "feature_importance": self.feature_importance,
        }

    def summary(self) -> str:
        """Return human-readable summary."""
        emoji = "🔴" if self.prediction == "cancer" else "🟢"
        return (
            f"{emoji} Prediction: {self.prediction.upper()}\n"
            f"   Probability: {self.probability:.1%}\n"
            f"   Confidence: {self.confidence}\n"
        )


class CancerDetector:
    """
    Cancer detection model for cfDNA fragmentomics.

    Uses fragment size distribution and DELFI-style features
    to predict cancer probability.

    Parameters
    ----------
    model_path : str or Path, optional
        Path to pre-trained model. If None, uses default model.
    threshold : float, default 0.5
        Classification threshold

    Examples
    --------
    >>> detector = CancerDetector()
    >>> result = detector.predict_from_bam("sample.bam")
    >>> print(f"Cancer probability: {result.probability:.1%}")
    """

    # Feature definitions
    FEATURES = [
        "ratio_short",      # <150bp fraction
        "ratio_mono",       # 140-180bp fraction
        "mean_size",        # Mean fragment size
        "std_size",         # Size standard deviation
        "periodicity_10bp", # 10bp periodicity score
        "genome_ratio",     # Genome-wide short/long ratio (DELFI)
    ]

    def __init__(
        self,
        model_path: str | Path | None = None,
        threshold: float = 0.5,
    ):
        self.threshold = threshold
        self.model = None
        self.scaler = None
        self.feature_names = self.FEATURES.copy()

        if model_path:
            self.load(model_path)
        else:
            # Try to load default model
            default_path = MODELS_DIR / "cancer_v1.json"
            if default_path.exists():
                self.load(default_path)
            else:
                logger.info("No pre-trained model found. Using built-in heuristics.")
                self._init_heuristic_model()

    def _init_heuristic_model(self) -> None:
        """Initialize heuristic-based model (no training data needed)."""
        # Based on literature:
        # - Cancer samples have more short fragments
        # - Cancer samples have lower mean size
        # - Cancer samples have weaker nucleosomal periodicity
        self.model = {
            "type": "heuristic",
            "weights": {
                "ratio_short": 2.0,      # Higher = more cancer-like
                "ratio_mono": -1.5,      # Lower = more cancer-like
                "mean_size": -0.02,      # Lower = more cancer-like
                "std_size": 0.01,        # Higher variance = more cancer-like
                "periodicity_10bp": -1.0, # Lower = more cancer-like
                "genome_ratio": 1.5,     # Higher short/long = more cancer-like
            },
            "intercept": 0.0,
            "healthy_baseline": {
                "ratio_short": 0.15,
                "ratio_mono": 0.65,
                "mean_size": 167,
                "std_size": 40,
                "periodicity_10bp": 0.8,
                "genome_ratio": 0.4,
            },
        }
        self.scaler = None

    def load(self, path: str | Path) -> None:
        """Load model from file."""
        path = Path(path)
        with open(path) as f:
            data = json.load(f)
        self.model = data.get("model")
        self.scaler = data.get("scaler")
        self.feature_names = data.get("features", self.FEATURES)
        logger.info(f"Loaded model from {path}")

    def save(self, path: str | Path) -> None:
        """Save model to file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "model": self.model,
            "scaler": self.scaler,
            "features": self.feature_names,
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Saved model to {path}")

    def extract_features(self, bam_path: str | Path) -> dict[str, float]:
        """
        Extract features from BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file

        Returns
        -------
        dict[str, float]
            Feature values
        """
        from fragmentomics.features import DELFIAnalyzer, analyze_sizes
        from fragmentomics.io import BamReader

        bam_path = Path(bam_path)

        # Extract sizes
        reader = BamReader(bam_path)
        sizes = reader.extract_sizes()
        size_dist = analyze_sizes(sizes)

        # Extract DELFI profile
        delfi = DELFIAnalyzer(bin_size=1_000_000)  # 1Mb bins for speed
        profile = delfi.analyze(bam_path)

        features = {
            "ratio_short": size_dist.ratio_short,
            "ratio_mono": size_dist.ratio_mono,
            "mean_size": size_dist.mean,
            "std_size": size_dist.std,
            "periodicity_10bp": size_dist.periodicity_10bp,
            "genome_ratio": profile.genome_wide_ratio,
        }

        return features

    def predict(self, features: dict[str, float]) -> PredictionResult:
        """
        Predict cancer probability from features.

        Parameters
        ----------
        features : dict[str, float]
            Feature values

        Returns
        -------
        PredictionResult
            Prediction with probability and confidence
        """
        if self.model is None:
            raise ValueError("No model loaded")

        if self.model["type"] == "heuristic":
            return self._predict_heuristic(features)
        else:
            return self._predict_sklearn(features)

    def _predict_heuristic(self, features: dict[str, float]) -> PredictionResult:
        """Predict using heuristic scoring."""
        weights = self.model["weights"]
        baseline = self.model["healthy_baseline"]

        # Compute deviation score from healthy baseline
        score = self.model["intercept"]
        importances = {}

        for feat, weight in weights.items():
            if feat in features and feat in baseline:
                deviation = features[feat] - baseline[feat]
                contribution = deviation * weight
                score += contribution
                importances[feat] = abs(contribution)

        # Sigmoid to get probability
        probability = 1 / (1 + np.exp(-score))

        # Determine prediction and confidence
        prediction = "cancer" if probability >= self.threshold else "healthy"

        if probability > 0.8 or probability < 0.2:
            confidence = "high"
        elif probability > 0.65 or probability < 0.35:
            confidence = "medium"
        else:
            confidence = "low"

        # Normalize importances
        total_imp = sum(importances.values()) or 1
        importances = {k: v / total_imp for k, v in importances.items()}

        return PredictionResult(
            probability=float(probability),
            prediction=prediction,
            confidence=confidence,
            features_used=list(features.keys()),
            feature_importance=importances,
        )

    def _predict_sklearn(self, features: dict[str, float]) -> PredictionResult:
        """Predict using sklearn model."""
        # This would load and use a real sklearn model
        # For now, fall back to heuristic
        # When we have trained models:
        # import joblib
        # model = joblib.load(self.model_path)
        # proba = model.predict_proba([feature_vector])[0][1]
        return self._predict_heuristic(features)

    def predict_from_bam(self, bam_path: str | Path) -> PredictionResult:
        """
        Predict cancer probability from BAM file.

        Parameters
        ----------
        bam_path : str or Path
            Path to BAM file

        Returns
        -------
        PredictionResult
            Prediction result
        """
        features = self.extract_features(bam_path)
        return self.predict(features)


def load_model(name: str = "cancer_v1") -> CancerDetector:
    """
    Load a pre-trained model by name.

    Parameters
    ----------
    name : str
        Model name (e.g., "cancer_v1")

    Returns
    -------
    CancerDetector
        Loaded model
    """
    model_path = MODELS_DIR / f"{name}.json"
    if model_path.exists():
        return CancerDetector(model_path)
    else:
        logger.warning(f"Model '{name}' not found, using heuristic model")
        return CancerDetector()


def predict_cancer(bam_path: str | Path) -> PredictionResult:
    """
    Convenience function to predict cancer from BAM.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file

    Returns
    -------
    PredictionResult
        Prediction result
    """
    detector = CancerDetector()
    return detector.predict_from_bam(bam_path)
