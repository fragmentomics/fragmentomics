# GC Bias Correction

GC content affects sequencing efficiency. This module corrects for GC bias.

## Overview

Fragments with extreme GC content may be under- or over-represented due to:

- PCR amplification bias
- Library preparation effects
- Sequencing efficiency

## Basic Usage

```python
from fragmentomics.features import compute_gc_bias, GCCorrector

# Compute bias profile
profile = compute_gc_bias("sample.bam", "reference.fa")
print(profile.summary())

# Apply correction
corrector = GCCorrector()
gc_values = corrector.compute_fragment_gc("sample.bam", "reference.fa")
profile = corrector.compute_bias_from_fragments(gc_values)
```

## Key Metrics

| Metric | Description |
|--------|-------------|
| `r_squared` | GC vs coverage correlation |
| `gc_dropout` | Fraction of bins with low coverage |
| `mean_gc` | Mean GC content |
| `correction_factors` | Per-bin correction multipliers |
