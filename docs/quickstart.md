# Quick Start

This guide will get you analyzing cfDNA in minutes.

## Basic Fragment Size Analysis

```python
from fragmentomics import FragMentor

# Load your BAM file
fm = FragMentor("sample.bam")

# Analyze fragment sizes
dist = fm.sizes()

# View summary
print(dist.summary())

# Key metrics
print(f"Median size: {dist.median:.0f} bp")
print(f"Short fragment ratio: {dist.ratio_short:.1%}")
print(f"Mononucleosome peak: {dist.peak_mono} bp")
```

## Visualization

```python
# Plot the distribution
fig, ax = fm.plot_sizes()
fig.savefig("fragment_sizes.png")
```

## Command Line Interface

```bash
# Analyze sizes
fragmentomics sizes sample.bam -o results/

# Analyze end motifs (requires reference)
fragmentomics motifs sample.bam -r hg38.fa -o results/

# Get BAM info
fragmentomics info sample.bam
```

## Feature Extraction for ML

```python
from fragmentomics.features import analyze_sizes

# Get features as dictionary
features = dist.to_dict()

# Use in pandas
import pandas as pd
df = pd.DataFrame([features])
```

## Next Steps

- [Fragment Size Analysis Guide](guide/sizes.md)
- [End Motif Analysis](guide/motifs.md)
- [API Reference](api/index.md)
