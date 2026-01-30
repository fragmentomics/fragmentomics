# Nucleosome Positioning

Infer nucleosome positions from cfDNA fragmentation patterns.

## Overview

The Windowed Protection Score (WPS) measures the difference between:

- **Protection**: Fragments spanning a position (nucleosome-bound)
- **Cutting**: Fragments ending at a position (accessible DNA)

High WPS indicates nucleosome protection; low WPS indicates linker regions.

## Basic Usage

```python
from fragmentomics.features import compute_wps

profile = compute_wps("sample.bam", "chr1", 1000000, 1001000)
print(profile.summary())

# Access nucleosome peaks
print(f"Found {len(profile.peak_positions)} nucleosomes")
print(f"Periodicity: {profile.periodicity:.0f} bp")
```

## Key Outputs

| Output | Description |
|--------|-------------|
| `wps` | Raw WPS values |
| `smoothed_wps` | Smoothed WPS for visualization |
| `peak_positions` | Nucleosome center positions |
| `trough_positions` | Linker positions |
| `periodicity` | Estimated nucleosome repeat length |

## Expected Values

- **Nucleosome periodicity**: ~167 bp (repeat length)
- **Linker size**: ~20 bp between nucleosomes
