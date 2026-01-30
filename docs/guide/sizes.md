# Fragment Size Analysis

Fragment size distribution is the most fundamental cfDNA analysis.

## Overview

Cell-free DNA fragments have characteristic sizes reflecting nucleosome structure:

- **Mononucleosome peak**: ~167 bp (DNA wrapped around one histone octamer)
- **Dinucleosome peak**: ~334 bp (two nucleosomes + linker)
- **Short fragments**: <150 bp (often tumor-enriched)

## Basic Usage

```python
from fragmentomics import FragMentor

fm = FragMentor("sample.bam")
dist = fm.sizes()

print(f"Median: {dist.median:.0f} bp")
print(f"Short ratio: {dist.ratio_short:.1%}")
```

## Key Features

| Feature | Description |
|---------|-------------|
| `ratio_short` | Fraction of fragments <150 bp |
| `ratio_mono` | Fraction in 140-180 bp range |
| `ratio_di` | Fraction in 280-360 bp range |
| `peak_mono` | Position of mononucleosome peak |
| `peak_di` | Position of dinucleosome peak |
| `periodicity_10bp` | 10-bp periodicity score |

## CLI Usage

```bash
fragmentomics sizes sample.bam -o results/
```
