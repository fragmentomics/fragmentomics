# FragMentor Documentation

**🧬 The definitive open-source toolkit for cfDNA fragmentomics analysis**

---

## What is FragMentor?

FragMentor analyzes cell-free DNA (cfDNA) fragmentation patterns from sequencing data. These patterns contain valuable information about tissue of origin and can detect cancer with high accuracy.

**From BAM to biological insight in minutes.**

## Quick Start

```bash
# Install
pip install fragmentomics

# Analyze fragment sizes
fragmentomics sizes sample.bam -o results/

# Python API
from fragmentomics import FragMentor
fm = FragMentor("sample.bam")
dist = fm.sizes()
print(f"Short fragment ratio: {dist.ratio_short:.1%}")
```

## Features

- **Fragment Size Analysis** — Distribution, peaks, ratios, periodicity
- **End Motif Profiling** — K-mer frequencies at fragment ends
- **Publication-Ready Plots** — Beautiful visualizations out of the box
- **ML-Ready** — Built-in feature extraction for machine learning
- **Fast** — Optimized with streaming I/O and efficient algorithms

## Documentation

- [Installation](installation.md)
- [Quick Start Tutorial](quickstart.md)
- [API Reference](api/index.md)
- [CLI Reference](cli.md)

## Links

- [GitHub Repository](https://github.com/fragmentomics/fragmentomics)
- [Issue Tracker](https://github.com/fragmentomics/fragmentomics/issues)
- [Contributing Guide](https://github.com/fragmentomics/fragmentomics/blob/main/CONTRIBUTING.md)

---

**FragMentor** — See what others miss. 🧬
