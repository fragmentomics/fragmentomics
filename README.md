# 🧬 FragMentor

> **The definitive open-source toolkit for cfDNA fragmentomics analysis**

[![Tests](https://github.com/fragmentomics/fragmentomics/actions/workflows/test.yml/badge.svg)](https://github.com/fragmentomics/fragmentomics/actions/workflows/test.yml)
[![Documentation](https://github.com/fragmentomics/fragmentomics/actions/workflows/docs.yml/badge.svg)](https://fragmentomics.github.io/fragmentomics/)
[![Docker](https://github.com/fragmentomics/fragmentomics/actions/workflows/docker.yml/badge.svg)](https://github.com/fragmentomics/fragmentomics/pkgs/container/fragmentomics)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**From BAM to biological insight in minutes.**

---

## 🎯 What is FragMentor?

FragMentor analyzes cell-free DNA (cfDNA) fragmentation patterns from sequencing data. These patterns contain valuable information about tissue of origin and can detect cancer with high accuracy.

Cell-free DNA fragments in blood carry signatures of their origin through:
- **Fragment size** — Different tissues produce characteristic size distributions
- **End motifs** — Enzymatic cleavage patterns vary by tissue  
- **Nucleosome positioning** — Chromatin accessibility differs across cell types

FragMentor extracts these features and makes them accessible for research and clinical applications.

---

## ✨ Features

| Feature | Description |
|---------|-------------|
| 📏 **Fragment Sizes** | Distribution analysis, peaks, ratios, periodicity |
| 🔬 **DELFI Profile** | Genome-wide fragmentation (Cristiano et al. 2019) |
| 🎯 **Cancer Detection** | ML-based cancer probability from cfDNA patterns |
| 🧩 **End Motifs** | 4-mer frequencies at fragment ends |
| 📊 **Coverage** | Copy number estimation with GC correction |
| 📋 **Clinical Reports** | Self-contained HTML reports with visualizations |
| 🛡️ **Nucleosomes** | Windowed Protection Score (WPS) |
| 🤖 **ML-Ready** | Built-in feature extraction for machine learning |
| ⚡ **Fast** | Optimized with Polars — 10x faster than pandas |
| 🐳 **Containerized** | Docker support for reproducibility |

---

## 🚀 Quick Start

### Installation

```bash
pip install fragmentomics
```

With conda:
```bash
conda install -c conda-forge fragmentomics
```

With Docker:
```bash
docker pull ghcr.io/fragmentomics/fragmentomics:latest
```

### Command Line

```bash
# Generate clinical report (QC, sizes, DELFI, cancer prediction)
fragmentomics report sample.bam -o report.html

# Full report with all sections (requires reference for motifs)
fragmentomics report sample.bam -r hg38.fa -o report.html --sections all

# Select specific analyses
fragmentomics report sample.bam --sections sizes,delfi,prediction

# Individual analyses
fragmentomics sizes sample.bam -o results/
fragmentomics delfi sample.bam -o results/
fragmentomics coverage sample.bam -o results/ --bin-size 100000
fragmentomics motifs sample.bam -r hg38.fa -o results/

# Batch processing (parallel)
fragmentomics batch *.bam -o results/ --threads 8

# Extract ALL features for ML
fragmentomics extract sample.bam -r hg38.fa -o features/
```

**35x faster than FinaleToolkit** — see [benchmarks](benchmarks/).

### Python API

```python
from fragmentomics import FragMentor, ReportGenerator

# Initialize analyzer
fm = FragMentor("sample.bam", reference="hg38.fa")

# Fragment size analysis
sizes = fm.sizes()
print(f"Median: {sizes.median:.0f} bp")
print(f"Short fragment ratio: {sizes.ratio_short:.1%}")

# DELFI fragmentation profile
delfi = fm.delfi()
print(f"Genome-wide S/L ratio: {delfi.genome_wide_ratio:.3f}")

# Cancer prediction
prediction = fm.predict()
print(f"Cancer probability: {prediction.probability:.1%}")
print(f"Prediction: {prediction.prediction} ({prediction.confidence})")

# Generate clinical report
rg = ReportGenerator("sample.bam", reference="hg38.fa")
rg.add_sizes().add_delfi().add_prediction().add_coverage()
rg.generate("report.html")
```

---

## 📖 Documentation

Full documentation: [fragmentomics.github.io/fragmentomics](https://fragmentomics.github.io/fragmentomics)

- [Installation Guide](docs/installation.md)
- [Quick Start Tutorial](docs/quickstart.md)  
- [API Reference](docs/api/)
- [Example Notebooks](notebooks/)

---

## 🔬 Scientific Background

### Why Fragmentomics?

When cells die, they release DNA into the bloodstream as cell-free DNA (cfDNA). The fragmentation patterns of this DNA are non-random — they reflect:

1. **Nucleosome positioning** — DNA wraps around histones in ~147bp units
2. **Enzymatic cleavage** — Different nucleases cut at different motifs
3. **Chromatin state** — Open vs. closed chromatin affects fragmentation

Cancer cells have altered chromatin landscapes, which produces distinctive fragmentation signatures detectable in blood.

### Key References

- Cristiano S, et al. (2019). "Genome-wide cell-free DNA fragmentation in patients with cancer." *Nature*. [doi:10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6)

- Mouliere F, et al. (2018). "Enhanced detection of circulating tumor DNA by fragment size analysis." *Science Translational Medicine*. [doi:10.1126/scitranslmed.aat4921](https://doi.org/10.1126/scitranslmed.aat4921)

- Snyder MW, et al. (2016). "Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin." *Cell*. [doi:10.1016/j.cell.2015.11.050](https://doi.org/10.1016/j.cell.2015.11.050)

---

## 🤝 Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

```bash
# Development setup
git clone https://github.com/fragmentomics/fragmentomics
cd fragmentomics
pip install -e ".[dev]"
pre-commit install

# Run tests
pytest

# Build docs locally
mkdocs serve
```

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgments

FragMentor is built on the shoulders of giants:
- [pysam](https://github.com/pysam-developers/pysam) for BAM/CRAM handling
- [polars](https://github.com/pola-rs/polars) for blazing-fast data processing
- The cfDNA research community for foundational science

---

## ⚠️ Disclaimer

**Research Use Only** — FragMentor is intended for research purposes. It is not approved for clinical diagnosis or treatment decisions.

---

<p align="center">
  <b>FragMentor</b> — See what others miss. 🧬
</p>
