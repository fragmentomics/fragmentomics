# 🧬 FragMentor

> **The definitive open-source toolkit for cfDNA fragmentomics analysis**

[![PyPI version](https://badge.fury.io/py/fragmentomics.svg)](https://badge.fury.io/py/fragmentomics)
[![Tests](https://github.com/fragmentomics/fragmentomics/workflows/Tests/badge.svg)](https://github.com/fragmentomics/fragmentomics/actions)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://fragmentomics.github.io/fragmentomics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
| 🧩 **End Motifs** | 4-mer frequencies at fragment ends |
| 📊 **Coverage** | Copy number estimation with GC correction |
| 🛡️ **Nucleosomes** | Windowed Protection Score (WPS) |
| 🎯 **Regions** | Analyze custom genomic regions (BED) |
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
# Analyze fragment sizes
fragmentomics sizes sample.bam -o results/

# Extract end motifs
fragmentomics motifs sample.bam -r hg38.fa -o results/

# Extract ALL features
fragmentomics extract sample.bam -r hg38.fa -o features/

# Run cancer detection model
fragmentomics predict sample.bam --model cancer_v1
```

### Python API

```python
from fragmentomics import FragMentor

# Initialize analyzer
fm = FragMentor("sample.bam", reference="hg38.fa")

# Fragment size analysis
sizes = fm.sizes()
print(f"Median: {sizes.median:.0f} bp")
print(f"Mode: {sizes.mode:.0f} bp")

# End motif analysis
motifs = fm.end_motifs(k=4)
print(f"Top motif: {motifs.most_common(1)}")

# Extract all features for ML
features = fm.extract_features()
features.to_parquet("sample_features.parquet")

# Cancer prediction (with pre-trained model)
prediction = fm.predict("cancer_v1")
print(f"Cancer probability: {prediction.score:.1%}")
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
