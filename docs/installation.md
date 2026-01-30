# Installation

## Quick Install

```bash
pip install fragmentomics
```

## From Source

```bash
git clone https://github.com/fragmentomics/fragmentomics.git
cd fragmentomics
pip install -e .
```

## With Development Dependencies

```bash
pip install -e ".[dev]"
```

## Requirements

- Python 3.10 or higher
- pysam (for BAM/CRAM handling)
- numpy, scipy, scikit-learn
- matplotlib, plotly (for visualization)

## Verify Installation

```python
import fragmentomics
print(fragmentomics.__version__)
```

Or via CLI:

```bash
fragmentomics --version
```
