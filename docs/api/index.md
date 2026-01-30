# API Reference

This is the complete API reference for FragMentor.

## Main Module

::: fragmentomics

## Submodules

- [fragmentomics.features](features.md) - Feature extraction
- [fragmentomics.io](io.md) - I/O utilities  
- [fragmentomics.viz](viz.md) - Visualization

## Quick Example

```python
from fragmentomics import FragMentor

# High-level interface
fm = FragMentor("sample.bam")
dist = fm.sizes()
fig, ax = fm.plot_sizes()

# Or use lower-level functions
from fragmentomics.io import read_fragments
from fragmentomics.features import analyze_sizes

sizes = read_fragments("sample.bam")
dist = analyze_sizes(sizes)
```
