# End Motif Analysis

Fragment end motifs reveal enzymatic cleavage patterns.

## Overview

cfDNA fragment ends have non-random sequence compositions reflecting:

- DNase cleavage preferences
- Nuclease specificity
- Tissue-specific patterns

## Basic Usage

```python
from fragmentomics.features import analyze_end_motifs

profile = analyze_end_motifs("sample.bam", "reference.fa")
print(profile.summary())

# Get feature vector for ML
vec = profile.to_vector()  # 256-dim for 4-mers
```

## Key Features

| Feature | Description |
|---------|-------------|
| `entropy` | Shannon entropy of motif distribution |
| `diversity` | Effective number of motifs |
| `gc_content` | Average GC at fragment ends |
| `top_motifs` | Most frequent motifs |

## CLI Usage

```bash
fragmentomics motifs sample.bam -r hg38.fa -o results/
```
