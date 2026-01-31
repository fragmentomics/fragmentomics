# FragMentor Benchmarks

Comparative benchmarks against other cfDNA fragmentomics tools.

## Competitors Tested

| Tool | Version | Install |
|------|---------|---------|
| FragMentor | 0.1.0-dev | `pip install fragmentomics` |
| FinaleToolkit | latest | `pip install finaletoolkit` |
| cfDNApipe | latest | conda only |

## Benchmark Categories

### 1. Speed Benchmarks
- Fragment size extraction (BAM → sizes)
- End motif analysis
- DELFI-style genome-wide profiling
- Batch processing (10, 50, 100 samples)

### 2. Memory Benchmarks
- Peak memory usage
- Memory scaling with file size

### 3. Accuracy Benchmarks
- Reproduce Cristiano et al. 2019 results
- Compare feature values on same input

## Test Data

### Synthetic Data
Generated with `tests/data/generate_test_bam.py`:
- `healthy_sample.bam` - 500 fragments, healthy-like distribution
- `cancer_sample.bam` - 500 fragments, cancer-like distribution

### Public Data (for full benchmarks)
- Snyder et al. 2016 (Cell) - downsampled BAMs available
- EGA datasets (pending access)

## Running Benchmarks

```bash
# Install dependencies
pip install fragmentomics finaletoolkit

# Run speed benchmarks
python benchmarks/speed_benchmark.py

# Run accuracy benchmarks
python benchmarks/accuracy_benchmark.py

# Generate report
python benchmarks/generate_report.py
```

## Results

See `benchmarks/results/` for latest benchmark results.
