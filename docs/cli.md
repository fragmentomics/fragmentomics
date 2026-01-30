# CLI Reference

FragMentor provides a powerful command-line interface.

## Global Options

```bash
fragmentomics --version  # Show version
fragmentomics --help     # Show help
```

## Commands

### `sizes` - Fragment Size Analysis

Analyze fragment size distribution from a BAM file.

```bash
fragmentomics sizes INPUT.bam [OPTIONS]
```

**Options:**
- `-o, --output PATH` - Output directory (default: ./output)
- `--min-mapq INT` - Minimum mapping quality (default: 30)
- `--min-size INT` - Minimum fragment size (default: 50)
- `--max-size INT` - Maximum fragment size (default: 500)
- `-r, --region TEXT` - Genomic region (e.g., chr1:1000-2000)
- `--plot/--no-plot` - Generate plots (default: --plot)
- `--json` - Output JSON format

**Example:**
```bash
fragmentomics sizes sample.bam -o results/ --min-mapq 20
```

### `motifs` - End Motif Analysis

Analyze k-mer motifs at fragment ends.

```bash
fragmentomics motifs INPUT.bam -r REFERENCE.fa [OPTIONS]
```

**Options:**
- `-r, --reference PATH` - Reference genome FASTA (required)
- `-o, --output PATH` - Output directory
- `-k, --kmer-size INT` - K-mer size (default: 4)
- `--region TEXT` - Genomic region
- `--max-fragments INT` - Max fragments to analyze
- `--json` - Output JSON format

**Example:**
```bash
fragmentomics motifs sample.bam -r hg38.fa -o results/
```

### `extract` - Batch Feature Extraction

Extract multiple feature types at once.

```bash
fragmentomics extract INPUT.bam [OPTIONS]
```

**Options:**
- `-r, --reference PATH` - Reference genome (required for motifs)
- `-o, --output PATH` - Output directory
- `-f, --features TEXT` - Features to extract: sizes, motifs, all

**Example:**
```bash
fragmentomics extract sample.bam -r hg38.fa -f all -o features/
```

### `info` - BAM File Information

Display information about a BAM file.

```bash
fragmentomics info INPUT.bam
```

**Example:**
```bash
fragmentomics info sample.bam
```

## Output Files

| Command | Output Files |
|---------|-------------|
| `sizes` | `size_stats.json`, `size_distribution.png` |
| `motifs` | `motif_stats.json`, `motif_vector.npy` |
| `extract` | `features.json` |
