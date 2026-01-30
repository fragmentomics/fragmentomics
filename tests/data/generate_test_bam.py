#!/usr/bin/env python3
"""Generate synthetic BAM files for testing fragmentomics."""

import pysam
import random
from pathlib import Path

# Output directory
DATA_DIR = Path(__file__).parent


def generate_test_bam(
    output_path: Path,
    n_fragments: int = 1000,
    size_distribution: str = "healthy",
    seed: int = 42,
) -> dict:
    """Generate a synthetic BAM file with known characteristics.
    
    Args:
        output_path: Path to write BAM file
        n_fragments: Number of fragments to generate
        size_distribution: "healthy" (peak ~167bp) or "cancer" (more short fragments)
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with ground truth statistics
    """
    random.seed(seed)
    
    # Define size distributions
    if size_distribution == "healthy":
        # Healthy: strong nucleosomal peak around 167bp
        sizes = []
        for _ in range(n_fragments):
            if random.random() < 0.8:
                # Mononucleosomal peak
                sizes.append(int(random.gauss(167, 15)))
            elif random.random() < 0.5:
                # Dinucleosomal
                sizes.append(int(random.gauss(334, 20)))
            else:
                # Short fragments
                sizes.append(int(random.gauss(100, 30)))
    else:  # cancer
        # Cancer: more short fragments, weaker nucleosomal structure
        sizes = []
        for _ in range(n_fragments):
            if random.random() < 0.5:
                # Mononucleosomal (weaker)
                sizes.append(int(random.gauss(167, 25)))
            elif random.random() < 0.3:
                # Short fragments (more prevalent)
                sizes.append(int(random.gauss(100, 30)))
            else:
                # Very short (tumor-derived)
                sizes.append(int(random.gauss(80, 20)))
    
    # Clamp sizes to reasonable range
    sizes = [max(30, min(600, s)) for s in sizes]
    
    # Create BAM header
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 248956422},
            {"SN": "chr2", "LN": 242193529},
        ],
        "RG": [{"ID": "test", "SM": "test_sample"}],
    }
    
    # Generate synthetic sequence
    bases = "ACGT"
    
    # Collect all reads first, then sort by position
    all_reads = []
    pos = 10000  # Starting position
    
    for i, frag_size in enumerate(sizes):
        # Use shorter read length for small fragments
        read_len = min(100, frag_size - 10)  # Ensure fragment > 2*read overlap
        read_len = max(50, read_len)  # Minimum read length
        
        # Calculate positions ensuring read2 starts after read1 for proper sorting
        read1_start = pos
        read2_start = pos + frag_size - read_len
        
        # Read 1 (forward)
        read1 = pysam.AlignedSegment()
        read1.query_name = f"read_{i:06d}"
        read1.query_sequence = "".join(random.choices(bases, k=read_len))
        read1.flag = 99  # paired, proper pair, mate reverse, first in pair
        read1.reference_id = 0  # chr1
        read1.reference_start = read1_start
        read1.mapping_quality = 60
        read1.cigar = [(0, read_len)]
        read1.next_reference_id = 0
        read1.next_reference_start = read2_start
        read1.template_length = frag_size
        read1.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        read1.set_tag("RG", "test")
        
        # Read 2 (reverse)
        read2 = pysam.AlignedSegment()
        read2.query_name = f"read_{i:06d}"
        read2.query_sequence = "".join(random.choices(bases, k=read_len))
        read2.flag = 147  # paired, proper pair, reverse, second in pair
        read2.reference_id = 0
        read2.reference_start = read2_start
        read2.mapping_quality = 60
        read2.cigar = [(0, read_len)]
        read2.next_reference_id = 0
        read2.next_reference_start = read1_start
        read2.template_length = -frag_size
        read2.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        read2.set_tag("RG", "test")
        
        all_reads.append((read1_start, read1))
        all_reads.append((read2_start, read2))
        
        # Advance position (ensure no overlap between fragment pairs)
        pos += frag_size + 500
    
    # Sort reads by position and write
    all_reads.sort(key=lambda x: x[0])
    
    with pysam.AlignmentFile(str(output_path), "wb", header=header) as outf:
        for _, read in all_reads:
            outf.write(read)
    
    # Index the BAM
    pysam.index(str(output_path))
    
    # Calculate ground truth statistics
    short_count = sum(1 for s in sizes if s < 150)
    mono_count = sum(1 for s in sizes if 150 <= s < 250)
    
    stats = {
        "n_fragments": n_fragments,
        "size_distribution": size_distribution,
        "mean_size": sum(sizes) / len(sizes),
        "median_size": sorted(sizes)[len(sizes) // 2],
        "short_ratio": short_count / n_fragments,
        "mononucleosomal_ratio": mono_count / n_fragments,
        "sizes": sizes,
    }
    
    return stats


def main():
    """Generate test BAM files."""
    # Generate healthy sample
    healthy_stats = generate_test_bam(
        DATA_DIR / "healthy_sample.bam",
        n_fragments=500,
        size_distribution="healthy",
        seed=42,
    )
    print(f"Generated healthy_sample.bam:")
    print(f"  Mean size: {healthy_stats['mean_size']:.1f}")
    print(f"  Median size: {healthy_stats['median_size']}")
    print(f"  Short ratio: {healthy_stats['short_ratio']:.3f}")
    print(f"  Mononucleosomal ratio: {healthy_stats['mononucleosomal_ratio']:.3f}")
    
    # Generate cancer sample
    cancer_stats = generate_test_bam(
        DATA_DIR / "cancer_sample.bam",
        n_fragments=500,
        size_distribution="cancer",
        seed=123,
    )
    print(f"\nGenerated cancer_sample.bam:")
    print(f"  Mean size: {cancer_stats['mean_size']:.1f}")
    print(f"  Median size: {cancer_stats['median_size']}")
    print(f"  Short ratio: {cancer_stats['short_ratio']:.3f}")
    print(f"  Mononucleosomal ratio: {cancer_stats['mononucleosomal_ratio']:.3f}")


if __name__ == "__main__":
    main()
