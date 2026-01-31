#!/usr/bin/env python3
"""
Speed benchmarks comparing FragMentor to other tools.

Usage:
    python benchmarks/speed_benchmark.py [--bam PATH] [--iterations N]
"""

import argparse
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np

# Results storage
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def benchmark_fragmentomics_sizes(bam_path: Path, iterations: int = 3) -> dict:
    """Benchmark FragMentor size extraction."""
    from fragmentomics.io import BamReader
    from fragmentomics.features import analyze_sizes

    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        reader = BamReader(bam_path)
        sizes = reader.extract_sizes()
        dist = analyze_sizes(sizes)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "tool": "fragmentomics",
        "operation": "size_extraction",
        "mean_time": np.mean(times),
        "std_time": np.std(times),
        "min_time": np.min(times),
        "iterations": iterations,
        "n_fragments": dist.n_fragments,
    }


def benchmark_finaletoolkit_sizes(bam_path: Path, iterations: int = 3) -> dict:
    """Benchmark FinaleToolkit size extraction."""
    try:
        import finaletoolkit as ft
    except ImportError:
        return {"tool": "finaletoolkit", "operation": "size_extraction", "error": "not installed"}

    times = []
    n_frags = 0
    for _ in range(iterations):
        start = time.perf_counter()
        try:
            # FinaleToolkit API: contig=None for genome-wide
            frags = list(ft.frag.frag_generator(str(bam_path), contig=None))
            sizes = [f[2] - f[1] for f in frags]  # stop - start = length
            n_frags = len(sizes)
        except Exception as e:
            return {"tool": "finaletoolkit", "operation": "size_extraction", "error": str(e)}
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "tool": "finaletoolkit",
        "operation": "size_extraction",
        "mean_time": np.mean(times),
        "std_time": np.std(times),
        "min_time": np.min(times),
        "iterations": iterations,
        "n_fragments": n_frags,
    }


def benchmark_fragmentomics_delfi(bam_path: Path, iterations: int = 3) -> dict:
    """Benchmark FragMentor DELFI-style analysis."""
    from fragmentomics.features import DELFIAnalyzer

    times = []
    n_bins = 0
    for _ in range(iterations):
        start = time.perf_counter()
        analyzer = DELFIAnalyzer(bin_size=100_000)  # Smaller bins for test data
        profile = analyzer.analyze(bam_path)
        n_bins = len(profile.bins)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "tool": "fragmentomics",
        "operation": "delfi_profile",
        "mean_time": np.mean(times),
        "std_time": np.std(times),
        "min_time": np.min(times),
        "iterations": iterations,
        "n_bins": n_bins,
    }


def benchmark_fragmentomics_wps(bam_path: Path, iterations: int = 3) -> dict:
    """Benchmark FragMentor WPS analysis."""
    from fragmentomics.features import NucleosomeAnalyzer

    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        analyzer = NucleosomeAnalyzer()
        # Use a small region for benchmarking
        profile = analyzer.compute_wps(bam_path, "chr1", 10000, 50000)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "tool": "fragmentomics",
        "operation": "wps_analysis",
        "mean_time": np.mean(times),
        "std_time": np.std(times),
        "min_time": np.min(times),
        "iterations": iterations,
        "region_size": 40000,
    }


def run_all_benchmarks(bam_path: Path, iterations: int = 3) -> list[dict]:
    """Run all benchmarks."""
    results = []

    print(f"Running benchmarks on: {bam_path}")
    print(f"Iterations per benchmark: {iterations}")
    print("-" * 50)

    # FragMentor benchmarks
    print("→ FragMentor size extraction...")
    results.append(benchmark_fragmentomics_sizes(bam_path, iterations))
    print(f"  {results[-1].get('mean_time', 'N/A'):.3f}s")

    print("→ FragMentor DELFI profile...")
    results.append(benchmark_fragmentomics_delfi(bam_path, iterations))
    print(f"  {results[-1].get('mean_time', 'N/A'):.3f}s")

    print("→ FragMentor WPS analysis...")
    results.append(benchmark_fragmentomics_wps(bam_path, iterations))
    print(f"  {results[-1].get('mean_time', 'N/A'):.3f}s")

    # FinaleToolkit benchmarks
    print("→ FinaleToolkit size extraction...")
    results.append(benchmark_finaletoolkit_sizes(bam_path, iterations))
    if "error" in results[-1]:
        print(f"  Error: {results[-1]['error']}")
    else:
        print(f"  {results[-1].get('mean_time', 'N/A'):.3f}s")

    return results


def generate_report(results: list[dict], output_path: Path) -> None:
    """Generate benchmark report."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "results": results,
    }

    # Save JSON
    json_path = output_path / "benchmark_results.json"
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nResults saved to: {json_path}")

    # Print summary table
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    print(f"{'Tool':<20} {'Operation':<25} {'Mean Time':>10}")
    print("-" * 60)
    for r in results:
        if "error" in r:
            print(f"{r['tool']:<20} {r['operation']:<25} {'ERROR':>10}")
        else:
            print(f"{r['tool']:<20} {r['operation']:<25} {r['mean_time']:>10.3f}s")


def main():
    parser = argparse.ArgumentParser(description="Run FragMentor speed benchmarks")
    parser.add_argument(
        "--bam",
        type=Path,
        default=Path(__file__).parent.parent / "tests/data/healthy_sample.bam",
        help="Path to BAM file for benchmarking",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=3,
        help="Number of iterations per benchmark",
    )
    args = parser.parse_args()

    if not args.bam.exists():
        print(f"Error: BAM file not found: {args.bam}")
        print("Generate test data first: python tests/data/generate_test_bam.py")
        return 1

    results = run_all_benchmarks(args.bam, args.iterations)
    generate_report(results, RESULTS_DIR)

    return 0


if __name__ == "__main__":
    exit(main())
