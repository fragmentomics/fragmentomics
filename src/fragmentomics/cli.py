"""
FragMentor CLI - Command-line interface for cfDNA fragmentomics analysis.

See what others miss. 🧬
"""

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn
from rich.table import Table

app = typer.Typer(
    name="fragmentomics",
    help="🧬 FragMentor — The definitive toolkit for cfDNA fragmentomics analysis",
    add_completion=False,
)
console = Console()


def version_callback(value: bool):
    if value:
        from fragmentomics import __brand__, __version__

        console.print(f"[bold blue]🧬 {__brand__}[/bold blue] v{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit",
    ),
):
    """
    🧬 FragMentor — The definitive toolkit for cfDNA fragmentomics analysis.

    From BAM to biological insight in minutes. See what others miss.
    """
    pass


@app.command()
def sizes(
    bam_path: Path = typer.Argument(..., help="Input BAM/CRAM file"),
    output: Path = typer.Option("./output", "-o", "--output", help="Output directory"),
    min_mapq: int = typer.Option(30, "--min-mapq", help="Minimum mapping quality"),
    min_size: int = typer.Option(50, "--min-size", help="Minimum fragment size"),
    max_size: int = typer.Option(500, "--max-size", help="Maximum fragment size"),
    plot: bool = typer.Option(True, "--plot/--no-plot", help="Generate plots"),
    region: str | None = typer.Option(
        None, "-r", "--region", help="Genomic region (chr:start-end)"
    ),
    json_output: bool = typer.Option(False, "--json", help="Output JSON format"),
):
    """
    Analyze fragment size distribution from BAM file.

    Extracts fragment sizes from properly-paired reads and computes
    fragmentomics features including size ratios, peaks, and periodicity.
    """
    console.print("[bold blue]🧬 FragMentor[/bold blue] — Fragment Size Analysis")
    console.print()

    # Validate input
    if not bam_path.exists():
        console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
        raise typer.Exit(1)

    # Create output directory
    output.mkdir(parents=True, exist_ok=True)

    # Import here to avoid slow startup
    from fragmentomics.features.sizes import FragmentSizeAnalyzer
    from fragmentomics.io.bam import BamReader
    from fragmentomics.viz.plots import plot_size_distribution, save_figure

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        console=console,
    ) as progress:
        # Extract fragments
        task = progress.add_task("Reading BAM file...", total=None)

        reader = BamReader(
            bam_path,
            min_mapq=min_mapq,
            min_size=min_size,
            max_size=max_size,
        )
        raw_sizes = reader.extract_sizes(region=region)
        progress.update(task, description=f"Extracted {len(raw_sizes):,} fragments")

        # Analyze
        progress.update(task, description="Computing features...")
        analyzer = FragmentSizeAnalyzer(min_size=min_size, max_size=max_size)
        dist = analyzer.analyze(raw_sizes)

        progress.update(task, completed=True, description="Analysis complete!")

    console.print()

    # Display results
    if json_output:
        print(json.dumps(dist.to_dict(), indent=2))
    else:
        console.print(dist.summary())

        # Show read stats
        console.print()
        console.print("[dim]Read Statistics:[/dim]")
        console.print(f"[dim]{reader.stats}[/dim]")

    # Save outputs
    stats_file = output / "size_stats.json"
    with open(stats_file, "w") as f:
        json.dump(dist.to_dict(), f, indent=2)
    console.print(f"\n[green]✓[/green] Stats saved to: {stats_file}")

    if plot:
        plot_file = output / "size_distribution.png"
        fig, ax = plot_size_distribution(dist, title=bam_path.stem)
        save_figure(fig, plot_file)
        console.print(f"[green]✓[/green] Plot saved to: {plot_file}")


@app.command()
def motifs(
    bam_path: Path = typer.Argument(..., help="Input BAM/CRAM file"),
    reference: Path = typer.Option(
        ..., "-r", "--reference", help="Reference genome FASTA"
    ),
    output: Path = typer.Option("./output", "-o", "--output", help="Output directory"),
    kmer_size: int = typer.Option(
        4, "-k", "--kmer-size", help="K-mer size for motif analysis"
    ),
    region: str | None = typer.Option(None, "--region", help="Genomic region"),
    max_fragments: int | None = typer.Option(
        None, "--max-fragments", help="Max fragments to analyze"
    ),
    json_output: bool = typer.Option(False, "--json", help="Output JSON format"),
):
    """
    Analyze fragment end motifs.

    Extracts k-mer sequences at fragment ends to compute motif
    frequencies and derived features like entropy and GC content.
    """
    console.print("[bold blue]🧬 FragMentor[/bold blue] — End Motif Analysis")
    console.print()

    # Validate inputs
    if not bam_path.exists():
        console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
        raise typer.Exit(1)
    if not reference.exists():
        console.print(f"[red]Error:[/red] Reference not found: {reference}")
        raise typer.Exit(1)

    output.mkdir(parents=True, exist_ok=True)

    from fragmentomics.features.motifs import EndMotifAnalyzer

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Analyzing end motifs...", total=None)

        analyzer = EndMotifAnalyzer(k=kmer_size)
        profile = analyzer.analyze_bam(
            bam_path,
            reference,
            region=region,
            max_fragments=max_fragments,
        )

        progress.update(task, description="Analysis complete!")

    console.print()

    if json_output:
        print(json.dumps(profile.to_dict(), indent=2))
    else:
        console.print(profile.summary())

    # Save outputs
    stats_file = output / "motif_stats.json"
    with open(stats_file, "w") as f:
        json.dump(profile.to_dict(), f, indent=2)
    console.print(f"\n[green]✓[/green] Stats saved to: {stats_file}")

    # Save feature vector
    vector_file = output / "motif_vector.npy"
    import numpy as np

    np.save(vector_file, profile.to_vector())
    console.print(f"[green]✓[/green] Feature vector saved to: {vector_file}")


@app.command()
def extract(
    bam_path: Path = typer.Argument(..., help="Input BAM/CRAM file"),
    reference: Path = typer.Option(
        None, "-r", "--reference", help="Reference genome FASTA"
    ),
    output: Path = typer.Option("./output", "-o", "--output", help="Output directory"),
    features: str = typer.Option(
        "sizes",
        "-f",
        "--features",
        help="Features to extract: sizes, motifs, all (motifs requires reference)",
    ),
):
    """
    Extract fragmentomic features from BAM file.

    Runs specified analyses and saves results for downstream ML/analysis.
    """
    console.print("[bold blue]🧬 FragMentor[/bold blue] — Feature Extraction")
    console.print()

    if not bam_path.exists():
        console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
        raise typer.Exit(1)

    output.mkdir(parents=True, exist_ok=True)

    feature_list = [f.strip().lower() for f in features.split(",")]
    if "all" in feature_list:
        feature_list = ["sizes", "motifs"]

    if "motifs" in feature_list and not reference:
        console.print("[red]Error:[/red] Reference required for motif analysis")
        raise typer.Exit(1)

    results = {}

    if "sizes" in feature_list:
        console.print("[bold]→ Extracting fragment sizes...[/bold]")
        from fragmentomics.features.sizes import FragmentSizeAnalyzer
        from fragmentomics.io.bam import BamReader

        reader = BamReader(bam_path)
        raw_sizes = reader.extract_sizes()
        analyzer = FragmentSizeAnalyzer()
        dist = analyzer.analyze(raw_sizes)
        results["sizes"] = dist.to_dict()
        console.print(f"  [green]✓[/green] {dist.n_fragments:,} fragments analyzed")

    if "motifs" in feature_list:
        console.print("[bold]→ Extracting end motifs...[/bold]")
        from fragmentomics.features.motifs import EndMotifAnalyzer

        analyzer = EndMotifAnalyzer()
        profile = analyzer.analyze_bam(bam_path, reference)
        results["motifs"] = profile.to_dict()
        console.print(f"  [green]✓[/green] {profile.n_ends:,} motifs analyzed")

    # Save combined results
    output_file = output / "features.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    console.print(f"\n[green]✓[/green] Features saved to: {output_file}")


@app.command()
def info(
    bam_path: Path = typer.Argument(..., help="Input BAM/CRAM file"),
):
    """
    Show information about a BAM file.
    """
    console.print("[bold blue]🧬 FragMentor[/bold blue] — BAM Info")
    console.print()

    if not bam_path.exists():
        console.print(f"[red]Error:[/red] File not found: {bam_path}")
        raise typer.Exit(1)

    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        # Basic stats
        table = Table(title=f"BAM File: {bam_path.name}")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("Format", "CRAM" if bam.is_cram else "BAM")
        table.add_row("References", str(bam.nreferences))

        # Count some reads
        count = 0
        proper_pairs = 0
        for read in bam.fetch():
            count += 1
            if read.is_proper_pair:
                proper_pairs += 1
            if count >= 10000:
                break

        table.add_row("Reads sampled", f"{count:,}")
        table.add_row(
            "Proper pairs",
            f"{proper_pairs:,} ({100*proper_pairs/count:.1f}%)" if count else "N/A",
        )

        console.print(table)


@app.command()
def batch(
    bam_files: list[Path] = typer.Argument(..., help="Input BAM/CRAM files"),
    output: Path = typer.Option("./batch_output", "-o", "--output", help="Output directory"),
    threads: int = typer.Option(1, "-t", "--threads", help="Number of parallel threads"),
    features: str = typer.Option("sizes", "-f", "--features", help="Features to extract"),
):
    """
    Batch analyze multiple BAM files.

    Example:
        fragmentomics batch sample1.bam sample2.bam -o results/
        fragmentomics batch *.bam -o results/ --threads 8
    """
    import json
    from concurrent.futures import ProcessPoolExecutor, as_completed

    console.print("[bold blue]🧬 FragMentor[/bold blue] — Batch Analysis")
    console.print()

    # Validate inputs
    valid_files = []
    for f in bam_files:
        if f.exists():
            valid_files.append(f)
        else:
            console.print(f"[yellow]Warning:[/yellow] File not found: {f}")

    if not valid_files:
        console.print("[red]Error:[/red] No valid BAM files found")
        raise typer.Exit(1)

    console.print(f"Processing {len(valid_files)} files with {threads} thread(s)\n")

    output.mkdir(parents=True, exist_ok=True)

    def analyze_file(bam_path: Path) -> dict:
        """Analyze a single BAM file."""
        from fragmentomics.features.sizes import FragmentSizeAnalyzer
        from fragmentomics.io.bam import BamReader

        try:
            reader = BamReader(bam_path)
            sizes = reader.extract_sizes()

            if len(sizes) == 0:
                return {"sample": bam_path.stem, "error": "No fragments extracted"}

            analyzer = FragmentSizeAnalyzer()
            dist = analyzer.analyze(sizes)

            result = dist.to_dict()
            result["sample"] = bam_path.stem
            result["path"] = str(bam_path)
            return result

        except Exception as e:
            return {"sample": bam_path.stem, "error": str(e)}

    results = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        console=console,
    ) as progress:
        task = progress.add_task("Analyzing...", total=len(valid_files))

        if threads == 1:
            # Sequential processing
            for bam_path in valid_files:
                result = analyze_file(bam_path)
                results.append(result)
                progress.advance(task)
                if "error" not in result:
                    console.print(f"  [green]✓[/green] {bam_path.name}")
                else:
                    console.print(f"  [red]✗[/red] {bam_path.name}: {result['error']}")
        else:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(analyze_file, f): f for f in valid_files}
                for future in as_completed(futures):
                    bam_path = futures[future]
                    result = future.result()
                    results.append(result)
                    progress.advance(task)
                    if "error" not in result:
                        console.print(f"  [green]✓[/green] {bam_path.name}")
                    else:
                        console.print(f"  [red]✗[/red] {bam_path.name}: {result['error']}")

    # Save combined results
    output_file = output / "batch_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    # Also save as CSV
    csv_file = output / "batch_results.csv"
    successful = [r for r in results if "error" not in r]
    if successful:
        import csv
        fieldnames = list(successful[0].keys())
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(successful)
        console.print(f"\n[green]✓[/green] Results saved to: {csv_file}")

    console.print(f"[green]✓[/green] JSON saved to: {output_file}")

    # Summary
    n_success = len([r for r in results if "error" not in r])
    n_error = len(results) - n_success
    console.print(f"\n[bold]Summary:[/bold] {n_success} succeeded, {n_error} failed")


if __name__ == "__main__":
    app()
