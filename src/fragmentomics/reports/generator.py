"""
Report generator for cfDNA fragmentomics analysis.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Self

from jinja2 import Environment, PackageLoader

import fragmentomics
from fragmentomics import FragMentor


class ReportGenerator:
    """
    Generate HTML/PDF reports for cfDNA fragmentomics analysis.

    Parameters
    ----------
    bam_path : str or Path
        Path to BAM file
    reference : str or Path, optional
        Path to reference genome FASTA (required for motifs)
    sample_name : str, optional
        Sample name for report header (defaults to BAM filename)

    Examples
    --------
    >>> report = ReportGenerator("sample.bam", reference="hg38.fa")
    >>> report.add_sizes().add_motifs().generate("report.html")

    >>> # Or generate everything at once
    >>> report.add_all().generate("report.html")
    """

    def __init__(
        self,
        bam_path: str | Path,
        reference: str | Path | None = None,
        sample_name: str | None = None,
    ) -> None:
        self.bam_path = Path(bam_path)
        self.reference = Path(reference) if reference else None
        self.sample_name = sample_name or self.bam_path.stem

        # Initialize FragMentor for analysis
        self._fm = FragMentor(str(self.bam_path), reference=str(self.reference) if self.reference else None)

        # Sections to include
        self._sections: list[dict[str, Any]] = []
        self._include_qc = True

        # Cached analysis results
        self._sizes = None
        self._motifs = None
        self._coverage = None

        # Template environment
        self._env = Environment(
            loader=PackageLoader("fragmentomics.reports", "templates"),
            autoescape=True,
        )

    def add_sizes(self, **kwargs: Any) -> Self:
        """
        Add fragment size distribution section.

        Parameters
        ----------
        **kwargs
            Passed to FragMentor.sizes()

        Returns
        -------
        Self
            For method chaining
        """
        self._sizes = self._fm.sizes(**kwargs)
        self._sections.append({
            "type": "sizes",
            "title": "Fragment Size Distribution",
            "data": self._sizes,
        })
        return self

    def add_motifs(self, k: int = 4, **kwargs: Any) -> Self:
        """
        Add end motif analysis section.

        Parameters
        ----------
        k : int, default 4
            K-mer length for motif analysis
        **kwargs
            Passed to FragMentor.end_motifs()

        Returns
        -------
        Self
            For method chaining

        Raises
        ------
        ValueError
            If no reference genome was provided
        """
        if not self.reference:
            raise ValueError("Reference genome required for motif analysis")

        self._motifs = self._fm.end_motifs(k=k, **kwargs)
        self._sections.append({
            "type": "motifs",
            "title": "End Motif Analysis",
            "data": self._motifs,
        })
        return self

    def add_coverage(self, bin_size: int = 100_000, **kwargs: Any) -> Self:
        """
        Add coverage/CNV section.

        Parameters
        ----------
        bin_size : int, default 100000
            Bin size for coverage analysis
        **kwargs
            Passed to FragMentor.coverage()

        Returns
        -------
        Self
            For method chaining
        """
        self._coverage = self._fm.coverage(bin_size=bin_size, **kwargs)
        self._sections.append({
            "type": "coverage",
            "title": "Genome-Wide Coverage",
            "data": self._coverage,
        })
        return self

    def add_delfi(self, bin_size: int = 5_000_000, **kwargs: Any) -> Self:
        """
        Add DELFI fragmentation section.

        Parameters
        ----------
        bin_size : int, default 5_000_000
            Bin size for DELFI analysis (5 Mb recommended)
        **kwargs
            Passed to FragMentor.delfi()

        Returns
        -------
        Self
            For method chaining
        """
        self._delfi = self._fm.delfi(bin_size=bin_size, **kwargs)
        self._sections.append({
            "type": "delfi",
            "title": "DELFI Fragmentation Profile",
            "data": self._delfi,
        })
        return self

    def add_prediction(self, **kwargs: Any) -> Self:
        """
        Add ML cancer prediction section.

        Parameters
        ----------
        **kwargs
            Passed to FragMentor.predict()

        Returns
        -------
        Self
            For method chaining
        """
        self._prediction = self._fm.predict(**kwargs)
        self._sections.append({
            "type": "prediction",
            "title": "Cancer Risk Assessment",
            "data": self._prediction,
        })
        return self

    def add_all(self) -> Self:
        """
        Add all available analysis sections.

        Adds sizes, DELFI, prediction, motifs (if reference), and coverage.

        Returns
        -------
        Self
            For method chaining
        """
        self.add_sizes()
        self.add_delfi()
        self.add_prediction()
        if self.reference:
            self.add_motifs()
        self.add_coverage()
        return self

    def generate(
        self,
        output: str | Path,
        format: str = "html",
        title: str | None = None,
    ) -> Path:
        """
        Generate the report.

        Parameters
        ----------
        output : str or Path
            Output file path
        format : str, default "html"
            Output format ("html" or "pdf")
        title : str, optional
            Report title (defaults to "FragMentor Report")

        Returns
        -------
        Path
            Path to generated report
        """
        output = Path(output)

        # Build context for template
        context = self._build_context(title)

        # Render template
        template = self._env.get_template("report.html")
        html_content = template.render(**context)

        if format == "html":
            output.write_text(html_content)
        elif format == "pdf":
            self._generate_pdf(html_content, output)
        else:
            raise ValueError(f"Unsupported format: {format}")

        return output

    def _build_context(self, title: str | None) -> dict[str, Any]:
        """Build template context with all sections."""
        # Lazy import to avoid circular dependency
        from fragmentomics.reports.sections import (
            build_coverage_section,
            build_delfi_section,
            build_motifs_section,
            build_prediction_section,
            build_qc_section,
            build_sizes_section,
        )

        # Build section content
        rendered_sections = []

        # Always include QC section first
        if self._include_qc:
            qc_section = build_qc_section(self._fm, self._sizes)
            rendered_sections.append(qc_section)

        # Add user-selected sections
        for section in self._sections:
            if section["type"] == "sizes":
                rendered_sections.append(
                    build_sizes_section(section["data"])
                )
            elif section["type"] == "delfi":
                rendered_sections.append(
                    build_delfi_section(section["data"])
                )
            elif section["type"] == "prediction":
                rendered_sections.append(
                    build_prediction_section(section["data"])
                )
            elif section["type"] == "motifs":
                rendered_sections.append(
                    build_motifs_section(section["data"])
                )
            elif section["type"] == "coverage":
                rendered_sections.append(
                    build_coverage_section(section["data"])
                )

        return {
            "title": title or "FragMentor Report",
            "sample_name": self.sample_name,
            "bam_path": str(self.bam_path),
            "reference": str(self.reference) if self.reference else "N/A",
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "version": fragmentomics.__version__,
            "sections": rendered_sections,
        }

    def _generate_pdf(self, html_content: str, output: Path) -> None:
        """Generate PDF from HTML content."""
        try:
            from weasyprint import HTML
        except ImportError:
            raise ImportError(
                "PDF generation requires weasyprint. "
                "Install with: pip install weasyprint"
            )

        HTML(string=html_content).write_pdf(output)
