# Changelog

All notable changes to FragMentor will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of FragMentor cfDNA fragmentomics toolkit
- Fragment size distribution analysis with nucleosomal peak detection
- Fragment end motif (4-mer) analysis
- GC bias correction module
- Nucleosome positioning / Windowed Protection Score (WPS)
- Command-line interface with `sizes`, `motifs`, `extract`, and `info` commands
- Visualization module with publication-ready plots
- Comprehensive documentation with MkDocs
- Docker support with multi-stage builds
- GitHub Actions CI/CD for testing, documentation, and releases
- 56 unit and integration tests with 48% code coverage

### Technical
- Streaming BAM processing for memory efficiency
- Polars-based data processing for performance
- Type hints throughout the codebase
- MIT licensed

## [0.1.0] - TBD

Initial public release.

---

[Unreleased]: https://github.com/fragmentomics/fragmentomics/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/fragmentomics/fragmentomics/releases/tag/v0.1.0
