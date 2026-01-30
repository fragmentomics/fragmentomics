# FragMentor Open Source Wedge Strategy

## Competitive Landscape Analysis

### Direct Competitors

| Tool | Stars | Last Active | Focus | Weakness |
|------|-------|-------------|-------|----------|
| **cfDNApipe** | 75⭐ | Dec 2025 | Full WGS/WGBS pipeline | Complex, heavyweight, WGBS-focused |
| **FinaleToolkit** | 35⭐ | Jan 2026 | Fragment features | Tied to FinaleDB, academic speed |
| **LIONHEART** | 11⭐ | Jan 2026 | Cancer detection | Requires custom mosdepth, complex setup |
| **ctDNAtools** | ~20⭐ | 2022 | R package | Stale, R ecosystem |
| **celfie** | ~15⭐ | 2023 | Tissue deconvolution | Narrow focus |

### Gap Analysis

**What's missing in the market:**

1. **Simple, focused tool** — cfDNApipe is overwhelming (full pipeline), others are too narrow
2. **Modern Python** — Most tools use older patterns, no type hints, poor DX
3. **Speed** — No tool emphasizes performance (we use Polars, streaming)
4. **ML-ready output** — Feature matrices ready for scikit-learn/PyTorch
5. **Docker-first** — Easy reproducibility without conda hell
6. **Great documentation** — Most have sparse/academic docs

---

## The Wedge Strategy

### Our Positioning

> **"FastQC for cfDNA"** — The first tool everyone runs, the one everyone cites.

Not a full pipeline. Not a complex framework. Just the **essential fragmentomics analysis** that every researcher needs, done **fast and right**.

### Target User Journey

```
Researcher gets cfDNA BAM file
    ↓
Runs: pip install fragmentomics
    ↓  
Runs: fragmentomics sizes sample.bam -o results/
    ↓
Gets: publication-ready plots + feature matrix
    ↓
Uses features in their own ML pipeline
    ↓
Cites FragMentor in paper
```

### Wedge #1: Speed (10x faster)

**Why it matters:** WGS files are 50-200GB. Researchers wait hours.

**Our advantage:**
- Streaming BAM processing (constant memory)
- Polars instead of pandas (10x faster)
- Parallel batch processing
- Skip intermediate files

**Marketing angle:** "Analyze 100 samples before lunch"

### Wedge #2: ML-Ready Output

**Why it matters:** Every lab wants to build classifiers. Current tools output plots, not features.

**Our advantage:**
- Direct NumPy/Polars output
- Scikit-learn compatible feature matrices
- JSON/CSV for any language
- Documented feature definitions

**Marketing angle:** "From BAM to classifier in one notebook"

### Wedge #3: Zero-Config Docker

**Why it matters:** Conda environments break. Reproducibility matters for papers.

**Our advantage:**
- Single Docker image, everything works
- No reference genome management
- Consistent results across systems

**Marketing angle:** "Same results, every time, everywhere"

### Wedge #4: Documentation & Tutorials

**Why it matters:** Academic tools have terrible docs. Researchers waste hours.

**Our advantage:**
- Beautiful MkDocs site
- Jupyter notebooks for every use case
- Copy-paste CLI examples
- Feature definitions with citations

**Marketing angle:** "Learn fragmentomics in an afternoon"

---

## Growth Playbook

### Phase 1: Establish Credibility (Months 1-3)

1. **Reproduce published results**
   - Cristiano et al. 2019 (Science) — DELFI
   - Jiang et al. 2021 — Nucleosome footprinting
   - Mathios et al. 2021 — LUCAS study

2. **Benchmark against competitors**
   - Speed comparison (we should win by 10x)
   - Feature accuracy (validate against published)
   - Create reproducible benchmark suite

3. **Target key labs**
   - Reach out to 10 academic groups doing cfDNA research
   - Offer free support for adoption
   - Get testimonials and use cases

### Phase 2: Community Building (Months 3-6)

1. **Content marketing**
   - Tutorial blog posts
   - Twitter/Bluesky fragmentomics content
   - YouTube walkthrough videos

2. **Academic presence**
   - Present at AACR, ASHG, AGBT
   - Submit to JOSS or Bioinformatics
   - Get cited in papers

3. **Discord/GitHub community**
   - Answer every issue within 24h
   - Feature request voting
   - Contributor recognition

### Phase 3: Network Effects (Months 6-12)

1. **Integrations**
   - Nextflow/Snakemake workflows
   - Galaxy tool wrapper
   - Terra/DNAnexus compatibility

2. **Data flywheel**
   - Public benchmark dataset
   - Community-contributed features
   - Leaderboard for methods

3. **Commercial adjacent**
   - Consulting for pharma/biotech
   - Training workshops
   - Priority support tier

---

## Specific Tactical Actions

### This Week
- [ ] Reproduce one key paper's results (Cristiano 2019)
- [ ] Create benchmark suite vs FinaleToolkit
- [ ] Set up Twitter/Bluesky presence
- [ ] Identify 10 target labs to reach out to

### This Month
- [ ] Submit to Bioconda
- [ ] Write comparison blog post
- [ ] Create benchmark results page
- [ ] First academic lab adoption
- [ ] Submit to JOSS

### This Quarter
- [ ] 100+ GitHub stars
- [ ] 3+ external users
- [ ] Conference presentation submitted
- [ ] Integration with one workflow system

---

## Positioning Messages

### Taglines
- "See what others miss." (current)
- "From BAM to insight in minutes."
- "The fragmentomics toolkit that just works."
- "cfDNA analysis, done right."

### Elevator Pitch
> FragMentor is the open-source toolkit for cfDNA fragmentomics analysis. 
> While other tools try to do everything, we focus on extracting the features 
> that matter — fragment sizes, end motifs, nucleosome positioning — and 
> delivering them fast, in ML-ready formats. Install with pip, get results 
> in minutes, publish with confidence.

### Differentiation Statement
> Unlike cfDNApipe (a heavyweight pipeline) or FinaleToolkit (tied to a specific database),
> FragMentor is lightweight, fast, and designed for researchers who want to build
> their own analyses. It's the foundation others build on, not the whole house.

---

## Key Metrics to Track

| Metric | 3 Month Target | 6 Month Target |
|--------|----------------|----------------|
| GitHub stars | 100 | 500 |
| PyPI downloads/month | 500 | 2,000 |
| Academic users | 3 | 10 |
| Papers citing | 0 | 3 |
| Contributors | 2 | 5 |

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| FinaleToolkit gets more traction | Move faster, better docs, more features |
| Big company releases competing tool | Focus on open source community, stay nimble |
| Can't reproduce published results | Get help from original authors, validate carefully |
| No academic adoption | Partner with labs, offer co-authorship |

---

*Last updated: 2026-01-30*
