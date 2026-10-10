# E. coli Olaquindox RNA-seq Reanalysis (PRJNA525186)

A fully reproducible reanalysis of a public bacterial RNA-seq dataset,
built end-to-end on a laptop: raw SRA reads → quality control →
transcript quantification → differential expression → publication-style figures.

> **Status:** complete — ran end-to-end 2026-10-09. Every number below is measured.

## Biological question

How does *Escherichia coli* ATCC 25922 (a standard antimicrobial-susceptibility
reference strain) remodel its transcriptome when exposed to sub-inhibitory
concentrations of the veterinary antibiotic olaquindox?

## Dataset

- NCBI BioProject PRJNA525186 / SRA study SRP187444
- 9 samples, biological triplicates: untreated control (B1–B3),
  1/10× MIC olaquindox (D1–D3), 1/2× MIC olaquindox (C1–C3)
- Run accessions SRR8664715–SRR8664723 (see data/samples.tsv; `data/SraRunInfo.csv` is produced in Stage 1)

## Pipeline

1. Download: SRA Toolkit (`prefetch` + `fasterq-dump --split-files`) — 18 paired-end FASTQ files
2. QC: FastQC, summarised with MultiQC
3. Quantification: Salmon (bias-corrected, validated mappings) against a
   transcriptome built with gffread from the strain-matched ATCC 25922
   assembly GCF_000743255.1 — **not** K-12 MG1655
4. Differential expression: PyDESeq2, design `~ condition`, contrasts
   1/10× MIC vs control and 1/2× MIC vs control — `scripts/run_deseq2.py`
5. Figures: PCA, volcano plots (`scripts/pca_volcano.py`), heatmap of top
   differentially expressed genes (`scripts/heatmap_only.py`)

## Key outputs (from the 2026-10-09 run)

- Reads: 9 paired-end samples, 9.3–11.7M reads/sample (SRR8664715–SRR8664723)
- Mapping rates: >70% per sample against ATCC 25922 (e.g. C3: 76.4%)
- Significant DEGs (padj < 0.05, |log2FC| ≥ 1): **329** at 1/10× MIC, **300** at 1/2× MIC vs control
- PCA: PC1 explains 59.7% of variance; triplicates cluster by condition
- Full tables: `results/deseq2/` · Figures: `results/figures/` · QC: `results/qc/multiqc_report.html`

## Reproduce it

```bash
git clone <this-repo> ~/rnaseq-olaquindox
cd ~/rnaseq-olaquindox
conda env create -f environment.yml
conda activate rnaseq-ola
bash scripts/download.sh     # Stage 1
bash scripts/qc.sh           # Stage 2
bash scripts/quantify.sh     # Stage 3
python scripts/run_deseq2.py # Stage 4
python scripts/make_figures.py  # Stage 5
```

Or follow `docs/runbook.md`, which walks through the same commands with
pass/fail checkpoints at every stage.

## What this demonstrates

End-to-end ownership of a standard differential-expression workflow:
public-data wrangling, strain-aware reference selection, QC judgement,
count-based statistics in Python, and figure generation — with a pass/fail
checkpoint at every stage and no step taken on trust.

## Limitations

- Reanalysis of public data (not newly generated sequencing); original
  experimental design and any batch structure are inherited as published.
- Enrichment analysis depends on gene-set library availability for *E. coli*;
  see Stage 6 note in docs/runbook.md.
