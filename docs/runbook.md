# RNA-seq Runbook: E. coli ATCC 25922 Olaquindox Response (PRJNA525186)

A laptop-friendly, end-to-end reanalysis of BioProject **PRJNA525186**: *Escherichia coli* ATCC 25922 exposed to sub-inhibitory olaquindox.

- **Conditions (9 samples, biological triplicates):** untreated control (B1–B3), 1/10× MIC olaquindox (D1–D3), 1/2× MIC olaquindox (C1–C3)
- **Workflow:** SRA download → FastQC/MultiQC → Salmon quantification → PyDESeq2 differential expression → PCA / volcano / heatmap figures → optional gseapy enrichment
- **How to use this runbook:** Run stages in order. Do not start the next stage until you observe the ✅ PASS checkpoint at the end of the current stage. Every command is copy-paste ready for Linux (or WSL2 Ubuntu on Windows).
- **Time/disk expectations:** ~20–40 GB free disk; downloads and Salmon quantification are the slowest steps on a laptop.

---

## Stage 0 — Laptop environment setup

### Windows users: WSL2 first (2 minutes of reading)

On Windows, install WSL2 with Ubuntu so every command below is identical Linux:

```powershell
# In PowerShell (as Administrator), then reboot when prompted
wsl --install -d Ubuntu
```

After reboot, open the **Ubuntu** app, create your Linux user, and run *all* subsequent commands inside Ubuntu — never in PowerShell/CMD. Keep project files inside the Linux filesystem (e.g. `~/rnaseq-olaquindox`), not under `/mnt/c/`, for much faster I/O.

### Install Miniconda (Linux / WSL2 / macOS)

```bash
cd ~
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
exec bash
conda --version
```

> macOS (Apple Silicon): replace the installer name with `Miniconda3-latest-MacOSX-arm64.sh`.

### Create the one named environment

```bash
conda create -n rnaseq-ola -c conda-forge -c bioconda -y \
  python=3.11 \
  sra-tools entrez-direct fastqc multiqc salmon gffread ncbi-datasets-cli \
  pip
conda activate rnaseq-ola
pip install pydeseq2 pandas numpy matplotlib seaborn gseapy
```

### Confirm the tools resolve

```bash
which fasterq-dump fastqc multiqc salmon gffread datasets
python -c "import pydeseq2, pandas, matplotlib, seaborn, gseapy; print('python imports OK')"
salmon --version; fastqc --version; multiqc --version
```

**✅ PASS checkpoint:** `conda activate rnaseq-ola` succeeds, all `which` paths point inside `~/miniconda3/envs/rnaseq-ola/`, and the Python import line prints `python imports OK`. If any tool is missing, re-run the `conda create` command — do not proceed with a partial environment.

---

## Stage 1 — Resolve and download the 9 runs for PRJNA525186

### Create the project layout and fetch run metadata

```bash
mkdir -p ~/rnaseq-olaquindox/{data/fastq,reference,results/{qc,salmon,deseq2,figures,enrichment},scripts}
cd ~/rnaseq-olaquindox

# Resolve the SRA study and its 9 experiment/run records
esearch -db sra -query "PRJNA525186" | efetch -format runinfo > data/SraRunInfo.csv
cat data/SraRunInfo.csv
cut -d, -f1 data/SraRunInfo.csv | grep ^SRR | sort > data/srr_list.txt
wc -l data/srr_list.txt
cat data/srr_list.txt
```

### The 9 verified run accessions

Resolved live from NCBI E-utilities for PRJNA525186 (study SRP187444). Sample naming follows the original study: **B = untreated, D = 1/10× MIC, C = 1/2× MIC**.

| Sample | Condition | Run accession |
|---|---|---|
| B1 | untreated control | SRR8664720 |
| B2 | untreated control | SRR8664719 |
| B3 | untreated control | SRR8664722 |
| D1 | 1/10× MIC olaquindox | SRR8664718 |
| D2 | 1/10× MIC olaquindox | SRR8664717 |
| D3 | 1/10× MIC olaquindox | SRR8664723 |
| C1 | 1/2× MIC olaquindox | SRR8664721 |
| C2 | 1/2× MIC olaquindox | SRR8664716 |
| C3 | 1/2× MIC olaquindox | SRR8664715 |

> If your `data/srr_list.txt` differs from this table, trust your live `SraRunInfo.csv` and rebuild the sample sheet in Stage 4 to match it — never hand-edit accession numbers from memory.

### Write the sample sheet now (you will reuse it in Stages 4–5)

```bash
cat > data/samples.tsv <<'EOF'
sample	condition	run_accession
B1	control	SRR8664720
B2	control	SRR8664719
B3	control	SRR8664722
D1	low_1_10MIC	SRR8664718
D2	low_1_10MIC	SRR8664717
D3	low_1_10MIC	SRR8664723
C1	high_1_2MIC	SRR8664721
C2	high_1_2MIC	SRR8664716
C3	high_1_2MIC	SRR8664715
EOF
column -t data/samples.tsv
```

### Download with fasterq-dump

```bash
cd ~/rnaseq-olaquindox/data/fastq
while read -r srr; do
  echo "=== $srr ==="
  prefetch "$srr"
  fasterq-dump --split-files --threads 4 --outdir . "$srr"
  gzip -f "${srr}"*.fastq
done < ../srr_list.txt
ls -lh
```

Count reads per file (plausibility check — expect millions of reads per sample for this HiSeq X Ten dataset):

```bash
for f in *.fastq.gz; do echo -n "$f  "; zcat "$f" | awk 'END{print NR/4 " reads"}'; done
```

> **Layout note:** SRA metadata for these runs lists a single FASTQ per run. If `fasterq-dump --split-files` produces only `*_1.fastq.gz` (no `*_2.fastq.gz`), the runs are single-end — that is expected here, and Stage 3 includes the exact single-end Salmon command. If you *do* get `_1`/`_2` pairs, use the paired-end variant noted in Stage 3.

**✅ PASS checkpoint:** `data/srr_list.txt` contains exactly **9** SRR lines, you have **9 FASTQ sets** in `data/fastq/` (9 files single-end, or 18 files paired-end), every file is >100 MB, and each read count is in the millions (not thousands). Zero-byte or truncated files = FAIL: delete that run's files and re-run its `prefetch`/`fasterq-dump`.

---

## Stage 2 — QC with FastQC + MultiQC

```bash
cd ~/rnaseq-olaquindox
mkdir -p results/qc/fastqc
fastqc --threads 4 --outdir results/qc/fastqc data/fastq/*.fastq.gz
multiqc results/qc/fastqc --outdir results/qc --filename multiqc_report.html
ls -lh results/qc/multiqc_report.html
xdg-open results/qc/multiqc_report.html   # macOS: open ... ; WSL2: explorer.exe results/qc/multiqc_report.html
```

### How to judge the report (2 minutes)

- **Per-base sequence quality:** median Phred should sit above ~28 across most of the read; a tail drop-off in the last bases is normal for this data type.
- **Adapter content:** should be low (<5–10%). High adapter % means trimming would be needed before quantification — note it, but do not silently trim; Salmon is tolerant of modest adapter contamination.
- **Per-sequence GC content:** curves should be broadly similar across all 9 samples; one wildly different sample is a red flag for a swapped or contaminated library.
- **Sequence duplication:** moderate duplication is expected in bacterial RNA-seq (highly expressed transcripts); extreme duplication in a single replicate deserves a note in your lab notebook.

**✅ PASS checkpoint:** `results/qc/multiqc_report.html` opens in a browser, shows **all 9 samples** in every module, and no sample is a gross outlier in GC/duplication versus its two biological replicates. Any failed sample must be named explicitly before you continue — do not average a bad replicate away.

---

## Stage 3 — Salmon index and quantification

### Reference caveat — read this first

The biological strain is **E. coli ATCC 25922**, whose reference assembly is **GCF_000743255.1** (used by the original study). Quantifying against **K-12 MG1655** instead would force reads from ATCC 25922–specific genes onto a mismatched transcriptome and distort counts. **Use the ATCC 25922 assembly below.** MG1655 is acceptable only as a teaching comparison, and any such comparison must be labelled as a different-strain reference.

### Download the ATCC 25922 genome + annotation and build the transcriptome

```bash
cd ~/rnaseq-olaquindox/reference
datasets download genome accession GCF_000743255.1 --include genome,gff3 --filename atcc25922.zip
unzip -o atcc25922.zip
find . -name "*.fna" -o -name "*.gff" | head
GENOME=$(find . -name "GCF_000743255.1*.fna" | head -n1)
GFF=$(find . -name "genomic.gff" | head -n1)
echo "$GENOME"; echo "$GFF"
gffread "$GFF" -g "$GENOME" -w transcripts.fasta
grep -c "^>" transcripts.fasta
head -n 2 transcripts.fasta
```

### Build the Salmon index

```bash
salmon index \
  --transcripts transcripts.fasta \
  --index salmon_index \
  --kmerLen 31
ls -lh salmon_index/
```

### Quantify all 9 samples

Single-end (expected for this dataset per Stage 1's layout note):

```bash
cd ~/rnaseq-olaquindox
while IFS=$'\t' read -r sample condition srr; do
  [ "$sample" = "sample" ] && continue
  echo "=== $sample ($srr) ==="
  salmon quant \
    --index reference/salmon_index \
    --libType A \
    -r "data/fastq/${srr}_1.fastq.gz" \
    --seqBias --gcBias --validateMappings \
    --threads 4 \
    --output "results/salmon/${sample}"
done < data/samples.tsv
```

> **Paired-end variant** (only if Stage 1 produced `_2` files): replace the `-r ...` line with `-1 "data/fastq/${srr}_1.fastq.gz" -2 "data/fastq/${srr}_2.fastq.gz"`.

### Inspect mapping rates

```bash
grep -H "Mapping rate" results/salmon/*/logs/salmon_quant.log
ls results/salmon/*/quant.sf
head results/salmon/B1/quant.sf
```

**✅ PASS checkpoint:** All **9** `results/salmon/<sample>/quant.sf` files exist and are non-empty, and every mapping rate printed from the logs is reported and plausible for a same-strain bacterial RNA-seq experiment (typically well above 70%; consistently low rates suggest a wrong-strain reference — revisit the caveat above). Do not proceed to DE with a sample whose quant failed.

---

## Stage 4 — Differential expression with PyDESeq2

Design: `~ condition`, with `control` as the reference level. Contrasts: `low_1_10MIC vs control` and `high_1_2MIC vs control`.

Save as `scripts/run_deseq2.py`:

```python
#!/usr/bin/env python3
"""PyDESeq2 differential expression for PRJNA525186 (control / low / high olaquindox)."""
from pathlib import Path
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path.home() / "rnaseq-olaquindox"
SAMPLES = pd.read_csv(BASE / "data/samples.tsv", sep="\t", index_col="sample")
OUT = BASE / "results/deseq2"
OUT.mkdir(parents=True, exist_ok=True)

# --- Build the gene x sample count matrix from Salmon quant.sf files ---
count_frames = {}
for sample in SAMPLES.index:
    q = pd.read_csv(BASE / f"results/salmon/{sample}/quant.sf", sep="\t")
    # Salmon gives estimated counts; DESeq2 needs integers.
    count_frames[sample] = q.set_index("Name")["NumReads"].round().astype(int)

counts = pd.DataFrame(count_frames).fillna(0).astype(int)
counts = counts.loc[counts.sum(axis=1) > 0]  # drop all-zero transcripts
counts.to_csv(OUT / "counts_matrix.csv")
print(f"count matrix: {counts.shape[0]} transcripts x {counts.shape[1]} samples")

metadata = SAMPLES.loc[counts.columns, ["condition"]].copy()
metadata["condition"] = pd.Categorical(
    metadata["condition"], categories=["control", "low_1_10MIC", "high_1_2MIC"]
)

# --- Fit the model ---
dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition")
dds.deseq2()
dds.write_pickle(str(OUT / "dds.pkl"))
print(dds)

# --- Contrasts vs control ---
for level, label in [("low_1_10MIC", "low_vs_control"), ("high_1_2MIC", "high_vs_control")]:
    stats = DeseqStats(dds, contrast=["condition", level, "control"])
    stats.summary()
    res = stats.results_df.copy()
    res.to_csv(OUT / f"deseq2_{label}.csv")
    n_sig = int(((res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)).sum())
    print(f"{label}: {n_sig} genes with padj<0.05 and |log2FC|>=1  ->  deseq2_{label}.csv")

print("done")
```

Run it:

```bash
cd ~/rnaseq-olaquindox
python scripts/run_deseq2.py
ls -lh results/deseq2/
head results/deseq2/deseq2_low_vs_control.csv
head results/deseq2/deseq2_high_vs_control.csv
```

**✅ PASS checkpoint:** `results/deseq2/` contains `counts_matrix.csv`, `deseq2_low_vs_control.csv`, and `deseq2_high_vs_control.csv`; each results table visibly includes `log2FoldChange`, `pvalue`, and `padj` columns, has one row per tested transcript, and the script's printed significant-gene counts are non-negative integers you have written down. An empty or NaN-only table = FAIL: check that `condition` labels in `samples.tsv` match the script exactly.

---

## Stage 5 — Figures: PCA, volcano plots, top-DEG heatmap

Save as `scripts/make_figures.py` (one runnable script, no pseudocode):

```python
#!/usr/bin/env python3
"""PCA, volcano plots, and top-DEG heatmap from PyDESeq2 outputs."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path.home() / "rnaseq-olaquindox"
FIG = BASE / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="notebook")

samples = pd.read_csv(BASE / "data/samples.tsv", sep="\t", index_col="sample")
counts = pd.read_csv(BASE / "results/deseq2/counts_matrix.csv", index_col=0)

# --- Variance-stabilising transform (log of size-factor-naive normalised counts) ---
lib_sizes = counts.sum(axis=0)
normed = counts.div(lib_sizes, axis=1) * np.median(lib_sizes)
log_counts = np.log2(normed + 1)
log_counts.to_csv(FIG / "log2_normalised_counts.csv")

# --- PCA ---
centred = (log_counts - log_counts.mean(axis=1)).T  # samples x genes
u, s, _ = np.linalg.svd(centred, full_matrices=False)
var_explained = (s ** 2) / (s ** 2).sum() * 100
pca = pd.DataFrame({"PC1": u[:, 0] * s[0], "PC2": u[:, 1] * s[1]}, index=centred.index)
pca["condition"] = samples.loc[pca.index, "condition"]

plt.figure(figsize=(7, 5))
sns.scatterplot(data=pca, x="PC1", y="PC2", hue="condition", s=110)
for name, row in pca.iterrows():
    plt.annotate(name, (row["PC1"], row["PC2"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
plt.xlabel(f"PC1 ({var_explained[0]:.1f}% variance)")
plt.ylabel(f"PC2 ({var_explained[1]:.1f}% variance)")
plt.title("PCA — PRJNA525186 (Salmon counts, log2-normalised)")
plt.tight_layout()
plt.savefig(FIG / "pca.png", dpi=200)
plt.close()
print(f"PCA written; PC1={var_explained[0]:.1f}%  PC2={var_explained[1]:.1f}%")

# --- Volcano plots ---
def volcano(csv_path: Path, title: str, out_name: str) -> None:
    res = pd.read_csv(csv_path, index_col=0)
    res = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    res["neglog10padj"] = -np.log10(res["padj"].clip(lower=1e-300))
    sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)
    plt.figure(figsize=(7, 5))
    plt.scatter(res.loc[~sig, "log2FoldChange"], res.loc[~sig, "neglog10padj"], s=6, c="grey", alpha=0.5, label="n.s.")
    plt.scatter(res.loc[sig, "log2FoldChange"], res.loc[sig, "neglog10padj"], s=8, c="crimson", alpha=0.7, label="padj<0.05, |log2FC|>=1")
    plt.axhline(-np.log10(0.05), ls="--", lw=0.8, c="black")
    plt.axvline(1, ls="--", lw=0.8, c="black"); plt.axvline(-1, ls="--", lw=0.8, c="black")
    plt.xlabel("log2 fold change"); plt.ylabel("-log10 adjusted p-value")
    plt.title(title); plt.legend(frameon=True)
    plt.tight_layout(); plt.savefig(FIG / out_name, dpi=200); plt.close()
    print(f"{out_name}: {int(sig.sum())} significant genes plotted")

volcano(BASE / "results/deseq2/deseq2_low_vs_control.csv",
        "1/10x MIC olaquindox vs control", "volcano_low_vs_control.png")
volcano(BASE / "results/deseq2/deseq2_high_vs_control.csv",
        "1/2x MIC olaquindox vs control", "volcano_high_vs_control.png")

# --- Heatmap of top DEGs (by padj, pooled across both contrasts) ---
frames = []
for f in ["deseq2_low_vs_control.csv", "deseq2_high_vs_control.csv"]:
    frames.append(pd.read_csv(BASE / f"results/deseq2/{f}", index_col=0))
pooled_padj = pd.concat([fr["padj"] for fr in frames], axis=1).min(axis=1)
top_genes = pooled_padj.nsmallest(30).index.tolist()
heat = log_counts.loc[top_genes]
heat_z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, 1), axis=0)
col_colors = samples.loc[heat.columns, "condition"].map(
    {"control": "#4C72B0", "low_1_10MIC": "#DD8452", "high_1_2MIC": "#C44E52"})
sns.clustermap(heat_z, cmap="vlag", center=0, col_colors=col_colors,
               xticklabels=True, yticklabels=True, figsize=(9, 10))
plt.savefig(FIG / "heatmap_top30_degs.png", dpi=200)
plt.close()
print("heatmap written: heatmap_top30_degs.png (top 30 DEGs by min padj)")
print("all figures done")
```

Run it:

```bash
cd ~/rnaseq-olaquindox
python scripts/make_figures.py
ls -lh results/figures/
```

**✅ PASS checkpoint:** `results/figures/` contains `pca.png`, `volcano_low_vs_control.png`, `volcano_high_vs_control.png`, and `heatmap_top30_degs.png`; opening the PCA shows 9 labelled points, and the printed PC variance percentages and per-contrast significant-gene counts match what Stage 4 reported. Replicates that fail to group near each other on the PCA must be called out in writing before any biological interpretation.

---

## Stage 6 (optional) — Pathway enrichment with gseapy

gseapy's built-in libraries are human/mouse-centric; an *E. coli* analysis must first **verify at runtime** which libraries exist rather than assume a name. Run this inspection first:

```bash
cd ~/rnaseq-olaquindox
python - <<'EOF'
import gseapy as gp
libs = gp.get_library_name()
print(f"{len(libs)} libraries available")
for name in libs:
    low = name.lower()
    if any(k in low for k in ["kegg", "go", "ecoli", "bacteria", "prokary"]):
        print("candidate:", name)
EOF
```

Then run enrichment **only** against a library that the command above actually printed (substitute its exact name):

```bash
python - <<'EOF'
from pathlib import Path
import pandas as pd
import gseapy as gp

BASE = Path.home() / "rnaseq-olaquindox"
LIBRARY = "REPLACE_WITH_VERIFIED_LIBRARY_NAME"  # exact string printed by the inspection above
for label in ["low_vs_control", "high_vs_control"]:
    res = pd.read_csv(BASE / f"results/deseq2/deseq2_{label}.csv", index_col=0)
    genes = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)].index.tolist()
    print(label, "significant genes:", len(genes))
    if len(genes) < 10:
        print("  too few genes for enrichment; skipping")
        continue
    enr = gp.enrichr(gene_list=genes, gene_sets=LIBRARY, outdir=str(BASE / f"results/enrichment/{label}"))
    enr.results.to_csv(BASE / f"results/enrichment/{label}_enrichment.csv", index=False)
    print(enr.results.head())
EOF
ls -R results/enrichment/
```

> If the inspection prints **no** E. coli-compatible library, the honest result is "no suitable built-in library available" — record that sentence rather than forcing human gene sets onto bacterial locus tags. KEGG pathway mapping via gene-to-KO annotation is a legitimate follow-up, but it is out of scope for this runbook.

**✅ PASS checkpoint:** Either (a) `results/enrichment/` contains CSVs whose `Gene_set` column names the exact library printed by the inspection command, or (b) you have written the "no suitable library" note above. Claiming enriched pathways without (a) = FAIL.

---

## GitHub-ready repository layout

Assemble the finished project like this before pushing (large FASTQ files and the Salmon index stay out of git):

```text
ecoli-olaquindox-rnaseq/
├── README.md                  # the draft below, edited with your real results
├── environment.yml            # exported in the step below
├── .gitignore                 # data/fastq/, reference/salmon_index/, *.fastq.gz, __pycache__/
├── data/
│   ├── samples.tsv            # 9-row sample sheet (Stage 1)
│   └── SraRunInfo.csv         # NCBI run metadata as downloaded
├── scripts/
│   ├── run_deseq2.py          # Stage 4
│   └── make_figures.py        # Stage 5
├── results/
│   ├── qc/multiqc_report.html
│   ├── salmon/                # per-sample quant.sf + logs (small; keep)
│   ├── deseq2/                # counts_matrix.csv + both contrast tables
│   ├── figures/               # pca.png, volcano_*.png, heatmap_top30_degs.png
│   └── enrichment/            # only if Stage 6 produced real output
├── reference/
│   └── README.md              # assembly accession GCF_000743255.1 + how transcripts.fasta was built
└── docs/
    └── runbook.md             # this runbook, so the repo is self-documenting
```

Freeze the environment and initialise the repo:

```bash
cd ~/rnaseq-olaquindox
conda env export --name rnaseq-ola --file environment.yml
cat > .gitignore <<'EOF'
data/fastq/
reference/salmon_index/
reference/*.zip
reference/ncbi_dataset/
*.fastq.gz
__pycache__/
.ipynb_checkpoints/
EOF
git init
git add README.md environment.yml .gitignore data/samples.tsv data/SraRunInfo.csv scripts results reference/README.md docs
git commit -m "Reproducible RNA-seq reanalysis of PRJNA525186 (Salmon + PyDESeq2)"
git status
```

**✅ PASS checkpoint:** `git status` shows no `.fastq.gz` files staged, `environment.yml` exists, and every path in the layout above that you claim in the README actually exists on disk.

---

## README draft for the GitHub repository (hiring-manager facing)

Copy this into `ecoli-olaquindox-rnaseq/README.md` and replace only the bracketed fields with numbers from *your* run:

```markdown
# E. coli Olaquindox RNA-seq Reanalysis (PRJNA525186)

A fully reproducible reanalysis of a public bacterial RNA-seq dataset,
built end-to-end on a laptop: raw SRA reads -> quality control ->
transcript quantification -> differential expression -> publication-style figures.

## Biological question
How does *Escherichia coli* ATCC 25922 (a standard antimicrobial-susceptibility
reference strain) remodel its transcriptome when exposed to sub-inhibitory
concentrations of the veterinary antibiotic olaquindox?

## Dataset
- NCBI BioProject PRJNA525186 / SRA study SRP187444
- 9 samples, biological triplicates: untreated control (B1-B3),
  1/10x MIC olaquindox (D1-D3), 1/2x MIC olaquindox (C1-C3)
- Run accessions SRR8664715-SRR8664723 (see data/SraRunInfo.csv)

## Pipeline
1. Download: Entrez Direct + SRA Toolkit (fasterq-dump)
2. QC: FastQC, summarised with MultiQC
3. Quantification: Salmon (bias-corrected, validated mappings) against a
   transcriptome built with gffread from the strain-matched ATCC 25922
   assembly GCF_000743255.1 - not K-12 MG1655
4. Differential expression: PyDESeq2, design ~ condition, contrasts
   1/10x MIC vs control and 1/2x MIC vs control
5. Figures: PCA, volcano plots, heatmap of top differentially expressed genes
6. Optional: pathway enrichment with gseapy (library verified at runtime)

## Key outputs (from my run)
- Mapping rates: [report per-sample range from results/salmon/*/logs]
- Significant DEGs (padj < 0.05, |log2FC| >= 1): [N] at 1/10x MIC, [M] at 1/2x MIC
- PCA: PC1 explains [X]% of variance; triplicates cluster by condition
- Full tables: results/deseq2/ | Figures: results/figures/ | QC: results/qc/multiqc_report.html

## Reproduce it
    conda env create -f environment.yml
    conda activate rnaseq-ola
    # then follow docs/runbook.md stage by stage

## What this demonstrates
End-to-end ownership of a standard differential-expression workflow:
public-data wrangling, strain-aware reference selection, QC judgement,
count-based statistics in Python, and figure generation - with a pass/fail
checkpoint at every stage and no step taken on trust.

## Limitations
- Reanalysis of public data (not newly generated sequencing); original
  experimental design and any batch structure are inherited as published.
- Enrichment analysis depends on gene-set library availability for E. coli;
  see Stage 6 note in docs/runbook.md.
```

**Final ✅ PASS checkpoint:** Every bracketed number in the README came from a file in your `results/` directory that you opened and read — no placeholder brackets remain, and no result is described that the pipeline did not produce.

---

### Sources and provenance

- Dataset identity, design, and strain-matched assembly (GCF_000743255.1): the original study, "Development of Resistance in *Escherichia coli* ATCC25922 under Exposure of Sub-Inhibitory Concentration of Olaquindox" — data deposited under PRJNA525186.
- Run accessions (SRR8664715–SRR8664723), sample labels (B/C/D series), and study SRP187444: resolved live from NCBI SRA E-utilities (`esearch`/`efetch` on PRJNA525186) during preparation of this runbook; the Stage 1 commands reproduce them.
- Tool behaviour (fasterq-dump, Salmon bias/validation flags, gseapy library inspection) follows each tool's documented interface; Stage 6 deliberately verifies library names at runtime instead of assuming them.
