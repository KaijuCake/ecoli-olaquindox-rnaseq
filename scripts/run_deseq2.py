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
dds = DeseqDataSet(counts=counts.T, metadata=metadata, design_factors="condition")
dds.deseq2()
import pickle as _pickle; open(OUT / "dds.pkl", "wb").write(_pickle.dumps(dds))
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
