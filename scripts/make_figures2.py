#!/usr/bin/env python3
"""Figures, standalone (bypasses make_figures.py)."""
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
lib_sizes = counts.sum(axis=0)
log_counts = np.log2(counts.div(lib_sizes, axis=1) * np.median(lib_sizes) + 1)

# PCA via Gram eigendecomposition
_lc = log_counts.replace([np.inf, -np.inf], np.nan).dropna()
_lc = _lc.loc[_lc.std(axis=1) > 0]
centred = (_lc - _lc.mean(axis=1)).T
G = centred.values @ centred.values.T
eigvals, eigvecs = np.linalg.eigh(G)
o = np.argsort(eigvals)[::-1]
eigvals, eigvecs = np.clip(eigvals[o], 0, None), eigvecs[:, o]
var_exp = eigvals / eigvals.sum() * 100
scores = eigvecs * np.sqrt(eigvals)
pca = pd.DataFrame({"PC1": scores[:, 0], "PC2": scores[:, 1]}, index=centred.index)
pca["condition"] = samples.loc[pca.index, "condition"]
plt.figure(figsize=(7, 5))
sns.scatterplot(data=pca, x="PC1", y="PC2", hue="condition", s=110)
for n, r in pca.iterrows():
    plt.annotate(n, (r["PC1"], r["PC2"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
plt.xlabel(f"PC1 ({var_exp[0]:.1f}%)"); plt.ylabel(f"PC2 ({var_exp[1]:.1f}%)")
plt.title("PCA — PRJNA525186"); plt.tight_layout()
plt.savefig(FIG / "pca.png", dpi=200); plt.close()
print(f"pca.png done (PC1={var_exp[0]:.1f}% PC2={var_exp[1]:.1f}%)", flush=True)

# Volcano plots
for csv, title, out in [
    ("deseq2_low_vs_control.csv", "1/10x MIC olaquindox vs control", "volcano_low_vs_control.png"),
    ("deseq2_high_vs_control.csv", "1/2x MIC olaquindox vs control", "volcano_high_vs_control.png")]:
    res = pd.read_csv(BASE / f"results/deseq2/{csv}", index_col=0).dropna(subset=["padj", "log2FoldChange"]).copy()
    res["nlp"] = -np.log10(res["padj"].clip(lower=1e-300))
    sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)
    plt.figure(figsize=(7, 5))
    plt.scatter(res.loc[~sig, "log2FoldChange"], res.loc[~sig, "nlp"], s=6, c="grey", alpha=0.5, label="n.s.")
    plt.scatter(res.loc[sig, "log2FoldChange"], res.loc[sig, "nlp"], s=8, c="crimson", alpha=0.7, label="padj<0.05, |log2FC|>=1")
    plt.axhline(-np.log10(0.05), ls="--", lw=0.8, c="black")
    plt.axvline(1, ls="--", lw=0.8, c="black"); plt.axvline(-1, ls="--", lw=0.8, c="black")
    plt.xlabel("log2 fold change"); plt.ylabel("-log10 padj"); plt.title(title); plt.legend()
    plt.tight_layout(); plt.savefig(FIG / out, dpi=200); plt.close()
    print(f"{out} done ({int(sig.sum())} sig)", flush=True)

# Heatmap top 30
frames = [pd.read_csv(BASE / f"results/deseq2/{f}", index_col=0) for f in ["deseq2_low_vs_control.csv", "deseq2_high_vs_control.csv"]]
top_genes = pd.concat([fr["padj"] for fr in frames], axis=1).min(axis=1).nsmallest(30).index.tolist()
print("heatmap genes:", top_genes[:3], flush=True)
heat = log_counts.reindex(top_genes).dropna()
hz = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, 1), axis=0)
cc = samples.loc[hz.columns, "condition"].map({"control": "#4C72B0", "low_1_10MIC": "#DD8452", "high_1_2MIC": "#C44E52"})
sns.clustermap(hz, cmap="vlag", center=0, col_colors=cc, xticklabels=True, yticklabels=True, figsize=(9, 10))
plt.savefig(FIG / "heatmap_top30_degs.png", dpi=200); plt.close()
print("heatmap done", flush=True)
