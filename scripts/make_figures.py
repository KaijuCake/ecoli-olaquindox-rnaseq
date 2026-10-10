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

lib_sizes = counts.sum(axis=0)
normed = counts.div(lib_sizes, axis=1) * np.median(lib_sizes)
log_counts = np.log2(normed + 1)
log_counts.to_csv(FIG / "log2_normalised_counts.csv")

# --- PCA (drop NaN/inf and zero-variance genes — they break SVD) ---
_lc = log_counts.replace([np.inf, -np.inf], np.nan).dropna()
_lc = _lc.loc[_lc.std(axis=1) > 0]
centred = (_lc - _lc.mean(axis=1)).T  # samples x genes
G = centred.values @ centred.values.T
eigvals, eigvecs = np.linalg.eigh(G)
_order = np.argsort(eigvals)[::-1]
eigvals, eigvecs = eigvals[_order], eigvecs[:, _order]
eigvals = np.clip(eigvals, 0, None)
var_explained = eigvals / eigvals.sum() * 100
scores = eigvecs * np.sqrt(eigvals)
pca = pd.DataFrame({"PC1": scores[:, 0], "PC2": scores[:, 1]}, index=centred.index)
pca["condition"] = samples.loc[pca.index, "condition"]
plt.figure(figsize=(7, 5))
sns.scatterplot(data=pca, x="PC1", y="PC2", hue="condition", s=110)
for name, row in pca.iterrows():
    plt.annotate(name, (row["PC1"], row["PC2"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
plt.xlabel(f"PC1 ({var_explained[0]:.1f}% variance)")
plt.ylabel(f"PC2 ({var_explained[1]:.1f}% variance)")
plt.title("PCA — PRJNA525186 (Salmon counts, log2-normalised)")
plt.tight_layout(); plt.savefig(FIG / "pca.png", dpi=200); plt.close()
print(f"PCA written; PC1={var_explained[0]:.1f}%  PC2={var_explained[1]:.1f}%")

def volcano(csv_path, title, out_name):
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

volcano(BASE / "results/deseq2/deseq2_low_vs_control.csv", "1/10x MIC olaquindox vs control", "volcano_low_vs_control.png")
volcano(BASE / "results/deseq2/deseq2_high_vs_control.csv", "1/2x MIC olaquindox vs control", "volcano_high_vs_control.png")

frames = [pd.read_csv(BASE / f"results/deseq2/{f}", index_col=0) for f in ["deseq2_low_vs_control.csv", "deseq2_high_vs_control.csv"]]
pooled_padj = pd.concat([fr["padj"] for fr in frames], axis=1).min(axis=1)
top_genes = pooled_padj.nsmallest(30).index.tolist()
print("DEBUG top5:", top_genes[:5])
print("DEBUG n:", len(top_genes))
heat = log_counts.loc[top_genes]
heat_z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, 1), axis=0)
col_colors = samples.loc[heat.columns, "condition"].map({"control": "#4C72B0", "low_1_10MIC": "#DD8452", "high_1_2MIC": "#C44E52"})
sns.clustermap(heat_z, cmap="vlag", center=0, col_colors=col_colors, xticklabels=True, yticklabels=True, figsize=(9, 10))
plt.savefig(FIG / "heatmap_top30_degs.png", dpi=200); plt.close()
print("all figures done")
