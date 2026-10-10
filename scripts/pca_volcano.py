from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
BASE = Path.home() / "rnaseq-olaquindox"
FIG = BASE / "results/figures"
sns.set_theme(style="whitegrid", context="notebook")
samples = pd.read_csv(BASE / "data/samples.tsv", sep="\t", index_col="sample")
counts = pd.read_csv(BASE / "results/deseq2/counts_matrix.csv", index_col=0)
lc = np.log2(counts.div(counts.sum(axis=0), axis=1) * np.median(counts.sum(axis=0)) + 1)
# PCA
_f = lc.replace([np.inf, -np.inf], np.nan).dropna()
_f = _f.loc[_f.std(axis=1) > 0]
c = (_f.sub(_f.mean(axis=1), axis=0)).T
G = c.values @ c.values.T
ev, ec = np.linalg.eigh(G)
o = np.argsort(ev)[::-1]
ev, ec = np.clip(ev[o], 0, None), ec[:, o]
ve = ev / ev.sum() * 100
sc = ec * np.sqrt(ev)
pca = pd.DataFrame({"PC1": sc[:, 0], "PC2": sc[:, 1]}, index=c.index)
pca["condition"] = samples.loc[pca.index, "condition"]
plt.figure(figsize=(7, 5))
sns.scatterplot(data=pca, x="PC1", y="PC2", hue="condition", s=110)
for n, r in pca.iterrows():
    plt.annotate(n, (r["PC1"], r["PC2"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
plt.xlabel(f"PC1 ({ve[0]:.1f}%)"); plt.ylabel(f"PC2 ({ve[1]:.1f}%)")
plt.title("PCA — PRJNA525186"); plt.tight_layout()
plt.savefig(FIG / "pca.png", dpi=200); plt.close()
print(f"pca done (PC1={ve[0]:.1f}% PC2={ve[1]:.1f}%)", flush=True)
# Volcano
for csv, title, out in [("deseq2_low_vs_control.csv","1/10x MIC olaquindox vs control","volcano_low_vs_control.png"),("deseq2_high_vs_control.csv","1/2x MIC olaquindox vs control","volcano_high_vs_control.png")]:
    res = pd.read_csv(BASE / f"results/deseq2/{csv}", index_col=0).dropna(subset=["padj","log2FoldChange"]).copy()
    res["nlp"] = -np.log10(res["padj"].clip(lower=1e-300))
    sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)
    plt.figure(figsize=(7, 5))
    plt.scatter(res.loc[~sig,"log2FoldChange"], res.loc[~sig,"nlp"], s=6, c="grey", alpha=0.5, label="n.s.")
    plt.scatter(res.loc[sig,"log2FoldChange"], res.loc[sig,"nlp"], s=8, c="crimson", alpha=0.7, label="sig")
    plt.axhline(-np.log10(0.05), ls="--", lw=0.8, c="black")
    plt.axvline(1, ls="--", lw=0.8, c="black"); plt.axvline(-1, ls="--", lw=0.8, c="black")
    plt.xlabel("log2 fold change"); plt.ylabel("-log10 padj"); plt.title(title); plt.legend()
    plt.tight_layout(); plt.savefig(FIG / out, dpi=200); plt.close()
    print(f"{out} done ({int(sig.sum())} sig)", flush=True)
