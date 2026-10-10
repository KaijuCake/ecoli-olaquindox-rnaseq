from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
BASE = Path.home() / "rnaseq-olaquindox"
FIG = BASE / "results/figures"
samples = pd.read_csv(BASE / "data/samples.tsv", sep="\t", index_col="sample")
counts = pd.read_csv(BASE / "results/deseq2/counts_matrix.csv", index_col=0)
lc = np.log2(counts.div(counts.sum(axis=0), axis=1) * np.median(counts.sum(axis=0)) + 1)
fr = [pd.read_csv(BASE / f"results/deseq2/{f}", index_col=0) for f in ["deseq2_low_vs_control.csv", "deseq2_high_vs_control.csv"]]
tg = pd.concat([x["padj"] for x in fr], axis=1).min(axis=1).nsmallest(30).index
print("n genes:", len(tg), flush=True)
print("missing:", [g for g in tg if g not in lc.index][:5], flush=True)
heat = lc.loc[[g for g in tg if g in lc.index]]
print("heat shape:", heat.shape, flush=True)
hz = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, 1), axis=0)
cc = samples.loc[hz.columns, "condition"].map({"control": "#4C72B0", "low_1_10MIC": "#DD8452", "high_1_2MIC": "#C44E52"})
sns.clustermap(hz, cmap="vlag", center=0, col_colors=cc, xticklabels=True, yticklabels=True, figsize=(9, 10))
plt.savefig(FIG / "heatmap_top30_degs.png", dpi=200); plt.close()
print("heatmap saved", flush=True)
