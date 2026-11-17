"""
Chapter 1, Hands-on 1 - Discrimination and agreement between typing methods on public Salmonella data.
Data: NCBI Pathogen Detection Salmonella table (typing_all.csv.gz built from the metadata and cluster list used in Chapters 9, 12 and 13), restricted to the 4,811 isolates that carry a
serotype, a PFGE XbaI pattern and a whole-genome SNP cluster. Output: table_discrimination.csv, table_wallace.csv, Figure 4 of the chapter (figs1/fig1_discrimination.png).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from typing1 import simpson_D, wallace, subsample_ci
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
import os
SRC = "typing_all.csv.gz" if os.path.exists("typing_all.csv.gz") else "typing_pfge_subset.csv.gz"
d = pd.read_csv(SRC, dtype={"country": str, "pfge_xbai": str, "pfge_blni": str}, low_memory=False)
q = d[d.pfge_xbai.notna() & d.snp_cluster.notna() & d.serotype.notna()].copy(); print(f"{len(q):,} isolates with serotype, PFGE XbaI pattern and SNP cluster; sectors {q.sector.value_counts().to_dict()}")
levels = {"Serotype": "serotype", "PFGE (XbaI)": "pfge_xbai", "SNP cluster": "snp_cluster"}; rows = []
for name, c in levels.items():
    lab = q[c]; vc = lab.value_counts(); D = simpson_D(lab); lo, hi = subsample_ci(q[[c]], lambda x: simpson_D(x[c]))
    rows.append(dict(level=name, types=lab.nunique(), singleton_types_pct=(vc == 1).mean() * 100, isolates_in_singletons_pct=(vc[vc == 1].sum() / len(q)) * 100, largest_type_pct=vc.iloc[0] / len(q) * 100, D=D, D_lo=lo, D_hi=hi))
T = pd.DataFrame(rows); T.round(4).to_csv("table_discrimination.csv", index=False); print(T.round(4).to_string(index=False))
names = list(levels); W = pd.DataFrame(index=names, columns=names, dtype=float); AW = W.copy(); WL = W.copy(); WH = W.copy()
for a in names:
    for b in names:
        if a == b: continue
        w, aw = wallace(q[levels[a]], q[levels[b]]); W.loc[a, b] = w; AW.loc[a, b] = aw
        lo, hi = subsample_ci(q[[levels[a], levels[b]]], lambda x: wallace(x[levels[a]], x[levels[b]])[0]); WL.loc[a, b], WH.loc[a, b] = lo, hi
out = pd.concat({"W": W, "W_lo": WL, "W_hi": WH, "adjustedW": AW}, axis=1); out.round(3).to_csv("table_wallace.csv"); print("Wallace W(row -> column):"); print(W.round(3).to_string()); print("adjusted Wallace:"); print(AW.round(3).to_string())
for a, b in [("Serotype", "PFGE (XbaI)"), ("PFGE (XbaI)", "Serotype"), ("PFGE (XbaI)", "SNP cluster"), ("SNP cluster", "PFGE (XbaI)"), ("SNP cluster", "Serotype"), ("Serotype", "SNP cluster")]: print(f"W({a} -> {b}) = {W.loc[a, b]:.3f} (95% {WL.loc[a, b]:.3f}-{WH.loc[a, b]:.3f})")
top = q.pfge_xbai.value_counts().head(10).index; nest = q[q.pfge_xbai.isin(top)].groupby("pfge_xbai").agg(isolates=("target_acc", "size"), snp_clusters=("snp_cluster", "nunique"), serotypes=("serotype", "nunique")).loc[top]; print(nest.to_string())
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5), gridspec_kw={"width_ratios": [0.9, 0.9, 1.2]}); cols = ["#a0aec0", "#2b6cb0", "#c53030"]
ax[0].bar(T.level, T.D, yerr=[T.D - T.D_lo, T.D_hi - T.D], color=cols, capsize=4); ax[0].set(ylim=(0.8, 1.003), ylabel="Simpson's index of diversity (D)", title=f"Discrimination ({len(q):,} isolates)"); ax[0].tick_params(axis="x", labelsize=9)
for i, r in enumerate(T.itertuples()): ax[0].text(i, r.D_hi + 0.002, f"{r.types} types", ha="center", fontsize=8.5)
im = ax[1].imshow(W.values.astype(float), vmin=0, vmax=1, cmap="Blues"); ax[1].set_xticks(range(3)); ax[1].set_yticks(range(3)); ax[1].set_xticklabels(["Serotype", "PFGE", "SNP\ncluster"], fontsize=9); ax[1].set_yticklabels(["Serotype", "PFGE", "SNP cluster"], fontsize=9); ax[1].set(xlabel="... also share this type", ylabel="Pairs that share this type", title="Wallace coefficient")
for i in range(3):
    for j in range(3):
        if i != j: ax[1].text(j, i, f"{W.values[i, j]:.2f}", ha="center", va="center", color="white" if W.values[i, j] > 0.55 else "black", fontsize=10)
ax[2].barh(range(len(nest)), nest.snp_clusters, color="#2b6cb0"); ax[2].set_yticks(range(len(nest))); ax[2].set_yticklabels([f"{p} (n={n})" for p, n in zip(nest.index, nest.isolates)], fontsize=8); ax[2].invert_yaxis(); ax[2].set(xlabel="SNP clusters within the PFGE pattern", title="One gel pattern, many genomic clusters")
fig.tight_layout(); fig.savefig("figs1/fig1_discrimination.png", dpi=200); plt.close(fig); print("figure written")
