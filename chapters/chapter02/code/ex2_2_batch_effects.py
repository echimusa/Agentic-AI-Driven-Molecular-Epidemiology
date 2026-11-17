"""
Chapter 2, Hands-on 2 - Batch effects: detection, consequences for design, and a real-data example.
Part 1 repeats the uploaded R simulation (100 samples, 20 assays, a 1.5-unit shift in batch 2) in Python and asks what happens when case status is spread across batches in three ways:
  randomised (half the cases in each batch), partly confounded (80% of cases in batch 2) and fully confounded (all cases in batch 2). The true case effect is zero, so every 'hit' is a false positive.
Part 2 measures how much of the variation in assembly quality of public Salmonella Enteritidis genomes is explained by the sequencing laboratory, and what that does to a comparison of human and non-human isolates.
Output: table_batch_design.csv, table_batch_labs.csv, Figure 2.
"""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(123); NS, NM, SHIFT, REPS = 100, 20, 1.5, 500; batch = np.r_[np.zeros(50, int), np.ones(50, int)]

def one(design, effect):
    X = rng.normal(5, 1, (NS, NM)) + SHIFT * batch[:, None]
    n2 = {"Randomised": 25, "Partly confounded": 40, "Fully confounded": 50}[design]
    case = np.zeros(NS, int); case[rng.choice(np.where(batch == 1)[0], n2, replace=False)] = 1; case[rng.choice(np.where(batch == 0)[0], 50 - n2, replace=False)] = 1
    X = X + effect * case[:, None]
    naive = stats.ttest_ind(X[case == 1], X[case == 0], axis=0).pvalue
    D = np.column_stack([np.ones(NS), case, batch]); corr = np.full(NM, np.nan)                  # adjust for batch as a covariate in a regression of each assay on case status
    if np.linalg.matrix_rank(D) == 3:
        XtXi = np.linalg.inv(D.T @ D); B = XtXi @ D.T @ X; res = X - D @ B; s2 = (res ** 2).sum(axis=0) / (NS - 3); tstat = B[1] / np.sqrt(s2 * XtXi[1, 1]); corr = 2 * stats.t.sf(np.abs(tstat), NS - 3)
    return naive, corr

rows = []
for design in ["Randomised", "Partly confounded", "Fully confounded"]:
    for effect in (0.0, 0.6):
        pn, pc = zip(*[one(design, effect) for _ in range(REPS)]); pn, pc = np.array(pn), np.array(pc)
        rows.append(dict(design=design, true_effect=effect, naive_hit_rate=(pn < 0.05).mean() * 100, naive_any_of_20=(pn < 0.05).any(axis=1).mean() * 100, corrected_hit_rate=np.nan if np.isnan(pc).all() else (pc < 0.05).mean() * 100, corrected_any_of_20=np.nan if np.isnan(pc).all() else (pc < 0.05).any(axis=1).mean() * 100))
T = pd.DataFrame(rows); T.round(2).to_csv("table_batch_design.csv", index=False); print(T.round(1).to_string(index=False))

# PCA picture of the uploaded simulation
X0 = rng.normal(5, 1, (NS, NM)) + SHIFT * batch[:, None]; Z = (X0 - X0.mean(0)) / X0.std(0, ddof=1); U, S, Vt = np.linalg.svd(Z, full_matrices=False); pcs = U * S; var = S ** 2 / (S ** 2).sum()
print(f"Python PCA of the uploaded design: PC1 {var[0] * 100:.1f}% and PC2 {var[1] * 100:.1f}% of the variance; PC1 against batch: t-test p = {stats.ttest_ind(pcs[batch == 0, 0], pcs[batch == 1, 0]).pvalue:.2g}")

# Part 2: real data, sequencing laboratory as the batch
import os
d = pd.read_pickle("asm_batch.pkl") if os.path.exists("asm_batch.pkl") else pd.read_csv("asm_batch_enteritidis.csv.gz", dtype={"country": str}); d["lab"] = d.sequenced_by.fillna(d.bioproject_center)
x = d[(d.plat == "ILLUMINA") & (d.serovar.str.lower() == "enteritidis") & d.lab.notna() & d.asm_stats_n_contig.notna() & d.asm_stats_contig_n50.notna() & d.sector.isin(["Human", "Food", "Animal", "Environment"])].copy()
big = x.lab.value_counts(); big = big[big >= 300].index; xb = x[x.lab.isin(big)].copy(); xb["lc"] = np.log10(xb.asm_stats_n_contig); xb["ln50"] = np.log10(xb.asm_stats_contig_n50); xb["len_mb"] = xb.asm_stats_length_bp / 1e6
print(f"{len(xb):,} Enteritidis Illumina isolates from {len(big)} laboratories with at least 300")
eta = {}
for c, lab in [("lc", "log10 contigs"), ("ln50", "log10 N50"), ("len_mb", "genome length"), ("number_amr_genes", "AMR genes")]:
    g = xb.groupby("lab")[c]; eta[lab] = (g.count() * (g.mean() - xb[c].mean()) ** 2).sum() / ((xb[c] - xb[c].mean()) ** 2).sum()
print("share of variance between laboratories:", {k: round(v, 3) for k, v in eta.items()})
x["lc"] = np.log10(x.asm_stats_n_contig); x["group"] = np.where(x.sector == "Human", "Human", "Non-human"); ct = x.groupby(["lab", "group"]).size().unstack(fill_value=0); both = ct[(ct.Human >= 30) & (ct["Non-human"] >= 30)]
m0 = smf.ols("lc ~ C(group)", x).fit(); b0 = m0.params["C(group)[T.Non-human]"]; xs = x[x.lab.isin(both.index)]; m1 = smf.ols("lc ~ C(group)", xs).fit(); m2 = smf.ols("lc ~ C(group) + C(lab)", xs).fit()
share_conf = x[x.group == "Non-human"].lab.map(lambda l: ct.loc[l, "Human"] == 0).mean()
print(f"all labs: crude difference in log10 contigs (non-human minus human) {b0:.3f}, ratio {10 ** b0:.2f}; {len(both)} of {len(ct)} laboratories sequenced at least 30 human and 30 non-human isolates; in those: crude {m1.params['C(group)[T.Non-human]']:.3f}, within-laboratory {m2.params['C(group)[T.Non-human]']:.3f} (ratio {10 ** m2.params['C(group)[T.Non-human]']:.2f}); {share_conf:.1%} of non-human isolates come from laboratories with no human isolates")
L = xb.groupby("lab").agg(isolates=("lc", "size"), median_contigs=("asm_stats_n_contig", "median"), median_n50_kb=("asm_stats_contig_n50", lambda s: s.median() / 1000), human_share=("sector", lambda s: (s == "Human").mean())).sort_values("isolates", ascending=False); L.head(15).round(3).to_csv("table_batch_labs.csv")

fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6), gridspec_kw={"width_ratios": [0.9, 1.1, 1.2]})
for b, c, lab in [(0, "#2b6cb0", "Batch 1"), (1, "#dd6b20", "Batch 2")]: ax[0].scatter(pcs[batch == b, 0], pcs[batch == b, 1], s=22, color=c, alpha=0.8, label=lab)
ax[0].set(xlabel=f"PC1 ({var[0] * 100:.1f}% of variance)", ylabel=f"PC2 ({var[1] * 100:.1f}%)", title="PCA exposes the batch (uploaded design)"); ax[0].legend(frameon=False, fontsize=8.5)
dd = T[T.true_effect == 0]; xs_ = np.arange(3); ax[1].bar(xs_ - 0.2, dd.naive_hit_rate, 0.4, color="#c53030", label="Naive test"); ax[1].bar(xs_ + 0.2, dd.corrected_hit_rate.fillna(0), 0.4, color="#2b6cb0", label="Adjusted for batch"); ax[1].axhline(5, color="black", ls="--", lw=1)
ax[1].text(2.2, 12, "cannot be\nestimated", ha="center", fontsize=8, color="#2b6cb0"); ax[1].set_xticks(xs_); ax[1].set_xticklabels(["Randomised", "Partly\nconfounded", "Fully\nconfounded"], fontsize=9); ax[1].set(ylabel="False positives per assay at P < 0.05 (%)", title="Only a balanced design survives a batch effect"); ax[1].legend(frameon=False, fontsize=8.3, loc="upper left")
top = L.head(15)[::-1]; sc = ax[2].scatter(top.median_contigs, range(len(top)), s=np.sqrt(top.isolates) * 2.2, c=top.human_share, cmap="viridis_r", vmin=0, vmax=1); ax[2].set_yticks(range(len(top))); short = {"Public Health Agency of Canada": "PHAC (Canada)", "United States Department of Agriculture|Food Safety and Inspection Service": "USDA FSIS", "Microbiological Diagnostic Unit Public Health Laboratory": "MDU PHL (Australia)", "New York State Department of Health": "NY State DoH"}; ax[2].set_yticklabels([short.get(s, s)[:24] for s in top.index], fontsize=7.5)
ax[2].set(xlabel="Median contigs per assembly (S. Enteritidis)", title=f"Laboratory explains {eta['log10 contigs'] * 100:.0f}% of contig-count variance"); cb = fig.colorbar(sc, ax=ax[2], pad=0.01); cb.set_label("Share of human isolates", fontsize=8.5)
fig.tight_layout(); fig.savefig("figs2/fig2_batch_effects.png", dpi=200); plt.close(fig); print("figure written")
