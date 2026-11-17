"""
Chapter 3, Hands-on 4 - Does assembly quality change the typing result? Public Salmonella data.
Two typing outcomes are examined: whether an isolate is placed in a SNP cluster by the Pathogen Detection pipeline, and the number of AMR genes detected in its assembly. Each is related to assembly
quality (number of contigs, N50 and the 'Moderate' rule of Hands-on 1), crude and then adjusted for sequencing laboratory, serovar and collection year, because laboratory and serovar are tied to both quality
and outcome (Chapter 2). A random sample of 100,000 Illumina isolates with a known laboratory is used. Output: table_quality_typing.csv, table_quality_deciles.csv, Figure 4.
"""
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
from asm3 import load_asm
d = load_asm(); d = d[(d.plat == "ILLUMINA") & d.lab.notna() & d.serovar.notna() & d.year.between(2012, 2024) & d.asm_stats_n_contig.notna() & d.asm_stats_contig_n50.notna() & d.number_amr_genes.notna()].copy()
d["lc"] = np.log10(d.asm_stats_n_contig); d["ln50"] = np.log10(d.asm_stats_contig_n50 / 1e3); d["len_mb"] = d.asm_stats_length_bp / 1e6
d["serovar"] = d.serovar.str.lower()
d["fail"] = ((d.asm_stats_n_contig > 300) | (d.asm_stats_contig_n50 < 50e3)).astype(int)
for c, k in [("lab", 40), ("serovar", 20)]: top = d[c].value_counts().head(k).index; d[c + "_g"] = d[c].where(d[c].isin(top), "other")
n_all = len(d); d = d.sample(100_000, random_state=7); d["yr"] = d.year.astype(int).astype(str)
print(f"{n_all:,} Illumina isolates with laboratory, serovar, year and AMR count; sample of {len(d):,}; in a SNP cluster {d.in_cluster.mean():.1%}; failing the Moderate rule {d.fail.mean():.1%}; AMR genes mean {d.number_amr_genes.mean():.2f}")
rows = []
def logit(f, term, label, outcome):
    m = smf.glm(f, d, family=sm.families.Binomial()).fit(); b = m.params[term]; se = m.bse[term]; return dict(outcome=outcome, exposure=label, adjustment="adjusted" if "C(" in f else "crude", estimate=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), kind="odds ratio")
def pois(f, term, label, outcome):
    m = smf.glm(f, d, family=sm.families.Poisson()).fit(cov_type="HC1"); b = m.params[term]; se = m.bse[term]; return dict(outcome=outcome, exposure=label, adjustment="adjusted" if "C(" in f else "crude", estimate=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), kind="rate ratio")
adj = " + C(lab_g) + C(serovar_g) + C(yr)"
for f, term, label in [("in_cluster ~ lc", "lc", "10-fold more contigs"), ("in_cluster ~ lc" + adj, "lc", "10-fold more contigs"), ("in_cluster ~ fail", "fail", "fails Moderate rule"), ("in_cluster ~ fail" + adj, "fail", "fails Moderate rule")]: rows.append(logit(f, term, label, "In a SNP cluster"))
for f, term, label in [("number_amr_genes ~ lc", "lc", "10-fold more contigs"), ("number_amr_genes ~ lc" + adj, "lc", "10-fold more contigs"), ("number_amr_genes ~ fail", "fail", "fails Moderate rule"), ("number_amr_genes ~ fail" + adj, "fail", "fails Moderate rule")]: rows.append(pois(f, term, label, "AMR genes detected"))
R = pd.DataFrame(rows); R.round(4).to_csv("table_quality_typing.csv", index=False); print(R.round(3).to_string(index=False))
d["dec"] = pd.qcut(d.asm_stats_n_contig, 10, duplicates="drop"); D = d.groupby("dec", observed=True).agg(contigs=("asm_stats_n_contig", "median"), isolates=("lc", "size"), in_cluster=("in_cluster", "mean"), amr=("number_amr_genes", "mean"), fail=("fail", "mean")); D["in_cluster"] *= 100; D.round(3).to_csv("table_quality_deciles.csv"); print(D.round(2).to_string())
# within-laboratory view: laboratories with both good and poor assemblies
lab = d.groupby("lab_g").apply(lambda g: pd.Series(dict(n=len(g), clu_pass=100 * g[g.fail == 0].in_cluster.mean(), clu_fail=100 * g[g.fail == 1].in_cluster.mean() if (g.fail == 1).sum() >= 100 else np.nan, nfail=(g.fail == 1).sum()))).dropna(); print(f"{len(lab)} laboratories with at least 100 failing isolates; in-cluster share among passing isolates {lab.clu_pass.median():.1f}% (median) against {lab.clu_fail.median():.1f}% among failing; failing lower in {int((lab.clu_fail < lab.clu_pass).sum())} of {len(lab)}")
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6))
ax[0].plot(range(len(D)), D.in_cluster, "o-", color="#2b6cb0"); ax[0].set_xticks(range(len(D))); ax[0].set_xticklabels([f"{int(c)}" for c in D.contigs], fontsize=8); ax[0].set(xlabel="Median contigs in the decile", ylabel="Isolates placed in a SNP cluster (%)", title="Cluster membership by contig count")
ax[1].plot(range(len(D)), D.amr, "o-", color="#dd6b20"); ax[1].set_xticks(range(len(D))); ax[1].set_xticklabels([f"{int(c)}" for c in D.contigs], fontsize=8); ax[1].set(xlabel="Median contigs in the decile", ylabel="AMR genes detected per isolate", title="AMR genes by contig count")
y0 = np.arange(len(R))[::-1]; col = {"crude": "#a0aec0", "adjusted": "#2b6cb0"}
for i, r in enumerate(R.itertuples()): ax[2].errorbar(r.estimate, y0[i], xerr=[[r.estimate - r.lo], [r.hi - r.estimate]], fmt="o", color=col[r.adjustment], capsize=3)
ax[2].axvline(1, color="black", ls="--", lw=0.8); ax[2].set_yticks(y0); ax[2].set_yticklabels([f"{r.outcome[:12]}: {r.exposure[:18]} ({r.adjustment[:3]})" for r in R.itertuples()], fontsize=7); ax[2].set(xlabel="Odds ratio or rate ratio", title="Crude (grey) and adjusted (blue)")
fig.tight_layout(); fig.savefig("figs3/fig4_quality_typing.png", dpi=200); plt.close(fig); print("figure written")
