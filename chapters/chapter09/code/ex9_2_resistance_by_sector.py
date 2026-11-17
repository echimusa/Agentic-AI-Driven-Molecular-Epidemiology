"""
Chapter 9, Hands-on 2 - Resistance determinants and genotypic multidrug resistance (MDR) across One Health sectors.
Input: sal_us_frame.csv.gz from ex9_1 and the SNP-cluster list. Output: Figures 4 and 5, table_or.csv, frame_for_r.csv (for the R companion).
Genotypic MDR = core determinants for 3 or more drug classes. It is a genotype, not a phenotype.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm, statsmodels.formula.api as smf
from amr9 import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
order = ["Human", "Food", "Animal", "Environment"]; COL = {"Human": "#2b6cb0", "Food": "#dd6b20", "Animal": "#2f855a", "Environment": "#805ad5"}

f = pd.read_csv("sal_us_frame.csv.gz", dtype={"state": str}); f["genes"] = f.genes.fillna(""); f["classes"] = f.classes.fillna("")
cat = load_catalog(); ceph = {g for g, (c, s) in cat.items() if c == "BETA-LACTAM" and "CEPHALOSPORIN" in s}
f["ceph"] = f.genes.map(lambda s: int(any(g in ceph for g in s.split(";") if g)))             # a determinant for third-generation cephalosporin resistance (ESBL or AmpC)
cl = pd.read_csv("salmonella_cluster_list.tsv", sep="\t", usecols=["PDS_acc", "target_acc"])
f = f.merge(cl, on="target_acc", how="left"); f["grp"] = pd.factorize(f.PDS_acc.fillna("single_" + f.target_acc))[0]   # isolates in the same SNP cluster share a group
print(f"{len(f):,} isolates; {f.PDS_acc.notna().mean():.1%} belong to a SNP cluster; {f.grp.nunique():,} groups")

# ---------------------------------------------------------------- Figure 3: prevalence of each drug class by sector
cls = f.assign(c=f.classes.str.split(";")).explode("c"); cls = cls[cls.c.notna() & (cls.c != "")]
prev = (cls.groupby(["sector", "c"]).size().unstack(fill_value=0).reindex(order) .div(f.sector.value_counts().reindex(order), axis=0) * 100)
top = prev.max().sort_values(ascending=False).head(11).index.tolist(); prev = prev[top]; prev.round(1).T.to_csv("table_class_prevalence.csv"); print(prev.round(1).T.to_string())
fig, ax = plt.subplots(figsize=(7.6, 4.6)); im = ax.imshow(prev.T.values, cmap="Blues", vmin=0, vmax=max(40, prev.values.max()), aspect="auto")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.set_yticks(range(len(top))); ax.set_yticklabels([t.title() for t in top])
for i in range(len(top)):
    for j in range(len(order)): v = prev.T.values[i, j]; ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=9, color="white" if v > 25 else "black")
ax.set(title="Isolates with a determinant for each drug class (%)"); fig.colorbar(im, ax=ax, fraction=0.04); fig.tight_layout(); fig.savefig("figs9/fig4_class_prevalence.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- crude and adjusted odds of genotypic MDR, cluster-robust
topsv = f.serovar_clean.value_counts().head(25).index; f["serovar_grp"] = np.where(f.serovar_clean.isin(topsv), f.serovar_clean, "Other")
f[["target_acc", "sector", "serovar_grp", "year", "mdr_geno", "ceph", "grp"]].to_csv("frame_for_r.csv.gz", index=False)
fm = f[["mdr_geno", "sector", "serovar_grp", "year", "grp"]].copy()                  # the models need five columns; the full frame holds long text columns
def odds(formula, cluster=True):
    mod = smf.glm(formula, data=fm, family=sm.families.Binomial())
    r = mod.fit(cov_type="cluster", cov_kwds={"groups": fm.grp}) if cluster else mod.fit()
    ci = r.conf_int(); keep = [k for k in r.params.index if k.startswith("C(sector")]
    return pd.DataFrame({"term": [k.split("[T.")[1].rstrip("]") for k in keep], "OR": np.exp(r.params[keep]).values, "lo": np.exp(ci.loc[keep, 0]).values, "hi": np.exp(ci.loc[keep, 1]).values})
S = "C(sector, Treatment('Human'))"
crude = odds(f"mdr_geno ~ {S}"); adj = odds(f"mdr_geno ~ {S} + C(serovar_grp) + C(year)"); adj_plain = odds(f"mdr_geno ~ {S} + C(serovar_grp) + C(year)", cluster=False)
res = crude.merge(adj, on="term", suffixes=("_crude", "_adj")).merge(adj_plain.rename(columns={"OR": "OR_adj_plain", "lo": "lo_plain", "hi": "hi_plain"}), on="term"); res = res.set_index("term").loc[["Food", "Animal", "Environment"]]
res.round(3).to_csv("table_or.csv"); print(res.round(2).to_string())

# ---------------------------------------------------------------- Figure 4: forest plot and trends
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2), gridspec_kw={"width_ratios": [1.15, 1, 1]})
yy = np.arange(len(res))[::-1]
for k, (suffix, lab, off, c) in enumerate([("_crude", "Crude", 0.16, "#a0aec0"), ("_adj", "Adjusted for serovar and year", -0.16, "#2b6cb0")]):
    ax[0].errorbar(res[f"OR{suffix}"], yy + off, xerr=[res[f"OR{suffix}"] - res[f"lo{suffix}"], res[f"hi{suffix}"] - res[f"OR{suffix}"]], fmt="o", color=c, capsize=3, label=lab)
ax[0].axvline(1, color="black", lw=0.8, ls="--"); ax[0].set_xscale("log"); ax[0].set_xticks([0.5, 1, 2, 4, 8]); ax[0].get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter()); ax[0].minorticks_off(); ax[0].set_yticks(yy); ax[0].set_yticklabels(res.index); ax[0].set_xlabel("Odds ratio for genotypic MDR vs human (log scale)")
ax[0].set_title("Sector and genotypic MDR"); ax[0].legend(frameon=False, fontsize=8.3, loc="lower right")
for a, col, ttl in [(ax[1], "mdr_geno", "Genotypic MDR (%)"), (ax[2], "ceph", "Cephalosporin determinant (%)")]:
    t = f.groupby(["year", "sector"])[col].mean().unstack()[order] * 100; n = f.groupby(["year", "sector"]).size().unstack()[order]
    for s in order:
        ok = n[s] >= 100; a.plot(t.index[ok], t[s][ok], marker="o", ms=3.5, color=COL[s], label=s)
    a.set(title=ttl, xlabel="Collection year", xticks=range(2015, 2025, 3)); a.set_ylim(bottom=0)
ax[1].legend(frameon=False, fontsize=8.3, ncol=2); fig.tight_layout(); fig.savefig("figs9/fig5_mdr_forest_trends.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- facts used in the text
print("\ncrude genotypic MDR % by sector:", (f.groupby("sector").mdr_geno.mean().reindex(order) * 100).round(1).to_dict())
print("cephalosporin determinant % by sector:", (f.groupby("sector").ceph.mean().reindex(order) * 100).round(1).to_dict())
g = f.assign(g=f.genes.str.split(";")).explode("g"); g = g[g.g.notna() & (g.g != "")]
for s in order: print(s, "top determinants:", (g[g.sector == s].g.value_counts().head(5) / (f.sector == s).sum() * 100).round(1).to_dict())
print("CTX-M-65 % by sector:", {s: round(100 * f[f.sector == s].genes.str.contains("blaCTX-M-65").mean(), 2) for s in order})
