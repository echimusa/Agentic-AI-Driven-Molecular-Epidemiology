"""
Chapter 9, Hands-on 1 - From public Salmonella isolates to a One Health sampling frame.
Downloads the NCBI Pathogen Detection Salmonella metadata (about 880 MB) if it is not already here, assigns each isolate to a sector,
builds the analysis frame (USA, 2015-2024) and draws Figures 1 and 2.
Needs: pandas numpy matplotlib scipy; the AMRFinderPlus reference gene catalogue (downloaded by ex9_5 or by this script).
"""
import os, urllib.request
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from amr9 import *

PD = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/latest_snps"
FILES = {"salmonella_metadata.tsv": f"{PD}/Metadata/PDG000000002.4290.metadata.tsv",
         "salmonella_cluster_list.tsv": f"{PD}/Clusters/PDG000000002.4290.reference_target.cluster_list.tsv"}
for f, u in FILES.items():
    if not os.path.exists(f) and not os.path.exists("sal_meta_raw.pkl" if "metadata" in f else f): print("downloading", f); urllib.request.urlretrieve(u, f)
os.makedirs("figs9", exist_ok=True); plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})

m = add_sector(load_metadata()); cat = load_catalog()
print(f"{len(m):,} Salmonella isolates in Pathogen Detection; {m.AST_phenotypes.notna().sum():,} with susceptibility results")
print("sector of every isolate:", m.sector.value_counts().to_dict())
print("source_type missing:", f"{m.source_type.isna().mean():.0%}", "| collection year missing:", f"{m.year.isna().mean():.1%}")

# ---------------------------------------------------------------- the analysis frame: USA, 2015-2024, four sectors
f = m[(m.country == "USA") & m.year.between(2015, 2024) & m.sector.isin(["Human", "Food", "Animal", "Environment"])].copy()
print(f"analysis frame (USA, 2015-2024, classified): {len(f):,} isolates; excluded unclassified in USA 2015-2024: {((m.country == 'USA') & m.year.between(2015, 2024) & (m.sector == 'Unclassified')).sum():,}")
f["state"] = f.geo_loc_name.str.extract(r"USA:\s*([A-Za-z]{2})\b")[0]
f["genes_list"] = f.AMR_genotypes_core.map(parse_genotype)                       # no core determinant = empty list, not missing
f["classes_list"] = f.genes_list.map(lambda g: [c for c in classes_of(g, cat) if c not in ("EFFLUX", "MULTIDRUG")])
f["n_classes"] = f.classes_list.map(len); f["mdr_geno"] = (f.n_classes >= 3).astype(int)
f["genes"] = f.genes_list.map(";".join); f["classes"] = f.classes_list.map(";".join)
f["has_ast"] = f.AST_phenotypes.notna().astype(int)
cols = ["target_acc", "asm_acc", "biosample_acc", "sector", "year", "state", "serovar_clean", "bioproject_center", "epi_type", "isolation_source", "genes", "classes", "n_classes", "mdr_geno", "has_ast", "AST_phenotypes", "minsame", "mindiff"]
f[cols].to_csv("sal_us_frame.csv.gz", index=False); print("wrote sal_us_frame.csv.gz")

tab = f.groupby("sector").agg(isolates=("target_acc", "size"), with_ast=("has_ast", "sum"), pct_ast=("has_ast", "mean"), median_genes=("genes_list", lambda s: s.map(len).median()), mdr_geno=("mdr_geno", "mean"),
                              top_centre=("bioproject_center", lambda s: s.value_counts().index[0]), top_serovar=("serovar_clean", lambda s: s.value_counts().index[0])).loc[["Human", "Food", "Animal", "Environment"]]
tab["pct_ast"] = (100 * tab.pct_ast).round(1); tab["mdr_geno"] = (100 * tab.mdr_geno).round(1); print(tab.to_string()); tab.to_csv("table_sampling_frame.csv")

# ---------------------------------------------------------------- Figure 1: who is sampled, when, and with what data
COL = {"Human": "#2b6cb0", "Food": "#dd6b20", "Animal": "#2f855a", "Environment": "#805ad5"}; order = ["Human", "Food", "Animal", "Environment"]
fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.1), gridspec_kw={"width_ratios": [1.5, 1, 1]})
y = f.groupby(["year", "sector"]).size().unstack(fill_value=0)[order]; bottom = np.zeros(len(y))
for s in order: ax[0].bar(y.index, y[s], bottom=bottom, color=COL[s], label=s, width=0.78); bottom += y[s].values
ax[0].set(title="Isolates sequenced per year, by sector", ylabel="Isolates", xticks=range(2015, 2025)); ax[0].set_xticklabels([str(i)[2:] for i in range(2015, 2025)]); ax[0].set_xlabel("Collection year (20xx)"); ax[0].legend(frameon=False, fontsize=8.5)
ax[1].bar(order, tab.isolates, color=[COL[s] for s in order]); ax[1].set(title="Isolates in the frame", ylabel="Isolates")
for i, v in enumerate(tab.isolates): ax[1].text(i, v, f"{v/1000:.0f}k" if v >= 1000 else str(v), ha="center", va="bottom", fontsize=8.5)
ax[2].bar(order, tab.pct_ast, color=[COL[s] for s in order]); ax[2].set(title="With a susceptibility test result", ylabel="% of isolates", ylim=(0, tab.pct_ast.max() * 1.28))
for i, (v, n) in enumerate(zip(tab.pct_ast, tab.with_ast)): ax[2].text(i, v, f"{v:.1f}%\n(n={n:,})", ha="center", va="bottom", fontsize=8)
for a in ax[1:]: a.tick_params(axis="x", rotation=25)
fig.tight_layout(); fig.savefig("figs9/fig1_sampling_frame.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 2: the serovar mix differs by sector (a confounder for every comparison)
top = f.serovar_clean.value_counts().head(9).index.tolist(); f["serovar_grp"] = np.where(f.serovar_clean.isin(top), f.serovar_clean, "Other")
mix = pd.crosstab(f.sector, f.serovar_grp, normalize="index").loc[order[::-1]]; mix = mix[top + ["Other"]]
fig, ax = plt.subplots(figsize=(9.6, 3.6)); left = np.zeros(len(mix)); cm = plt.get_cmap("tab20")
for i, s in enumerate(mix.columns):
    ax.barh(mix.index, mix[s], left=left, color="#cbd5e0" if s == "Other" else cm(i), label=s, edgecolor="white", linewidth=0.6); left += mix[s].values
ax.set(xlim=(0, 1), xlabel="Share of isolates in the sector", title="Serovar mix by sector"); ax.xaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(1.0))
ax.legend(frameon=False, fontsize=8.5, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5)); fig.tight_layout(); fig.savefig("figs9/fig2_serovar_mix.png", dpi=200); plt.close(fig)
mix.round(3).to_csv("table_serovar_mix.csv"); print(mix.round(2).to_string())
