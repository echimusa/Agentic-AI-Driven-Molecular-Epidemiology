"""
Chapter 9 - the figures of the extension, drawn from the result tables of ex9_3b, ex9_5b, ex9_7 and from the stages of the Quadripartite guidance.
Figure 10 two tools on eight genomes (table_two_tools_summary.csv, table_two_tools.csv, table_two_tools_loose.csv, table_two_tools_background.csv; hits per genome are recomputed from rgi_outputs/)
Figure 7  genotype and phenotype disagree (table_discordance.csv, table_azithromycin_period.csv)
Figure 3  completeness of the public metadata against the core fields (table_metadata_completeness.csv)
Diagram   the stages of integrated surveillance and where a genomic study fits (drawn, no data)
(Figure 11, the SNP thresholds, is drawn by ex9_6_snp_thresholds.py.)
"""
import glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs9", exist_ok=True)
COL = {"Human": "#2b6cb0", "Food": "#dd6b20", "Animal": "#2f855a", "Environment": "#805ad5"}; order = ["Human", "Food", "Animal", "Environment"]; INK = "#1a202c"; GREY = "#718096"
# ---------------------------------------------------------------- Figure 9
S = pd.read_csv("table_two_tools_summary.csv").iloc[0]; C = pd.read_csv("table_two_tools.csv"); Lz = pd.read_csv("table_two_tools_loose.csv")
accs = sorted(os.path.basename(p)[:15] for p in glob.glob("amrfinder_outputs/*.amrfinder.tsv"))
def n_amr(a):
    d = pd.read_csv(f"amrfinder_outputs/{a}.amrfinder.tsv", sep="\t"); d = d[(d.Type == "AMR") & (d.Scope == "core") & ~d.Method.str.startswith("PARTIAL") & ~d.Method.str.startswith("INTERNAL_STOP")]; return d["Element symbol"].nunique()
amr = [n_amr(a) for a in accs]; rgi = [len(pd.read_csv(f"rgi_outputs/{a}.rgi.txt", sep="\t")) for a in accs]
fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.5), gridspec_kw={"width_ratios": [1.35, 1]}); x = np.arange(8); w = 0.38
ax[0].bar(x - w / 2, amr, w, color="#1B998B", label="AMRFinderPlus (core elements)"); ax[0].bar(x + w / 2, rgi, w, color="#E4572E", label="RGI (Perfect and Strict hits)")
for i, (u, v) in enumerate(zip(amr, rgi)): ax[0].text(i - w / 2, u + 0.6, str(u), ha="center", fontsize=8.5); ax[0].text(i + w / 2, v + 0.6, str(v), ha="center", fontsize=8.5)
ax[0].set_xticks(x); ax[0].set_xticklabels([a[4:13] for a in accs], rotation=35, ha="right", fontsize=8.5); ax[0].set(ylabel="Elements or hits reported", title="Same eight genomes, two answers", ylim=(0, 48)); ax[0].legend(frameon=False, fontsize=9, loc="upper left")
cats = ["same name", "same family, different allele", "AMRFinderPlus only", "RGI only"]; cnt = [int((C.result == c).sum()) for c in cats]; cc = ["#2f855a", "#a0aec0", "#1B998B", "#E4572E"]
ax[1].barh(range(4), cnt, color=cc); ax[1].set_yticks(range(4)); ax[1].set_yticklabels(["Same name", "Same family,\ndifferent allele", "AMRFinderPlus only", "RGI only"]); ax[1].invert_yaxis()
for i, v in enumerate(cnt): ax[1].text(v + 0.3, i, str(v), va="center")
ax[1].set(xlabel="Elements in the four genomes with a call", title=f"After setting aside {int(S.rgi_background)} background hits"); ax[1].set_xlim(0, 21)
ax[1].text(0.98, 0.04, f"all {len(Lz)} AMRFinderPlus-only elements\nappear as Loose RGI hits", transform=ax[1].transAxes, ha="right", va="bottom", fontsize=9, color=GREY)
fig.tight_layout(); fig.savefig("figs9/fig10_two_tools.png", dpi=200); plt.close(fig)
# ---------------------------------------------------------------- Figure 10
D = pd.read_csv("table_discordance.csv"); D = D[D.fn > 0].sort_values("fn"); P = pd.read_csv("table_azithromycin_period.csv")
fig, ax = plt.subplots(1, 2, figsize=(11.8, 5.0), gridspec_kw={"width_ratios": [1.5, 1]}); y = np.arange(len(D))
ax[0].barh(y, D.fn_other_determinants_only, color="#3D7EA6", label="carries determinants of other classes only"); ax[0].barh(y, D.fn_no_core_determinant, left=D.fn_other_determinants_only, color="#F3A712", label="carries no core determinant")
ax[0].set_yticks(y); ax[0].set_yticklabels(D.drug, fontsize=9); ax[0].set(xlabel="Resistant isolates called susceptible by the rule", title=f"{int(D.fn.sum())} false-susceptible calls, by what the isolate carries"); ax[0].legend(frameon=False, fontsize=9, loc="lower right")
ax[1].bar(range(3), 100 * P.sensitivity, color="#2E294E"); ax[1].set_xticks(range(3)); ax[1].set_xticklabels(P.period)
for i, r in enumerate(P.itertuples()): ax[1].text(i, 100 * r.sensitivity + 2, f"{int(r.tp)} of {int(r.n_R)}", ha="center", fontsize=9)
ax[1].set(ylabel="Resistant isolates with a macrolide determinant (%)", title="Azithromycin sensitivity by period", ylim=(0, 85)); fig.tight_layout(); fig.savefig("figs9/fig7_discordance.png", dpi=200); plt.close(fig)
# ---------------------------------------------------------------- Figure 11
M = pd.read_csv("table_metadata_completeness.csv"); V = M[order].values; fig, ax = plt.subplots(figsize=(8.6, 4.4)); im = ax.imshow(V, cmap="Blues", vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(order); ax.set_yticks(range(len(M))); ax.set_yticklabels(M.field, fontsize=9); ax.xaxis.tick_top(); [ax.spines[s].set_visible(False) for s in ax.spines]
for i in range(V.shape[0]):
    for j in range(V.shape[1]): ax.text(j, i, f"{V[i, j]:.1f}", ha="center", va="center", color="white" if V[i, j] > 55 else INK, fontsize=9)
ax.set_title("Isolates carrying each core field (%), United States 2015-2024", fontsize=10.5, pad=34); fig.tight_layout(); fig.savefig("figs9/fig3_metadata_completeness.png", dpi=200); plt.close(fig)
# ---------------------------------------------------------------- diagram: stages of integrated surveillance
fig, ax = plt.subplots(figsize=(12.4, 4.9)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.9)
ax.text(0.15, 4.58, "Integrated surveillance runs in five stages; a genomic study touches every one", fontsize=13.5, fontweight="bold", color=INK, va="center")
stages = [("1  Situation analysis", "Who samples already?\nWhich objective?\nWhich resources?"), ("2  Sampling and data", "Sectors, interfaces\nIndicator organism\nCore data fields"), ("3  Laboratory and QC", "Susceptibility tests\nOptional: genomes,\nmetagenomes"),
          ("4  Analysis, reporting", "Integrate sectors\nSNP thresholds, links\nShare, communicate"), ("5  Evaluation", "Did the data change\na decision?\nWhat did it cost?")]
for i, (t, b) in enumerate(stages):
    x0 = 0.15 + i * 2.43; ax.add_patch(FancyBboxPatch((x0, 2.0), 2.2, 2.05, boxstyle="round,pad=0.02,rounding_size=0.12", fc="white", ec="#c8c5bb", lw=1.4)); ax.text(x0 + 0.14, 3.8, t, fontsize=9.6, fontweight="bold", color=INK, va="center"); ax.text(x0 + 0.14, 2.95, b, fontsize=9.0, color="#4a5568", va="center", linespacing=1.55)
    if i < 4: ax.annotate("", xy=(x0 + 2.4, 3.0), xytext=(x0 + 2.22, 3.0), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.2))
objs = [("1  Patterns and trends", 0.15), ("2  Risk analysis", 3.25), ("3  Alert on emerging resistance", 6.35), ("4  Inform and judge interventions", 9.45)]
for t, x0 in objs: ax.add_patch(FancyBboxPatch((x0, 0.55), 2.85, 0.62, boxstyle="round,pad=0.02,rounding_size=0.1", fc="#E4F0F7", ec="#3D7EA6", lw=1.2)); ax.text(x0 + 1.43, 0.86, t, ha="center", va="center", fontsize=9.6, color=INK)
ax.text(0.15, 1.55, "Four objectives the stages serve", fontsize=9.5, color=GREY, va="center"); ax.text(0.15, 0.2, "Genomic additions in the core data: resistance gene profile, mobilome, virulome, phylogeny with working definitions, BioProject.", fontsize=9.3, color=GREY, va="center")
fig.savefig("figs9/diagram_integrated_surveillance.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("figures written")
