"""
Chapter 5, Hands-on 1 - Distances and distance trees on real alignments.
Data: 195 Ebola virus genomes (18,961 positions, Sierra Leone 2014; Chapter 13) and 60 influenza A (H3N2) hemagglutinin sequences from NCBI (12 per year, 2014-2018; 1,762 aligned positions).
Steps: (1) p-distance, Jukes-Cantor (JC69) and Kimura (K80) distances written from the formulas of the chapter; (2) a check of the uploaded Chapter_5_2.py calls: Biopython's 'identity' model is the p-distance and
'jc69' is not a model; (3) neighbour-joining (NJ) and UPGMA trees, compared with the maximum-likelihood (ML) tree of IQ-TREE by the normalised Robinson-Foulds distance.
Output: table_distances.csv, table_tree_comparison.csv, Figure 1. Needs h3n2.treefile and ebola.treefile (written by run_iqtree.sh).
"""
import numpy as np, pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo_utils import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
sets = {"H3N2 HA (60 sequences, 1,762 sites)": ("h3n2_ha_aln.fasta", "h3n2.treefile"), "Ebola virus (195 genomes, 18,961 sites)": ("ebola_aln.fasta", "ebola.treefile")}
rows, comp, store = [], [], {}
for name, (aln, tf) in sets.items():
    names, X = read_aln(aln); E = encode(X); p, P, Q = pairwise(E); iu = np.triu_indices(len(names), 1); dJ = jc(p); dK = k80(P, Q)
    rows.append(dict(dataset=name, taxa=len(names), sites=X.shape[1], mean_p=p[iu].mean(), max_p=p[iu].max(), max_jc=dJ[iu].max(), jc_over_p_max=(dJ[iu] / np.where(p[iu] > 0, p[iu], np.nan)).max(), k80_minus_jc_max=(dK[iu] - dJ[iu]).max(), ts_share=(P[iu].sum() / (P[iu].sum() + Q[iu].sum())), identical_pairs_pct=100 * (p[iu] == 0).mean()))
    ml = Phylo.read(tf, "newick"); L = X.shape[1]; sp_ml = splits(ml, names, min_len=0.5 / L); res = {}
    for label, D in [("NJ (JC69)", dJ), ("NJ (p-distance)", p)]:
        t = DistanceTreeConstructor().nj(to_dm(names, np.nan_to_num(D))); res[label] = t
    res["UPGMA (JC69)"] = DistanceTreeConstructor().upgma(to_dm(names, np.nan_to_num(dJ)))
    for label, t in res.items(): s = splits(t, names, min_len=0.5 / L); comp.append(dict(dataset=name, method=label, splits=len(s), ml_splits=len(sp_ml), shared=len(s & sp_ml), rf=rf(s, sp_ml)))
    store[name] = (names, p, dJ, dK, res, ml)
    if name.startswith("H3N2"):
        bc = DistanceCalculator("identity").get_distance(__import__("Bio").AlignIO.read(aln, "fasta")); B = np.array(bc.matrix + [[]] * 0, dtype=object); M = np.zeros((len(names), len(names)))
        for i in range(len(names)):
            for j in range(i): M[i, j] = M[j, i] = bc.matrix[i][j]
        print(f"Biopython 'identity' distance against the p-distance written here: largest absolute difference {np.abs(M - p).max():.4f}; against JC69: {np.abs(M - dJ).max():.4f}")
try: DistanceCalculator("jc69")
except ValueError as e: print("DistanceCalculator('jc69') raises ValueError:", str(e)[:60])
D = pd.DataFrame(rows); D.round(5).to_csv("table_distances.csv", index=False); print(D.round(4).T.to_string()); C = pd.DataFrame(comp); C.round(3).to_csv("table_tree_comparison.csv", index=False); print(C.round(3).to_string(index=False))
# ------------------------------------------------ figure
h3 = list(sets)[0]; names, p, dJ, dK, res, ml = store[h3]; yrs = pd.read_csv("h3n2_meta.csv").set_index("accession").year; yr = {n: int(yrs[n]) for n in names}; pal = {2014: "#2b6cb0", 2015: "#2f855a", 2016: "#d69e2e", 2017: "#dd6b20", 2018: "#c53030"}
fig, ax = plt.subplots(1, 3, figsize=(14.0, 4.8), gridspec_kw={"width_ratios": [1, 1.05, 0.9]})
x = np.linspace(0, 0.7, 300); ax[0].plot(x, x, color="black", ls=":", lw=1, label="p (no correction)"); ax[0].plot(x, jc(x), color="#2b6cb0", lw=2, label="JC69")
for r, c in [(0.25, "#2f855a"), (0.7, "#dd6b20")]: Pp = x * r; Qq = x - Pp; ok = (1 - 2 * Pp - Qq > 0) & (1 - 2 * Qq > 0); ax[0].plot(x[ok], k80(Pp[ok], Qq[ok]), color=c, lw=1.6, ls="--", label=f"K80, {int(r * 100)}% of differences are transitions")
iu = np.triu_indices(len(names), 1); ax[0].scatter(p[iu], dJ[iu], s=3, color="#c53030", alpha=0.3, label="H3N2 pairs"); ax[0].set(xlim=(0, 0.7), ylim=(0, 2.0), xlabel="Observed proportion of differences, p", ylabel="Estimated substitutions per site, d", title="Correcting for multiple hits"); ax[0].legend(frameon=False, fontsize=7.4, loc="upper left")
tree = res["NJ (JC69)"]; draw_tree(ax[1], tree, colors={n: pal[yr[n]] for n in names}); ax[1].set(title="H3N2 HA: neighbour-joining tree (JC69)", xlabel="Substitutions per site", yticks=[]); ax[1].spines["left"].set_visible(False)
for y, c in pal.items(): ax[1].plot([], [], "o", color=c, label=str(y))
ax[1].legend(frameon=False, fontsize=7.6, loc="lower right", title="Year", title_fontsize=7.6)
cc = C.copy(); cc["short"] = cc.dataset.str.split(" ").str[0]; w = 0.26
for k, m in enumerate(["NJ (JC69)", "NJ (p-distance)", "UPGMA (JC69)"]): s = cc[cc.method == m]; ax[2].bar(np.arange(len(s)) + (k - 1) * w, s.rf, w, label=m, color=["#2b6cb0", "#63b3ed", "#dd6b20"][k])
ax[2].set_xticks(range(2)); ax[2].set_xticklabels(cc[cc.method == "NJ (JC69)"].short); ax[2].set(ylabel="Robinson-Foulds distance to the ML tree", title="Distance trees against the ML tree", ylim=(0, 1)); ax[2].legend(frameon=False, fontsize=7.6)
fig.tight_layout(); fig.savefig("figs5/fig1_distances_nj.png", dpi=200); plt.close(fig); print("figure written")
