"""
Chapter 4, Hands-on 1 - Distance matrices, thresholds and linkage rules.
Part A repairs and reruns the uploaded Chapter_4_3.py (ten simulated isolates, 50 SNP loci). As uploaded it stops with a TypeError because the genotype array is floating point and the
script applies the bitwise XOR operator to it; here the genotypes are integers. Part B applies the same steps to real data: all pairwise SNP distances of two clusters of Listeria monocytogenes
from the NCBI Pathogen Detection pipeline (111 and 220 isolates), clustered with single, average and complete linkage at thresholds of 5, 10, 20 and 50 SNPs.
Output: table_linkage_rules.csv, Figure 1.
"""
import numpy as np, pandas as pd
import seaborn as sns
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
# ---------------------------------------------------------------- Part A: the uploaded simulation, repaired
np.random.seed(88); n, L = 10, 50; labels = [f"Pat_{i:02d}" for i in range(1, n + 1)]; G = np.zeros((n, L), dtype=int)
G[0] = np.random.binomial(1, 0.02, L)
for i in range(1, 4): G[i] = np.bitwise_xor(G[0], np.random.binomial(1, 0.04, L))
G[4] = np.random.binomial(1, 0.35, L)
for i in range(5, 7): G[i] = np.bitwise_xor(G[4], np.random.binomial(1, 0.04, L))
G[7] = np.random.binomial(1, 0.40, L); G[8] = np.random.binomial(1, 0.45, L); G[9] = np.random.binomial(1, 0.50, L)
D = np.abs(G[:, None, :] - G[None, :, :]).sum(axis=2); dfm = pd.DataFrame(D, index=labels, columns=labels); pairs = [(labels[i], labels[j], int(D[i, j])) for i in range(n) for j in range(i + 1, n) if D[i, j] <= 12]
print("Part A: pairs within 12 SNPs:", pairs); LA = linkage(squareform(D), "average"); print("average-linkage clusters at 12 SNPs:", fcluster(LA, 12, "distance").tolist())
# ---------------------------------------------------------------- Part B: real clusters
P = pd.read_csv("listeria_pairs_selected.csv.gz"); rows = []; mats = {}
for pds, g in P.groupby("PDS_acc"):
    ids = sorted(set(g.target_acc_1) | set(g.target_acc_2)); ix = {k: i for i, k in enumerate(ids)}; m = np.zeros((len(ids), len(ids)))
    for a, b, x in zip(g.target_acc_1, g.target_acc_2, g.delta_positions_unambiguous): m[ix[a], ix[b]] = m[ix[b], ix[a]] = x
    mats[pds] = m; iu = np.triu_indices(len(ids), 1); print(f"{pds}: {len(ids)} isolates, {len(iu[0]):,} pairs, median distance {np.median(m[iu]):.0f}, diameter {m.max():.0f}, pairs within 10 SNPs {(m[iu] <= 10).mean():.1%}, within 50 SNPs {(m[iu] <= 50).mean():.1%}")
    for method in ("single", "average", "complete"):
        Z = linkage(squareform(m, checks=False), method)
        for T in (5, 10, 20, 50):
            lab = fcluster(Z, T, "distance"); sizes = np.bincount(lab)[1:]; big = np.argmax(sizes) + 1; sub = m[np.ix_(lab == big, lab == big)]
            rows.append(dict(cluster=pds, isolates=len(ids), linkage=method, threshold=T, clusters=len(sizes), largest=int(sizes.max()), singletons=int((sizes == 1).sum()), diameter_of_largest=float(sub.max())))
R = pd.DataFrame(rows); R.to_csv("table_linkage_rules.csv", index=False); print(R[R.cluster == "PDS000041936.59"].pivot_table(index="threshold", columns="linkage", values=["clusters", "largest", "diameter_of_largest"]).round(0).to_string())
pds = "PDS000041936.59"; m = mats[pds]; Z = linkage(squareform(m, checks=False), "average"); order = dendrogram(Z, no_plot=True)["leaves"]
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [0.95, 1, 1]})
sns.heatmap(dfm, annot=True, fmt="d", cmap="YlOrRd", cbar=False, ax=ax[0], annot_kws={"size": 6}, linewidths=0.3); ax[0].set_title("Repaired simulation: SNP distances (10 isolates)", fontsize=10); ax[0].tick_params(labelsize=7)
im = ax[1].imshow(m[np.ix_(order, order)], cmap="YlOrRd_r", vmin=0, vmax=60); ax[1].set(title=f"Real Listeria cluster, {len(m)} isolates", xticks=[], yticks=[]); fig.colorbar(im, ax=ax[1], pad=0.01).set_label("SNP distance (capped at 60)", fontsize=8.5)
for method, c in [("single", "#c53030"), ("average", "#2b6cb0"), ("complete", "#2f855a")]:
    s = R[(R.cluster == pds) & (R.linkage == method)]; ax[2].plot(s.threshold, s.clusters, "o-", color=c, label=method)
ax[2].set_xscale("log"); ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[2].set_xticks([5, 10, 20, 50]); ax[2].set_xticklabels(["5", "10", "20", "50"]); ax[2].set(xlabel="Threshold (SNPs)", ylabel="Clusters in the 220 isolates", title="The same threshold, three rules"); ax[2].legend(frameon=False, fontsize=8.5, title="Linkage", title_fontsize=8.5)
fig.tight_layout(); fig.savefig("figs4/fig1_distance_matrices.png", dpi=200); plt.close(fig); print("figure written")
