"""
Chapter 4, Hands-on 4 - Chaining, fragmentation and threshold sensitivity in real clusters (Listeria monocytogenes, NCBI Pathogen Detection).
(a) How wide are clusters defined by single linkage at 50 SNPs? (diameter = largest pairwise distance inside a cluster)
(b) How do the same 74,758 isolates divide if single linkage is applied at tighter thresholds? The SNP-distance pairs are streamed and connected components are computed with a sparse graph; at 50 SNPs this must reproduce
    the pipeline's own clusters, which checks the method.
(c) How often does an isolate have a neighbour of the same or of the opposite isolation type (clinical against environmental) within a threshold?
Inputs: listeria_SNP_distances.parquet, listeria_cluster_list.tsv, listeria_all_isolates.tsv, listeria_cluster_stats.csv. Output: table_cluster_diameter.csv, table_fragmentation.csv, table_neighbours.csv, Figure 4.
"""
import numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
import os, subprocess
BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Listeria/latest_snps/Clusters/PDG000000001.4985.reference_target"
def ensure(local, suffix):
    """Download an NCBI file if it is not present. NCBI re-runs the analysis and retires old versions, so later downloads may differ from the version used in the chapter (4985, 13 November 2026)."""
    if not os.path.exists(local): subprocess.run(["curl", "-s", "-L", "--retry", "4", "-o", local, f"{BASE}.{suffix}"], check=True)
ensure("listeria_SNP_distances.parquet", "SNP_distances.parquet"); ensure("listeria_cluster_list.tsv", "cluster_list.tsv"); ensure("listeria_all_isolates.tsv", "all_isolates.tsv")
S = pd.read_csv("listeria_cluster_stats.csv", index_col=0); big = S[S.isolates >= 10]
D = pd.DataFrame([dict(clusters_10plus=len(big), diameter_median=big.diameter.median(), diameter_p90=big.diameter.quantile(0.9), diameter_max=big.diameter.max(), share_over_50=(big.diameter > 50).mean() * 100, share_over_100=(big.diameter > 100).mean() * 100, pairs_le5_pct=100 * big.le5.sum() / big.pairs.sum(), pairs_le10_pct=100 * big.le10.sum() / big.pairs.sum(), pairs_le20_pct=100 * big.le20.sum() / big.pairs.sum(), pairs_le50_pct=100 * big.le50.sum() / big.pairs.sum(), spearman_size_diameter=big.isolates.corr(big.diameter, method="spearman"))]); D.round(3).to_csv("table_cluster_diameter.csv", index=False); print(D.round(2).T.to_string())
cl = pd.read_csv("listeria_cluster_list.tsv", sep="\t"); ids = pd.Index(cl.target_acc); N = len(ids); pds = cl.set_index("target_acc").PDS_acc; ei, ej, ed = [], [], []
for b in pq.ParquetFile("listeria_SNP_distances.parquet").iter_batches(batch_size=1_000_000, columns=["target_acc_1", "target_acc_2", "delta_positions_unambiguous"]):
    d = b.to_pandas(); k = d.delta_positions_unambiguous.values <= 50; d = d[k]; ei.append(ids.get_indexer(d.target_acc_1)); ej.append(ids.get_indexer(d.target_acc_2)); ed.append(d.delta_positions_unambiguous.values.astype(np.int16))
ei, ej, ed = np.concatenate(ei), np.concatenate(ej), np.concatenate(ed); ok = (ei >= 0) & (ej >= 0); ei, ej, ed = ei[ok], ej[ok], ed[ok]; print(f"{len(ed):,} pairs within 50 SNPs among {N:,} isolates")
rows = []; comp50 = None
for T in (0, 2, 5, 10, 20, 50):
    m = ed <= T; g = coo_matrix((np.ones(m.sum()), (ei[m], ej[m])), shape=(N, N)); nc, lab = connected_components(g, directed=False); sizes = np.bincount(lab)
    if T == 50: comp50 = lab
    sub = pd.Series(lab, index=ids).groupby(pds.reindex(ids).values).nunique(); rows.append(dict(threshold=T, groups=nc, singletons_pct=100 * (sizes == 1).sum() / N, isolates_in_groups_of_2plus_pct=100 * sizes[sizes > 1].sum() / N, largest_group=int(sizes.max()), median_group_size_2plus=float(np.median(sizes[sizes > 1])) if (sizes > 1).any() else np.nan, mean_subgroups_per_50snp_cluster=sub.mean()))
F = pd.DataFrame(rows); F.round(2).to_csv("table_fragmentation.csv", index=False); print(F.round(2).to_string(index=False))
agree = (pd.Series(comp50, index=ids).groupby(pds.reindex(ids).values).nunique() == 1).mean(); print(f"components at 50 SNPs reproduce the pipeline clusters: each pipeline cluster is one component in {agree:.1%} of {pds.nunique():,} clusters; components {F[F.threshold == 50].groups.iloc[0]:,} against {pds.nunique():,} clusters")
A = pd.read_csv("listeria_all_isolates.tsv", sep="\t", na_values=["NULL"]); NB = []
for T in (0, 2, 5, 10, 20, 50):
    NB.append(dict(threshold=T, same_type_neighbour_pct=100 * (A.min_dist_same <= T).sum() / A.min_dist_same.notna().sum(), opposite_type_neighbour_pct=100 * (A.min_dist_opp <= T).sum() / A.min_dist_opp.notna().sum(), same_n=int(A.min_dist_same.notna().sum()), opp_n=int(A.min_dist_opp.notna().sum()), opp_neighbours_per_100_isolates=100 * (A.min_dist_opp <= T).sum() / len(A)))
NBd = pd.DataFrame(NB); NBd.round(2).to_csv("table_neighbours.csv", index=False); print(NBd.round(2).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6)); x = np.sort(big.diameter.values); ax[0].plot(x, np.arange(1, len(x) + 1) / len(x) * 100, color="#2b6cb0", lw=2); ax[0].axvline(50, color="#c53030", ls="--", lw=1); ax[0].text(53, 8, "50-SNP threshold", color="#c53030", fontsize=8.5)
ax[0].set(xlabel="Cluster diameter (largest pairwise SNP distance)", ylabel="Clusters of 10 or more isolates (cumulative %)", title=f"Single linkage at 50 SNPs: {len(big):,} clusters"); ax[0].set_xlim(0, 130)
ax[1].bar(range(len(F)), F.groups, color="#2b6cb0"); ax[1].set_xticks(range(len(F))); ax[1].set_xticklabels([str(t) for t in F.threshold]); ax[1].set(xlabel="Single-linkage threshold (SNPs)", ylabel="Groups among 74,758 isolates", title="One cluster becomes many"); a2 = ax[1].twinx(); a2.plot(range(len(F)), F.singletons_pct, "o-", color="#c53030"); a2.set_ylabel("Isolates left alone (%)", color="#c53030"); a2.spines["top"].set_visible(False)
ax[2].plot(NBd.threshold, NBd.same_type_neighbour_pct, "o-", color="#2b6cb0", label="Nearest same-type neighbour"); ax[2].plot(NBd.threshold, NBd.opposite_type_neighbour_pct, "o-", color="#dd6b20", label="Nearest opposite-type neighbour"); ax[2].set(xlabel="Threshold (SNPs)", ylabel="Isolates with a neighbour within it (%)", title="Links between isolation types"); ax[2].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("figs4/fig4_real_clusters.png", dpi=200); plt.close(fig); print("figure written")
