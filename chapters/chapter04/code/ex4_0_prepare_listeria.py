"""
Chapter 4 data preparation: streams the NCBI Pathogen Detection Listeria SNP-distance table (17.3 million within-cluster pairs, 2.9 GB as text, 148 MB as parquet) and writes
  listeria_cluster_stats.csv  (one row per 50-SNP single-linkage cluster: isolates, pairs, diameter = largest pairwise SNP distance, pairs within 5, 10, 20 and 50 SNPs) and
  listeria_pairs_selected.csv.gz  (all pairs of two clusters chosen for the worked examples: the largest-diameter cluster with 60-120 isolates and the largest-diameter cluster with 150-400 isolates).
Input files (downloaded from https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Listeria/latest_snps/Clusters/): PDG000000001.4985.reference_target.SNP_distances.parquet and .cluster_list.tsv.
"""
import numpy as np, pandas as pd, pyarrow.parquet as pq
PQ, CL = "listeria_SNP_distances.parquet", "listeria_cluster_list.tsv"; COLS = ["PDS_acc", "target_acc_1", "target_acc_2", "delta_positions_unambiguous"]
import os, subprocess
BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Listeria/latest_snps/Clusters/PDG000000001.4985.reference_target"
def ensure(local, suffix):
    """Download an NCBI file if it is not present. NCBI re-runs the analysis and retires old versions, so later downloads may differ from the version used in the chapter (4985, 13 November 2026)."""
    if not os.path.exists(local): subprocess.run(["curl", "-s", "-L", "--retry", "4", "-o", local, f"{BASE}.{suffix}"], check=True)
ensure("listeria_SNP_distances.parquet", "SNP_distances.parquet"); ensure("listeria_cluster_list.tsv", "cluster_list.tsv")
f = pq.ParquetFile(PQ); parts = []
for b in f.iter_batches(batch_size=1_000_000, columns=COLS):
    d = b.to_pandas(); d["x"] = d.delta_positions_unambiguous.astype(float)
    g = d.groupby("PDS_acc").x.agg(pairs="size", diameter="max", total="sum"); g["le0"] = d.assign(k=d.x <= 0).groupby("PDS_acc").k.sum(); g["le5"] = d.assign(k=d.x <= 5).groupby("PDS_acc").k.sum(); g["le10"] = d.assign(k=d.x <= 10).groupby("PDS_acc").k.sum(); g["le20"] = d.assign(k=d.x <= 20).groupby("PDS_acc").k.sum(); g["le50"] = d.assign(k=d.x <= 50).groupby("PDS_acc").k.sum(); parts.append(g)
S = pd.concat(parts).groupby(level=0).agg(pairs=("pairs", "sum"), diameter=("diameter", "max"), total=("total", "sum"), le0=("le0", "sum"), le5=("le5", "sum"), le10=("le10", "sum"), le20=("le20", "sum"), le50=("le50", "sum"))
cl = pd.read_csv(CL, sep="\t"); n = cl.groupby("PDS_acc").size().rename("isolates"); S = S.join(n); S["mean_distance"] = S.total / S.pairs; S = S.drop(columns="total"); S.to_csv("listeria_cluster_stats.csv")
print(f"{len(S):,} clusters with pair data; {int(S.pairs.sum()):,} pairs; isolates per cluster median {S.isolates.median():.0f}, max {int(S.isolates.max())}")
big = S[S.isolates >= 10]; print(f"clusters of at least 10 isolates: {len(big):,}; diameter above 50 SNPs in {(big.diameter > 50).mean():.1%}, above 100 in {(big.diameter > 100).mean():.1%}, above 200 in {(big.diameter > 200).mean():.1%}; median diameter {big.diameter.median():.0f}; median mean pairwise distance {big.mean_distance.median():.1f}")
pick = []
for lo, hi in [(60, 120), (150, 400)]:
    c = S[(S.isolates >= lo) & (S.isolates <= hi)].sort_values("diameter", ascending=False).head(1); pick += list(c.index); print("selected:", c[["isolates", "pairs", "diameter", "mean_distance", "le10", "le50"]].round(1).to_dict("records"), c.index.tolist())
keep = []
for b in f.iter_batches(batch_size=1_000_000, columns=COLS):
    d = b.to_pandas(); k = d[d.PDS_acc.isin(pick)]
    if len(k): keep.append(k)
P = pd.concat(keep); P.to_csv("listeria_pairs_selected.csv.gz", index=False); print("pairs saved:", len(P), P.PDS_acc.value_counts().to_dict())
