"""
Chapter 9, Hands-on 6 - how many human / non-human pairs survive each SNP threshold?
The 89 clusters of Hands-on 4 (4 to 300 isolates in the frame, at least one human and one food or animal isolate, every isolate sharing determinants for three or more drug classes) are looked up in the
NCBI Pathogen Detection SNP trees (https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/latest_snps/SNP_trees/). Each tree has branch lengths in SNPs, so the distance between two isolates is the sum of the branch
lengths on the path between them (a parsimony path distance, which can exceed the true pairwise difference when a site changes more than once). For every threshold the script counts the human / non-human
isolate pairs of the frame that lie within it, and the clusters that hold at least one such pair.
Input: table_mixed_clusters_all.csv (ex9_4), sal_us_frame.csv.gz (ex9_1). Downloads the trees once into snp_trees/ (a few MB each) and keeps only the Newick file.
Output: table_snp_thresholds.csv, table_snp_min_distance.csv, table_snp_pairs.csv (every human / non-human pair and its distance), figs9/fig11_snp_thresholds.png.
"""
import os, re, io, tarfile, time, urllib.request
import numpy as np, pandas as pd
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs9", exist_ok=True); os.makedirs("snp_trees", exist_ok=True)
BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Salmonella/latest_snps/SNP_trees/"; THR = [0, 1, 2, 5, 10, 21, 50, 100]; MAXMB = 80
big = pd.read_csv("table_mixed_clusters_all.csv"); f = pd.read_csv("sal_us_frame.csv.gz", usecols=["target_acc", "sector"]); sec = dict(zip(f.target_acc, f.sector)); sec_nov = {k.split(".")[0]: v for k, v in sec.items()}
def listing():
    if os.path.exists("snp_trees/_listing.txt"): return open("snp_trees/_listing.txt").read().split()
    html = urllib.request.urlopen(BASE, timeout=120).read().decode(); names = sorted(set(re.findall(r'href="(PDS\d+\.\d+)\.tar\.gz"', html))); open("snp_trees/_listing.txt", "w").write("\n".join(names)); return names
names = listing(); latest = {n.split(".")[0]: n for n in names}; print(f"{len(names):,} trees listed")
def fetch(acc):
    out = f"snp_trees/{acc}.nwk"
    if os.path.exists(out): return out
    url = BASE + acc + ".tar.gz"
    for k in range(4):
        try:
            data = urllib.request.urlopen(url, timeout=120).read(); break
        except Exception: time.sleep(3 * (k + 1))
    else: return None
    if len(data) > MAXMB * 1e6: return None
    with tarfile.open(fileobj=io.BytesIO(data)) as t:
        m = [x for x in t.getmembers() if x.name.endswith(".newick")][0]; open(out, "wb").write(t.extractfile(m).read())
    return out
def dist_pairs(nwk):
    tree = Phylo.read(nwk, "newick"); par, blen, depth = {}, {}, {}
    def walk(c, p, d):
        par[id(c)] = p; depth[id(c)] = d
        for ch in c.clades: walk(ch, c, d + (ch.branch_length or 0.0))
    import sys; sys.setrecursionlimit(100000); walk(tree.root, None, 0.0); leaves = {}
    for t in tree.get_terminals():
        nm = (t.name or "").strip("'")
        s = sec.get(nm) or sec_nov.get(nm.split(".")[0])
        if s: leaves[nm] = (t, s)
    anc = {}
    for nm, (t, s) in leaves.items():
        path, c = [], t
        while c is not None: path.append(id(c)); c = par[id(c)]
        anc[nm] = path
    hum = [n for n, (t, s) in leaves.items() if s == "Human"]; non = [n for n, (t, s) in leaves.items() if s != "Human"]; out = []
    for h in hum:
        ah = {node: i for i, node in enumerate(anc[h])}
        for n in non:
            lca = next(node for node in anc[n] if node in ah); out.append((h, n, depth[id(leaves[h][0])] + depth[id(leaves[n][0])] - 2 * depth[lca], leaves[n][1]))
    return out, len(leaves), tree.count_terminals()
rows, mins, miss = [], [], []
for r in big.itertuples():
    base = r.cluster.split(".")[0]; acc = latest.get(base)
    if acc is None: miss.append(r.cluster); continue
    nwk = fetch(acc)
    if nwk is None: miss.append(r.cluster); continue
    pr, nleaf, ntip = dist_pairs(nwk)
    for h, n, d, s in pr: rows.append(dict(cluster=base, tree=acc, human=h, other=n, sector=s, snp_distance=d))
    if pr: mins.append(dict(cluster=base, tree=acc, serovar=r.serovar, frame_isolates=r.n, tree_tips=ntip, frame_leaves=nleaf, pairs=len(pr), min_distance=min(d for _, _, d, _ in pr), min_human_food=min([d for _, _, d, s in pr if s == "Food"] or [np.nan]), min_human_animal=min([d for _, _, d, s in pr if s == "Animal"] or [np.nan])))
P = pd.DataFrame(rows); P.to_csv("table_snp_pairs.csv", index=False); M = pd.DataFrame(mins); print(f"{len(M)} of {len(big)} clusters analysed; {len(miss)} not found or too large: {miss[:6]}; {len(P):,} human / non-human pairs")
T = []
for t in THR:
    q = P[P.snp_distance <= t]; T.append(dict(threshold_snps=t, pairs=len(q), pairs_human_food=int((q.sector == "Food").sum()), pairs_human_animal=int((q.sector == "Animal").sum()), pairs_human_environment=int((q.sector == "Environment").sum()),
                                           clusters_with_a_pair=int(q.cluster.nunique()), clusters_analysed=len(M), share_of_clusters=round(q.cluster.nunique() / max(len(M), 1), 4), share_of_pairs=round(len(q) / max(len(P), 1), 4)))
T = pd.DataFrame(T); T.to_csv("table_snp_thresholds.csv", index=False); M.to_csv("table_snp_min_distance.csv", index=False); print(T.to_string(index=False)); print(M.min_distance.describe().round(1).to_dict())
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3)); xs = np.arange(len(THR)); w = 0.38
ax[0].bar(xs - w / 2, [int(P[(P.snp_distance <= t) & (P.sector == "Food")].cluster.nunique()) for t in THR], w, color="#dd6b20", label="human and food"); ax[0].bar(xs + w / 2, [int(P[(P.snp_distance <= t) & (P.sector == "Animal")].cluster.nunique()) for t in THR], w, color="#2f855a", label="human and animal")
ax[0].set_xticks(xs); ax[0].set_xticklabels([f"≤{t}" for t in THR]); ax[0].set(xlabel="SNP threshold", ylabel=f"Clusters with a pair within the threshold (of {len(M)})", title="Cross-sector pairs depend on the threshold"); ax[0].legend(frameon=False)
v = M.min_distance.clip(lower=0.5); ax[1].hist(np.log10(v), bins=np.linspace(-0.4, np.log10(max(v.max(), 200)) + 0.1, 24), color="#2b6cb0")
for t, lab in [(5, "5"), (21, "21"), (100, "100")]: ax[1].axvline(np.log10(t), color="#c53030", ls=":"); ax[1].text(np.log10(t) + 0.03, ax[1].get_ylim()[1] * 0.92, lab, color="#c53030", fontsize=9)
ax[1].set_xlim(-0.5, np.log10(120)); ax[1].text(0.02, 0.97, "0 SNPs plotted at 0.5", transform=ax[1].transAxes, fontsize=8.5, color="#718096", va="top"); ax[1].set(xlabel="Closest human / non-human pair in the cluster (SNPs, log scale)", ylabel="Clusters", title="Distance of the closest cross-sector pair"); ax[1].set_xticks(np.log10([1, 3, 10, 30, 100])); ax[1].set_xticklabels(["1", "3", "10", "30", "100"])
fig.tight_layout(); fig.savefig("figs9/fig11_snp_thresholds.png", dpi=200); print("figure written")
