"""
Chapter 9, Hands-on 4 - (a) Do SNP clusters mix human and non-human isolates more than the sampling frame predicts?  (b) Which resistance genes travel together?
(a) Observed counts of mixed clusters are compared with a permutation null: sector labels are shuffled among isolates of the same serovar and collection year,
    keeping the cluster structure fixed. (b) Gene co-occurrence (phi coefficient) and Louvain communities.
Input: sal_us_frame.csv.gz (ex9_1) and salmonella_cluster_list.tsv. Output: Figures 8 and 9, table_mixed_clusters.csv, table_mixed_clusters_all.csv (all 89 clusters, for ex9_6), table_modules.csv.
"""
import numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
order = ["Human", "Food", "Animal", "Environment"]; COL = {"Human": "#2b6cb0", "Food": "#dd6b20", "Animal": "#2f855a", "Environment": "#805ad5"}
rng = np.random.default_rng(9); NPERM = 1000

f = pd.read_csv("sal_us_frame.csv.gz", dtype={"state": str}); f["genes"] = f.genes.fillna(""); f["classes"] = f.classes.fillna("")
cl = pd.read_csv("salmonella_cluster_list.tsv", sep="\t", usecols=["PDS_acc", "target_acc"]); f = f.merge(cl, on="target_acc", how="left")
c = f[f.PDS_acc.notna()].copy(); sizes = c.PDS_acc.map(c.PDS_acc.value_counts()); c = c[sizes >= 2].reset_index(drop=True)   # clusters with at least two isolates in the frame
c["cid"] = pd.factorize(c.PDS_acc)[0]; ncl = c.cid.max() + 1; size = np.bincount(c.cid)
print(f"{len(f):,} isolates in the frame; {len(c):,} of them in {ncl:,} clusters with 2 or more isolates in the frame (median size {int(np.median(size))}, largest {size.max()})")
sec = c.sector.map({s: i for i, s in enumerate(order)}).values; cid = c.cid.values

def counts(sec_vec):
    return np.stack([np.bincount(cid, weights=(sec_vec == k).astype(float), minlength=ncl) for k in range(4)])      # clusters x sectors
def metrics(n):
    H, F, A, E = n; nonh = F + A + E
    return dict(mixed_any=int(((H > 0) & (nonh > 0)).sum()), mixed_human_food=int(((H > 0) & (F > 0)).sum()), mixed_human_animal=int(((H > 0) & (A > 0)).sum()),
                human_in_mixed=int(H[(H > 0) & (nonh > 0)].sum()))
obs = metrics(counts(sec))

# permutation null: shuffle sector labels within serovar x year strata
strata = pd.factorize(c.serovar_clean.astype(str) + "|" + c.year.astype(int).astype(str))[0]; groups = [np.where(strata == s)[0] for s in range(strata.max() + 1)]
null = []
for _ in range(NPERM):
    p = sec.copy()
    for g in groups:
        if len(g) > 1: p[g] = sec[rng.permutation(g)]
    null.append(metrics(counts(p)))
null = pd.DataFrame(null); summ = pd.DataFrame({"observed": pd.Series(obs), "null_mean": null.mean(), "null_lo": null.quantile(0.025), "null_hi": null.quantile(0.975)})
summ["ratio"] = summ.observed / summ.null_mean; summ["p_ge"] = [(null[k] >= obs[k]).mean() for k in summ.index]; summ["p_le"] = [(null[k] <= obs[k]).mean() for k in summ.index]; print(summ.round(3).to_string()); summ.to_csv("table_cluster_null.csv")

# the same test for small clusters only (2 to 10 isolates), where chance mixing is least dominant
small = size[cid] <= 10; sub = np.where(small)[0]
def counts_sub(sv): return np.stack([np.bincount(cid[sub], weights=(sv[sub] == k).astype(float), minlength=ncl) for k in range(4)])
def m_small(n): H, F, A, E = n; keep = (H + F + A + E) > 0; nonh = F + A + E; return int((keep & (H > 0) & (nonh > 0)).sum())
o2 = m_small(counts_sub(sec)); n2 = []
for _ in range(NPERM // 2):
    p = sec.copy()
    for g in groups:
        if len(g) > 1: p[g] = sec[rng.permutation(g)]
    n2.append(m_small(counts_sub(p)))
print(f"clusters of 2-10 isolates: observed mixed human/non-human clusters {o2}, null mean {np.mean(n2):.1f} (95% {np.quantile(n2, .025):.0f}-{np.quantile(n2, .975):.0f}); n small clusters {len(set(cid[sub])):,}")

# ---------------------------------------------------------------- candidate cross-sector resistance lineages
# clusters of 4 to 300 isolates with at least one human and one food or animal isolate, in which EVERY isolate carries the same determinants for 3 or more drug classes
from amr9 import load_catalog, classes_of
cat = load_catalog(); rows = []
for k, g in c.groupby("cid"):
    if not (4 <= len(g) <= 300): continue
    s_ = g.sector.value_counts()
    if "Human" not in s_ or (s_.get("Food", 0) + s_.get("Animal", 0)) == 0: continue
    sets = g.genes.map(lambda x: frozenset(x.split(";")) if x else frozenset()); common = frozenset.intersection(*sets)
    cc = [x for x in classes_of(sorted(common), cat) if x not in ("EFFLUX", "MULTIDRUG")]
    if len(cc) >= 3: rows.append(dict(cluster=g.PDS_acc.iloc[0], n=len(g), serovar=g.serovar_clean.mode().iloc[0], **{x: int(s_.get(x, 0)) for x in order}, n_shared_genes=len(common), shared_classes=len(cc),
                                        years=f"{int(g.year.min())}-{int(g.year.max())}", shared=";".join(sorted(common))))
big = pd.DataFrame(rows).sort_values("n", ascending=False)
n_mixed4 = sum(1 for k, g in c.groupby("cid") if 4 <= len(g) <= 300 and "Human" in set(g.sector) and len(set(g.sector) - {"Human"}) > 0)
print(f"\nclusters of 4-300 isolates with a human and a non-human isolate: {n_mixed4:,}; with a shared genotype covering 3 or more drug classes in every isolate: {len(big):,}")
big.head(40).to_csv("table_mixed_clusters.csv", index=False); big.to_csv("table_mixed_clusters_all.csv", index=False); print(big.head(8)[["cluster", "n", "serovar", "Human", "Food", "Animal", "Environment", "shared_classes", "years", "shared"]].to_string(index=False))
print("serovars among these clusters:", big.serovar.value_counts().head(6).to_dict())

# ---------------------------------------------------------------- Figure 6
fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
labs = {"mixed_any": "Human + any\nnon-human", "mixed_human_food": "Human + food", "mixed_human_animal": "Human + animal"}
for i, k in enumerate(labs):
    r = summ.loc[k]; ax[0].errorbar(i, r.null_mean, yerr=[[r.null_mean - r.null_lo], [r.null_hi - r.null_mean]], fmt="s", color="#a0aec0", capsize=5, ms=7, label="Permutation null (95%)" if i == 0 else None)
    ax[0].scatter(i, r.observed, color="#c53030", s=60, zorder=4, label="Observed" if i == 0 else None); ax[0].text(i + 0.08, r.observed, f"{int(r.observed):,}", fontsize=8.5, va="center")
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labs.values()); ax[0].set_xlim(-0.5, 2.7); ax[0].set_ylim(0, 6800); ax[0].set(ylabel="Clusters containing both sectors", title="Mixed clusters: observed vs expected by chance"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper right")
ex = big.head(12).iloc[::-1]; left = np.zeros(len(ex))
for s in order: ax[1].barh(range(len(ex)), ex[s], left=left, color=COL[s], label=s); left += ex[s].values
ax[1].set_yticks(range(len(ex))); ax[1].set_yticklabels([f"{r.serovar} {r.cluster[3:].lstrip(chr(48))}, {r.years} ({r.shared_classes} classes)" for r in ex.itertuples()], fontsize=8)
ax[1].set(xlabel="Isolates in the cluster", title="Largest mixed clusters sharing 3 or more resistance classes"); ax[1].legend(frameon=False, fontsize=8.5, loc="lower right")
fig.tight_layout(); fig.savefig("figs9/fig8_mixed_clusters.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- co-occurrence of resistance genes
from collections import Counter                                  # str.get_dummies on 387,000 rows builds a dense table of every gene symbol (about 2 GB); counting first and keeping only the common genes needs about 12 MB
toks = f.genes.map(lambda t: set(t.split(";")) - {""}); cnt = Counter(g for ss in toks for g in ss)
X = pd.DataFrame({g: toks.map(lambda ss, g=g: g in ss).astype(np.uint8) for g in sorted(g for g, n in cnt.items() if n / len(f) >= 0.005)}); print(f"\n{X.shape[1]} determinants present in at least 0.5% of the {len(X):,} isolates")
phi = np.corrcoef(X.values.T.astype(np.float32)); np.fill_diagonal(phi, 0); G = nx.Graph()
for i, a in enumerate(X.columns):
    G.add_node(a, prev=float(X[a].mean()))
    for j in range(i + 1, len(X.columns)):
        if phi[i, j] >= 0.5: G.add_edge(a, X.columns[j], w=float(phi[i, j]))
G.remove_nodes_from([n for n in list(G) if G.degree(n) == 0]); comms = [sorted(cm) for cm in nx.community.louvain_communities(G, weight="w", seed=1, resolution=1.0) if len(cm) >= 3]
comms.sort(key=lambda cm: -sum(G.nodes[g]["prev"] for g in cm)); mod = []
for k, cm in enumerate(comms):
    have = X[cm].sum(axis=1) >= max(3, int(0.6 * len(cm))); sub = f[have.values]
    mod.append(dict(module=k + 1, genes=";".join(cm), n_genes=len(cm), carriers=int(have.sum()), pct_of_isolates=round(100 * have.mean(), 2), top_serovar=sub.serovar_clean.mode().iloc[0] if len(sub) else "", top_serovar_share=round(100 * (sub.serovar_clean == sub.serovar_clean.mode().iloc[0]).mean(), 1) if len(sub) else 0,
                    **{f"pct_{s.lower()}": round(100 * (sub.sector == s).mean(), 1) for s in order}))
mod = pd.DataFrame(mod); mod.to_csv("table_modules.csv", index=False); print(mod.drop(columns="genes").to_string(index=False)); [print(f"module {r.module}: {r.genes}") for r in mod.itertuples()]
pal = plt.get_cmap("tab10"); node_col = {g: pal(k) for k, cm in enumerate(comms) for g in cm}; keep = [n for n in G if n in node_col]; H = G.subgraph(keep)
centres = [(0, 0), (4.6, 0), (0, -4.4), (4.6, -4.4), (9.2, 0), (9.2, -4.4)]; pos, lab = {}, {}
for k, cm in enumerate(comms):                                    # each module on its own circle, so no label sits on another
    cx, cy = centres[k]; ang = np.linspace(0, 2 * np.pi, len(cm), endpoint=False) + np.pi / 2
    for g, t in zip(cm, ang): pos[g] = (cx + 1.25 * np.cos(t), cy + 1.25 * np.sin(t)); lab[g] = (cx + 1.62 * np.cos(t), cy + 1.62 * np.sin(t))
fig, ax = plt.subplots(figsize=(10.4, 6.6))
nx.draw_networkx_edges(H, pos, ax=ax, width=[0.5 + 2.4 * H[u][v]["w"] for u, v in H.edges], edge_color="#cbd5e0")
nx.draw_networkx_nodes(H, pos, ax=ax, node_color=[node_col[n] for n in H], node_size=[150 + 3600 * H.nodes[n]["prev"] for n in H], edgecolors="white", linewidths=0.8)
for g, (x, y_) in lab.items(): ax.text(x, y_, g, fontsize=8.4, ha="center", va="center")
for k, r in enumerate(mod.itertuples()):
    cx, cy = centres[k]; ax.text(cx, cy + 2.15, f"Module {r.module}: {r.pct_of_isolates:.1f}% of isolates, {r.top_serovar_share:.0f}% {r.top_serovar}", ha="center", fontsize=9, fontweight="bold", color=pal(k))
ax.set_xlim(-2.4, 11.6); ax.set_ylim(-6.9, 2.9); ax.set_aspect("equal"); ax.axis("off"); ax.set_title("Resistance determinants that co-occur (phi at least 0.5); node size = prevalence", fontsize=10.5)
fig.tight_layout(); fig.savefig("figs9/fig9_cooccurrence.png", dpi=200); plt.close(fig)
print("figures written")
