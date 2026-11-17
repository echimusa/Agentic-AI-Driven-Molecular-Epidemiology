"""
Chapter 5, Hands-on 4 - Rooting and the temporal signal in a tree (root-to-tip regression).
An ML tree from IQ-TREE is unrooted, and a root has to be placed before a tree can be read as a history. Every position along every branch is tried as the root; for each, the distance from the root to every tip is regressed on the
sampling date, d_i = mu * (t_i - t_0) + error, and the root with the largest correlation (positive slope) is kept, as in TempEst. The same regression is computed for the midpoint root, a root at a random tip, and for 1,000 date permutations
to ask whether the temporal signal could arise by chance. Data: the 195 Ebola virus genomes (collection dates) and the 60 H3N2 HA sequences (year of collection).
Output: table_root_to_tip.csv, Figure 4.
"""
import numpy as np, pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
def analyse(tf, dates, L, nperm=1000, seed=3):
    tree = Phylo.read(tf, "newick"); cl = list(tree.find_clades()); ids = {c: i for i, c in enumerate(cl)}; edges = []
    for c in cl:
        for ch in c.clades: edges.append((ids[c], ids[ch], max(ch.branch_length or 0.0, 0.0)))
    N = len(cl); r = [e[0] for e in edges]; q = [e[1] for e in edges]; w = [e[2] + 1e-12 for e in edges]; G = coo_matrix((w + w, (r + q, q + r)), shape=(N, N)).tocsr(); SP = shortest_path(G, directed=False)
    tips = [c for c in cl if c.is_terminal()]; tn = [t.name for t in tips]; ti = [ids[t] for t in tips]; t = np.array([dates[n] for n in tn]); cand = []; Rm = []; where = []
    for (u, v, l) in edges:
        side = SP[u, ti] < SP[v, ti]; steps = np.linspace(0, 1, 11) if l > 1e-9 else [0.0]
        for f in steps: Rm.append(np.where(side, SP[u, ti] + f * l, SP[v, ti] + (1 - f) * l)); where.append((u, v, f, l))
    Rm = np.array(Rm); Z = (Rm - Rm.mean(1, keepdims=True)) / np.maximum(Rm.std(1, keepdims=True), 1e-12); tz = (t - t.mean()) / t.std(); corr = Z @ tz / len(t)
    best = int(np.argmax(corr)); mid = int(np.argmin(Rm.max(1))); rng = np.random.default_rng(seed); rt = rng.choice(len(ti)); tipc = [k for k, (u, v, f, l) in enumerate(where) if v == ti[rt] and f == 1.0]
    def fitline(k): x = Rm[k]; b, a = np.polyfit(t, x, 1); return dict(r2=np.corrcoef(x, t)[0, 1] ** 2, slope=b, t0=-a / b if b != 0 else np.nan, n=len(t))
    P = rng.permuted(np.tile(tz, (nperm, 1)), axis=1).T; null = (Z @ P / len(t)).max(0); p = (1 + (null >= corr[best]).sum()) / (nperm + 1)
    out = {"best": fitline(best), "midpoint": fitline(mid), "random tip": fitline(tipc[0]) if tipc else None}
    for k, v in out.items():
        if v: v["subs_per_genome_year"] = v["slope"] * L
    return out, p, corr, null, Rm[best], t, tn
rows = []; res = {}
eb = pd.read_csv("ebola_dates.csv"); ed = {r.acc: pd.Timestamp(r.date).year + (pd.Timestamp(r.date).dayofyear - 1) / 365.25 for r in eb.itertuples()}
hm = pd.read_csv("h3n2_meta.csv"); hd = {r.accession: r.year + 0.5 for r in hm.itertuples()}
for name, tf, dates, L in [("Ebola virus (195 genomes)", "ebola.treefile", ed, 18961), ("H3N2 HA (60 sequences)", "h3n2.treefile", hd, 1762)]:
    out, p, corr, null, Rb, t, tn = analyse(tf, dates, L); res[name] = (out, p, corr, null, Rb, t)
    for k, v in out.items():
        if v: rows.append(dict(dataset=name, rooting=k, r2=v["r2"], slope_per_site_per_year=v["slope"], substitutions_per_genome_year=v["subs_per_genome_year"], root_date=v["t0"], permutation_p=p if k == "best" else np.nan, candidate_roots=len(corr)))
    print(f"{name}: candidate roots {len(corr)}; best-root R2 {out['best']['r2']:.3f}, slope {out['best']['slope']:.2e} per site per year = {out['best']['subs_per_genome_year']:.2f} per genome per year, root date {out['best']['t0']:.2f}; permutation p = {p:.4f}; null median max-R2 {np.median(null ** 2):.3f}, 95th percentile {np.percentile(null ** 2, 95):.3f}")
    print(f"   midpoint root: R2 {out['midpoint']['r2']:.3f}, rate {out['midpoint']['subs_per_genome_year']:.2f}, root date {out['midpoint']['t0']:.2f}" + (f"; random tip root: R2 {out['random tip']['r2']:.3f}, rate {out['random tip']['subs_per_genome_year']:.2f}, root date {out['random tip']['t0']:.2f}" if out['random tip'] else "") + f"; R2 range over all candidate roots {np.nanmin(np.where(corr > 0, corr ** 2, np.nan)):.3f} to {np.nanmax(corr ** 2):.3f}")
T = pd.DataFrame(rows); T.round(5).to_csv("table_root_to_tip.csv", index=False)
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.6)); n1 = "Ebola virus (195 genomes)"; out, p, corr, null, Rb, t = res[n1]; ax[0].scatter(t, Rb, s=14, color="#2b6cb0", alpha=0.7); xs = np.linspace(t.min(), t.max(), 50); f = out["best"]; ax[0].plot(xs, f["slope"] * (xs - f["t0"]), color="#c53030", lw=2)
ax[0].ticklabel_format(useOffset=False); ax[0].set(xlabel="Sampling date", ylabel="Root-to-tip distance (substitutions per site)", title=f"Ebola virus: R\u00b2 {f['r2']:.2f}, {f['subs_per_genome_year']:.0f} substitutions per genome-year")
for nm, c in [(n1, "#2b6cb0"), ("H3N2 HA (60 sequences)", "#dd6b20")]: o, pp, cc, nu, Rr, tt = res[nm]; ax[1].hist(nu ** 2, bins=30, color=c, alpha=0.45, density=True, label=f"{nm.split(' (')[0]}: permuted dates"); ax[1].axvline(o["best"]["r2"], color=c, lw=2)
ax[1].set(xlabel="R\u00b2 of the best root", ylabel="Density", title="Observed R\u00b2 (lines) against date permutations"); ax[1].legend(frameon=False, fontsize=8)
o, pp, cc, nu, Rr, tt = res["H3N2 HA (60 sequences)"]; ax[2].scatter(tt + np.random.default_rng(1).uniform(-0.25, 0.25, len(tt)), Rr, s=16, color="#dd6b20", alpha=0.8); xs = np.linspace(tt.min(), tt.max(), 50); f = o["best"]; ax[2].plot(xs, f["slope"] * (xs - f["t0"]), color="#c53030", lw=2); ax[2].set(xlabel="Year of collection (jittered)", ylabel="Root-to-tip distance (substitutions per site)", title=f"H3N2 HA: R\u00b2 {f['r2']:.2f}")
fig.tight_layout(); fig.savefig("figs5/fig4_root_to_tip.png", dpi=200); plt.close(fig); print("figure written")
