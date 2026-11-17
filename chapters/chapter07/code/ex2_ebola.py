"""
Chapter 7, Exercise 2 - What do genomes add? Ebola virus (Makona) genomes, Sierra Leone 2014.

Public data: 249 Ebola virus genomes in NCBI BioProject PRJNA257197, collected in Sierra
Leone between 25 May and 28 Sep 2014 with day-level collection dates. No true infectors
are recorded for these genomes, so this exercise measures how much the genomes narrow
down the candidate infectors; it cannot score accuracy (Exercise 1 does that).
"""
import io, os, json, time, urllib.parse, urllib.request
import numpy as np, pandas as pd, networkx as nx, mappy
from Bio import SeqIO
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = "figs2"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
C0, C1, C2 = "#2b6cb0", "#c05621", "#718096"

# ---------------------------------------------------------------- 1. download
def get(url):
    with urllib.request.urlopen(url, timeout=180) as r: return r.read().decode()

def download(project="PRJNA257197", cache="ebola.fasta", meta_cache="ebola_meta.csv"):
    if os.path.exists(cache) and os.path.exists(meta_cache): return
    ids = json.loads(get(E + "esearch.fcgi?db=nuccore&retmode=json&retmax=500&term="
                         + urllib.parse.quote(f"{project}[BioProject]")))["esearchresult"]["idlist"]
    rows, fasta = [], ""
    for i in range(0, len(ids), 50):
        ch = ids[i:i + 50]
        s = json.loads(get(E + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ch)))["result"]
        for u in s["uids"]:
            r = s[u]; keys, vals = r["subtype"].split("|"), r["subname"].split("|")
            f = dict(zip(keys, vals))
            rows.append(dict(acc=r["accessionversion"], isolate=f.get("isolate"), date=f.get("collection_date"),
                             place=f.get("country")))
        fasta += get(E + "efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ch)); time.sleep(0.5)
    open(cache, "w").write(fasta); pd.DataFrame(rows).to_csv(meta_cache, index=False)

download()
meta = pd.read_csv("ebola_meta.csv").set_index("acc")
seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse("ebola.fasta", "fasta")}
meta["date"] = pd.to_datetime(meta.date, format="%d-%b-%Y")
meta["pid"] = meta.isolate.str.extract(r"Makona-(G\d+)")[0]
meta = meta.dropna(subset=["pid"])
meta["n_frac"] = [sum(c not in "ACGT" for c in seqs[a]) / len(seqs[a]) for a in meta.index]
print(f"{len(meta)} genomes after dropping records without a patient identifier, collected {meta.date.min():%d %b} to {meta.date.max():%d %b %Y}")

# keep genomes with < 2% ambiguous bases, and the cleanest record per patient identifier
meta = meta[meta.n_frac < 0.02].sort_values(["pid", "n_frac"]).drop_duplicates("pid").sort_values("date")
print(f"analysis set: {len(meta)} genomes from {meta.pid.nunique()} distinct patient identifiers")

# ------------------------------------------- 2. align to a reference (reference-anchored)
ref_id = max(meta.index, key=lambda a: len(seqs[a])); ref = seqs[ref_id]; L = len(ref)
aligner = mappy.Aligner(seq=ref, preset="asm5")
code = {"A": 0, "C": 1, "G": 2, "T": 3}
def to_row(s):
    row = np.full(L, 4, dtype=np.int8)                        # 4 = missing
    for h in aligner.map(s):
        q, r = h.q_st, h.r_st
        for n, op in h.cigar:
            if op in (0, 7, 8):                               # match / mismatch
                row[r:r + n] = [code.get(c, 4) for c in s[q:q + n]]; q += n; r += n
            elif op == 1: q += n                              # insertion in the query
            elif op == 2: r += n                              # deletion in the query
    return row
A = np.vstack([to_row(seqs[a]) for a in meta.index])
called = A < 4
var = [j for j in range(L) if len(np.unique(A[:, j][called[:, j]])) > 1]
V = A[:, var]
print(f"alignment {A.shape[0]} x {L} nt, {len(var)} variable sites")

def snp_matrix(V):
    n = len(V); D = np.zeros((n, n), dtype=int)
    for i in range(n):
        ok = (V[i] < 4) & (V < 4)
        D[i] = ((V[i] != V) & ok).sum(axis=1)
    return D
D = snp_matrix(V)
iu = np.triu_indices(len(D), 1)
print(f"pairwise SNPs: median {np.median(D[iu]):.0f}, {np.mean(D[iu] == 0):.1%} of pairs identical, max {D.max()}")

# ----------------------------- 3. molecular clock from root-to-tip regression (Chapter 6)
t_yr = ((meta.date - meta.date.min()).dt.days / 365.25).to_numpy()
first = np.argsort(t_yr)[:5]
root = np.array([np.bincount(V[first, j][V[first, j] < 4], minlength=4).argmax() for j in range(V.shape[1])])
rtt = ((V != root) & (V < 4)).sum(axis=1)
fit = stats.linregress(t_yr, rtt)
rate_site = fit.slope / np.mean(called.sum(axis=1))
print(f"root-to-tip: {fit.slope:.1f} SNPs/genome/year (SE {fit.stderr:.1f}), "
      f"{rate_site:.2e} subs/site/year, R2 = {fit.rvalue**2:.2f}")

# ------------------------------------------------ 4. SNP-threshold clusters (sensitivity)
def clusters(D, thr):
    G = nx.Graph(); G.add_nodes_from(range(len(D)))
    G.add_edges_from(zip(*np.where(np.triu(D <= thr, 1))))
    return sorted((len(c) for c in nx.connected_components(G)), reverse=True)
thr_tab = pd.DataFrame([{"threshold (SNPs)": t, "clusters": len(clusters(D, t)),
                         "largest cluster": clusters(D, t)[0],
                         "singletons": sum(c == 1 for c in clusters(D, t))} for t in [0, 1, 2, 3, 5, 10]])
print(thr_tab.to_string(index=False))

# ------------- 5. who could have infected whom? time-only vs time + genetics (early window)
early = np.where((meta.date >= "2014-05-25") & (meta.date <= "2014-06-18"))[0]
print(f"early window 25 May-18 Jun: {len(early)} genomes")
mu = max(fit.slope, 1e-3) / 365.25                        # SNPs per genome per day (clock above)
si_mean, si_sd = 15.3, 9.3                                # Ebola serial interval, WHO Ebola Response Team 2014
shape, scale = (si_mean / si_sd) ** 2, si_sd ** 2 / si_mean
day = (meta.date - meta.date.min()).dt.days.to_numpy()
def posterior(j, genetics):
    dt = day[j] - day[early]                              # collection-date gap (proxy for onset gap)
    w = np.where(dt >= 0, stats.gamma.pdf(np.maximum(dt, 0.5), shape, scale=scale), 0.0)
    w[early == j] = 0.0
    if genetics:
        lam = mu * (np.abs(dt) + si_mean)                 # expected SNPs: both branches since the split
        w = w * stats.poisson.pmf(D[j, early], lam)
    return w / w.sum() if w.sum() > 0 else None
rows = []
for j in early:
    pt, pg = posterior(j, False), posterior(j, True)
    if pt is None or pg is None: continue
    rows.append(dict(case=meta.pid.iloc[j], date=meta.date.iloc[j],
                     neff_time=1 / (pt ** 2).sum(), neff_gen=1 / (pg ** 2).sum(),
                     top=meta.pid.iloc[early[pg.argmax()]], p_top=pg.max(), p_top_time=pt.max()))
rank = pd.DataFrame(rows)
print(f"{len(rank)} cases with at least one earlier candidate infector")
print(f"effective number of candidate infectors: time only median {rank.neff_time.median():.1f}, "
      f"time + genetics median {rank.neff_gen.median():.1f}; "
      f"cases with a top candidate above 50%: {(rank.p_top > 0.5).sum()} (time only: {(rank.p_top_time > 0.5).sum()})")
rank.to_csv("ranking.csv", index=False)

# ------------------------------------------------------------------ figures
order = np.argsort(day)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
im = ax[0].imshow(np.minimum(D[np.ix_(order, order)], 12), cmap="viridis_r", vmin=0, vmax=12)
ax[0].set(title="Pairwise SNP distances, ordered by collection date", xlabel="Genome", ylabel="Genome")
cb = fig.colorbar(im, ax=ax[0], fraction=0.046); cb.set_label("SNPs (capped at 12)")
ax[1].hist(D[iu], bins=np.arange(-0.5, D[iu].max() + 1.5, 1), color=C0)
ax[1].set(xlabel="SNPs between two genomes", ylabel="Pairs", title=f"{np.mean(D[iu] == 0):.0%} of pairs are identical")
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_snp_matrix.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.scatter(meta.date, rtt + np.random.default_rng(1).uniform(-0.15, 0.15, len(rtt)), s=16, color=C0, alpha=0.7)
xs = np.array([t_yr.min(), t_yr.max()])
ax.plot(meta.date.min() + pd.to_timedelta(xs * 365.25, unit="D"), fit.intercept + fit.slope * xs, color=C1, lw=2)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.set(ylabel="SNPs from the early consensus", title=f"Root-to-tip: {fit.slope:.0f} SNPs per genome per year ({rate_site:.1e} per site)")
fig.autofmt_xdate(); fig.tight_layout(); fig.savefig(f"{OUT}/fig7_clock.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(thr_tab["threshold (SNPs)"], thr_tab["clusters"], "o-", color=C0, label="number of clusters")
ax.plot(thr_tab["threshold (SNPs)"], thr_tab["largest cluster"], "s-", color=C1, label="size of largest cluster")
ax.set(xlabel="SNP threshold", ylabel="Count", title="Cluster membership depends on the threshold"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig8_threshold.png", dpi=200); plt.close(fig)

# minimum spanning network of genome groups that are identical at jointly called sites (early window)
De = D[np.ix_(early, early)]
Z = nx.Graph(); Z.add_nodes_from(range(len(early))); Z.add_edges_from(zip(*np.where(np.triu(De == 0, 1))))
groups = [sorted(c) for c in nx.connected_components(Z)]               # members are positions within `early`
Dg = np.array([[De[np.ix_(g, h)].min() if g is not h else 0 for h in groups] for g in groups])
G = nx.Graph(); G.add_nodes_from(range(len(groups)))
for i in range(len(groups)):
    for k in range(i + 1, len(groups)): G.add_edge(i, k, weight=Dg[i, k])
T = nx.minimum_spanning_tree(G)
pos = nx.kamada_kawai_layout(T, weight=None)
mean_day = [np.mean([day[early[m]] for m in g]) for g in groups]
fig, ax = plt.subplots(figsize=(8.5, 6))
nx.draw_networkx_edges(T, pos, ax=ax, edge_color="#cbd5e0", width=1.5)
sc = ax.scatter([pos[i][0] for i in range(len(groups))], [pos[i][1] for i in range(len(groups))],
                s=[40 + 45 * len(g) for g in groups], c=mean_day, cmap="plasma", edgecolor="white", zorder=3)
for i, j, d_ in T.edges(data="weight"):
    ax.text((pos[i][0] + pos[j][0]) / 2, (pos[i][1] + pos[j][1]) / 2, int(d_), fontsize=8, color=C2, ha="center")
cb = fig.colorbar(sc, ax=ax, fraction=0.04); cb.set_label("Mean collection day (0 = 25 May 2014)")
ax.set_axis_off()
ax.set_title(f"Minimum spanning network of early genomes\n{len(groups)} groups from {len(early)} genomes; dot size = genomes, edge label = SNPs")
fig.tight_layout(); fig.savefig(f"{OUT}/fig9_network.png", dpi=200); plt.close(fig)
print(f"network: {len(groups)} groups of identical genomes; largest has {max(len(g) for g in groups)} genomes")

fig, ax = plt.subplots(figsize=(8, 4.2))
for col, c, lab in [("neff_time", C2, "timing only"), ("neff_gen", C1, "timing + genetics")]:
    v = np.sort(rank[col].to_numpy())
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=c, lw=2.2,
            label=f"{lab} (median {np.median(v):.0f})")
ax.set(xlabel="Effective number of candidate infectors per case", ylabel="Share of cases at or below",
       title="How much the genomes narrow the candidate infectors", ylim=(0, 1.02))
ax.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(1.0))
ax.legend(frameon=False, loc="lower right"); fig.tight_layout(); fig.savefig(f"{OUT}/fig10_candidates.png", dpi=200); plt.close(fig)
print("figures written to", OUT)
