"""
Chapter 13, Hands-on 2 - Growth rate from genomes, agreement with case counts, and the reproduction number.
Public data: Ebola virus genomes of NCBI BioProject PRJNA257197 collected in Sierra Leone, May to September 2014. Steps: download, align (MAFFT), tree (FastTree GTR), time-scale (TreeTime),
then fit an exponential-growth coalescent to the dated tree and compare the growth rate with the case-based one from Hands-on 1.
Output: table_genome_growth.csv, Figure 2.   Needs mafft, FastTree, treetime, biopython, scipy.
"""
import os, io, json, time, subprocess, urllib.parse, urllib.request
import numpy as np, pandas as pd
from scipy import optimize, stats
from Bio import SeqIO, Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ebolacore import r_to_R
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for a in range(5):
        try:
            with urllib.request.urlopen(u, timeout=180) as r: return r.read().decode()
        except Exception: time.sleep(4 * (a + 1))
if not (os.path.exists("ebola.fasta") and os.path.exists("ebola_meta.csv")):
    ids = json.loads(get(E + "esearch.fcgi?db=nuccore&retmode=json&retmax=500&term=" + urllib.parse.quote("PRJNA257197[BioProject]")))["esearchresult"]["idlist"]; rows, fasta = [], ""
    for i in range(0, len(ids), 50):
        ch = ids[i:i + 50]; s = json.loads(get(E + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ch)))["result"]
        for u in s["uids"]:
            r = s[u]; f = dict(zip(r["subtype"].split("|"), r["subname"].split("|"))); rows.append(dict(acc=r["accessionversion"], isolate=f.get("isolate"), date=f.get("collection_date"), place=f.get("country")))
        fasta += get(E + "efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ch)); time.sleep(0.5)
    open("ebola.fasta", "w").write(fasta); pd.DataFrame(rows).to_csv("ebola_meta.csv", index=False)
meta = pd.read_csv("ebola_meta.csv").set_index("acc"); seqs = {r.id: r for r in SeqIO.parse("ebola.fasta", "fasta")}
meta["date"] = pd.to_datetime(meta.date, format="%d-%b-%Y", errors="coerce"); meta["pid"] = meta.isolate.str.extract(r"Makona-(G\d+)")[0]; meta = meta.dropna(subset=["pid", "date"])
meta["nfrac"] = [sum(c not in "ACGT" for c in str(seqs[a].seq).upper()) / len(seqs[a]) for a in meta.index]; meta = meta[meta.nfrac < 0.02].sort_values(["pid", "nfrac"]).drop_duplicates("pid").sort_values("date")
print(f"{len(meta)} genomes from {meta.pid.nunique()} patient identifiers; collected {meta.date.min():%d %b} to {meta.date.max():%d %b %Y}")
with open("ebola_clean.fasta", "w") as f: SeqIO.write([seqs[a] for a in meta.index], f, "fasta")
meta.reset_index()[["acc", "date"]].assign(date=lambda d: d.date.dt.strftime("%Y-%m-%d")).to_csv("ebola_dates.csv", index=False)

# ---------------------------------------------------------------- align, tree, time-scale
t0 = time.time()
if not os.path.exists("ebola_aln.fasta"): subprocess.run("mafft --auto --thread 1 ebola_clean.fasta > ebola_aln.fasta 2> mafft13.log", shell=True, check=True)
if not os.path.exists("ebola_ml.nwk"): subprocess.run("FastTree -nt -gtr -quiet ebola_aln.fasta > ebola_ml.nwk", shell=True, check=True)
print(f"alignment and tree: {time.time() - t0:.0f} s")
from treetime import TreeTime
from treetime.utils import parse_dates
if not os.path.exists("ebola_node_dates.csv"):
    dates = parse_dates("ebola_dates.csv", date_col="date", name_col="acc") if False else {r.acc: float(r.date[:4]) + (pd.Timestamp(r.date).dayofyear - 1) / 365.25 for r in pd.read_csv("ebola_dates.csv").itertuples()}
    tt = TreeTime(dates=dates, tree="ebola_ml.nwk", aln="ebola_aln.fasta", gtr="JC69", verbose=0); tt.clock_filter(reroot="least-squares", n_iqd=3, plot=False)
    tt.run(root="least-squares", infer_gtr=True, max_iter=2, resolve_polytomies=True, branch_length_mode="joint", time_marginal=False)
    tt.run(infer_gtr=False, max_iter=1, time_marginal=True, resolve_polytomies=False)                       # marginal node dates: posterior means, not ties at zero-length branches
    Phylo.write(tt.tree, "ebola_time.nwk", "newick"); pd.DataFrame([(n.name, (n.up.name if n.up is not None else ""), n.numdate, n.is_terminal(), len(n.clades)) for n in tt.tree.find_clades()], columns=["name", "parent", "numdate", "tip", "n_children"]).to_csv("ebola_node_dates.csv", index=False)
    json.dump(dict(rate=float(tt.date2dist.clock_rate), root=float(tt.tree.root.numdate)), open("ebola_time.json", "w")); print(f"clock rate {tt.date2dist.clock_rate:.2e} per site per year; root {tt.tree.root.numdate:.2f}")
print(f"time-scaling: {time.time() - t0:.0f} s")

# ---------------------------------------------------------------- exponential-growth coalescent on the dated tree
nd = pd.read_csv("ebola_node_dates.csv", keep_default_na=False).set_index("name"); root_date = json.load(open("ebola_time.json"))["root"]; DAY = 365.25
tips = nd[nd.tip].copy(); internal = nd[~nd.tip].copy(); tip_t = tips.numdate.values * DAY; Tmax = tip_t.max()
def events(floor_days=0.0):
    """Sampling events (+1) and coalescent events (-(children-1)) in time before the last sample; coalescences closer than floor_days are spread apart so that ties do not drive the fit."""
    ev = [(Tmax - x, +1) for x in tip_t]; co = sorted([(Tmax - r.numdate * DAY, r.n_children - 1) for r in internal.itertuples()])
    last = -1e9
    for t, m in co:
        for _ in range(m): t = max(t, last + floor_days); ev.append((t, -1)); last = t
    ev.sort(key=lambda e: (e[0], -e[1])); return ev
ev = events(0.0); ties = sum(1 for i in range(1, len(ev)) if ev[i][0] == ev[i - 1][0] and ev[i][1] < 0 and ev[i - 1][1] < 0)
print(f"{len(tips)} tips sampled over {(tip_t.max() - tip_t.min()):.0f} days; {len(internal)} internal nodes; most recent common ancestor {root_date:.2f}; coalescent events tied in time: {ties}")
def loglik(r, logN0, ev=None):
    ev = ev or globals()["ev"]; N0 = np.exp(logN0); k = 0; ll = 0.0; tprev = 0.0
    for t, d in ev:
        if k >= 2 and t > tprev:
            hz = k * (k - 1) / 2 / N0 * (((np.exp(r * t) - np.exp(r * tprev)) / r) if abs(r) > 1e-9 else (t - tprev)); ll -= hz
        if d < 0: ll += np.log(k * (k - 1) / 2) - logN0 + r * t; k -= 1
        else: k += d
        tprev = t
    return ll
FLOORS = (0, 0.5, 1, 2, 5); rr = []
for fl in FLOORS:
    e_ = events(fl); f_ = optimize.minimize(lambda p: -loglik(p[0], p[1], e_), x0=[0.02, np.log(50)], method="L-BFGS-B", bounds=[(-0.01, 0.12), (np.log(0.5), np.log(5000))]); rr.append(f_.x[0]); print(f"floor of {fl} days between coalescent events: r = {f_.x[0]:.4f} per day (N0 = {np.exp(f_.x[1]):.0f} days)")
rc, sec = np.load("r_cases.npy"); print(f"cases (1 Jun to 15 Sep 2014): r = {rc:.4f} (SE {sec:.4f}); genomic r ranges from {min(rr):.4f} to {max(rr):.4f} depending only on how tied coalescent events are spread")
mrca = pd.Timestamp(f"{int(root_date)}-01-01") + pd.Timedelta(days=(root_date - int(root_date)) * 365.25); cases = pd.read_csv("ebola_sl_linelist.csv", parse_dates=["date_of_onset"]); first = cases.date_of_onset.min()
upto = (cases.date_of_onset <= tips.numdate.max().__float__().__class__(0) and cases.date_of_onset <= pd.Timestamp("2014-09-28")).sum() if False else (cases.date_of_onset <= pd.Timestamp("2014-09-28")).sum()
print(f"time of the most recent common ancestor: {mrca:%d %b %Y} (decimal year {root_date:.2f}); first case onset in the line list: {first:%d %b %Y}; cases with onset up to 28 Sep 2014: {upto:,}; genomes: {len(tips)} ({100 * len(tips) / upto:.1f}% of those cases)")
print(f"resolution of the tree: {len(tips)} tips, {len(internal)} internal nodes (a fully resolved binary tree would have {len(tips) - 1})")
rows = [dict(coalescent_floor_days=fl, r_per_day=r_, doubling_days=np.log(2) / r_ if r_ > 0 else np.nan) for fl, r_ in zip(FLOORS, rr)]; T = pd.DataFrame(rows); T.round(4).to_csv("table_genome_growth.csv", index=False)

# ---------------------------------------------------------------- Figure 2
kids = {}; [kids.setdefault(r.parent, []).append(n) for n, r in nd.iterrows() if r.parent != ""]; leaves = []
def collect(n):
    if nd.loc[n, "tip"]: leaves.append(n)
    else: [collect(c) for c in kids.get(n, [])]
root_name = nd.index[nd.parent == ""][0]; collect(root_name); order = {n: i for i, n in enumerate(leaves)}; ypos = {}
def ycoord(n):
    ypos[n] = order[n] if nd.loc[n, "tip"] else np.mean([ycoord(c) for c in kids[n]]); return ypos[n]
ycoord(root_name); mon = {a: meta.loc[a, "date"].month for a in meta.index}; cm = plt.get_cmap("viridis")
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [1.1, 1, 1]})
for n, r in nd.iterrows():
    if r.parent != "": ax[0].plot([nd.loc[r.parent, "numdate"], r.numdate], [ypos[n]] * 2, color="#a0aec0", lw=0.6)
    if not r.tip: ys = [ypos[c] for c in kids[n]]; ax[0].plot([r.numdate] * 2, [min(ys), max(ys)], color="#a0aec0", lw=0.6)
for n in leaves: ax[0].scatter(nd.loc[n, "numdate"], ypos[n], s=7, color=cm((mon.get(n, 7) - 5) / 5))
ax[0].set(xlabel="Year", title=f"Time-scaled tree of {len(tips)} genomes", xlim=(2014.3, 2014.85)); ax[0].ticklabel_format(useOffset=False, axis="x"); ax[0].set_yticks([]); ax[0].spines["left"].set_visible(False)
ax[1].axhline(0, color="black", lw=0.8); ax[1].axhspan(rc - 1.96 * sec, rc + 1.96 * sec, color="#fbd38d", alpha=0.55, label="Cases: 95% interval, 1 Jun to 15 Sep"); ax[1].plot(range(len(FLOORS)), rr, "o-", color="#2b6cb0", label="Genomes: coalescent fit"); ax[1].set_xticks(range(len(FLOORS))); ax[1].set_xticklabels([str(f) for f in FLOORS]); ax[1].legend(frameon=False, fontsize=8.5, loc="upper right")
ax[1].set(xlabel="Minimum spacing between tied coalescent events (days)", ylabel="Growth rate r (per day)", title="The genomic estimate depends on how ties are treated")
dd = cases.groupby("date_of_onset").size().asfreq("D", fill_value=0)["2014-05-15":"2014-10-15"]; ax[2].bar(dd.index, dd.values, width=1, color="#a0aec0", label="Cases by onset date"); gd = meta.date.value_counts().sort_index(); ax[2].vlines(gd.index, 0, -gd.values * 1.5, color="#2b6cb0", lw=1.2, label="Genomes sampled (downward)"); ax[2].axvline(mrca, color="#c53030", ls="--", lw=1.2, label="Common ancestor of the genomes")
ax[2].set(ylabel="Cases per day (up), genomes per day x 1.5 (down)", title="Genomes cover the first four months"); ax[2].yaxis.label.set_size(8.5); ax[2].legend(frameon=False, fontsize=8, loc="upper left"); ax[2].tick_params(axis="x", rotation=30)
fig.tight_layout(); fig.savefig("figs13/fig2_genome_growth.png", dpi=200); plt.close(fig); print("figure written")
