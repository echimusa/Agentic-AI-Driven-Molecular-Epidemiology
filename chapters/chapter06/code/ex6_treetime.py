"""
Chapter 6, Hands-on 3 - Time-scaled trees with TreeTime on real genomes.
Data: the 195 Ebola virus genomes of Chapter 5 (alignment, collection dates, the ML tree of Chapter 5 with its support values removed) and the 327 SARS-CoV-2 genomes of Hands-on 1 plus the Wuhan-Hu-1 reference (aligned to the reference with minimap2, tree built with FastTree).
For each data set TreeTime (the command line of the Session 6 tutorial) is run several ways, and the root date (the common ancestor of the sampled genomes) is read from dates.tsv:
  default           treetime --tree T --aln A --dates D                                      (least-squares clock, outlier filter 3)
  covariation       ... --covariation --confidence                                           (rate standard error, 90% date intervals)
  no filter         ... --covariation --confidence --clock-filter 0                          (keep every genome)
  rate x0.5, x2     ... --clock-rate R --clock-std-dev 0.2R --confidence                    (an imposed rate, as when a rate is taken from another study)
  skyline           ... --coalescent skyline --n-skyline 8 --covariation                     (Ebola only)
and a date-randomisation test reruns `treetime clock` after shuffling the dates (50 shuffles for Ebola, 30 for SARS-CoV-2).
Output: table_treetime.csv, table_treetime_drt.csv, table_skyline_ebola.csv, Figure 3 (figs6/fig_ch6_timetree.png), the divergence trees of the default runs.
Needs: treetime (pip install phylo-treetime), FastTree, mappy, biopython, numpy, pandas, matplotlib. About 25 minutes from scratch; a second call reuses the runs in tt_runs and table_treetime_drt.csv unless RERUN=1.
"""
import os, re, sys, gzip, shutil, subprocess, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
NPERM = {"Ebola": 50, "SARS-CoV-2": 30}; L = {"Ebola": 18961, "SARS-CoV-2": 29903}
def make_sars2():
    """Reference-anchored alignment of the 327 genomes of Hands-on 1 plus the reference, written as sars2_aln.fasta.gz, with sars2_dates.csv."""
    import mappy
    from Bio import SeqIO
    tab = pd.read_csv("sars2_root_to_tip_table.csv"); seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse("sars2.fasta", "fasta")}; ref = seqs["NC_045512.2"]; code = {"A": 0, "C": 1, "G": 2, "T": 3}; aligner = mappy.Aligner(seq=ref, preset="asm5")
    def row(s):
        r_ = np.full(len(ref), 4, dtype=np.int8)
        for h in aligner.map(s):
            q, r = h.q_st, h.r_st
            for n, op in h.cigar:
                if op in (0, 7, 8): r_[r:r + n] = [code.get(c, 4) for c in s[q:q + n]]; q += n; r += n
                elif op == 1: q += n
                elif op == 2: r += n
        return "".join("ACGTN"[k] for k in r_)
    with gzip.open("sars2_aln.fasta.gz", "wt") as f:
        f.write(f">NC_045512.2\n{ref}\n")
        for a in tab.acc: f.write(f">{a}\n{row(seqs[a])}\n")
    pd.concat([pd.DataFrame({"acc": ["NC_045512.2"], "date": ["2019-12-26"]}), tab[["acc", "date"]]]).to_csv("sars2_dates.csv", index=False)
if not os.path.exists("sars2_aln.fasta.gz"): make_sars2()
if not os.path.exists("sars2.fasta.aln"):
    with gzip.open("sars2_aln.fasta.gz", "rt") as f, open("sars2.fasta.aln", "w") as g: shutil.copyfileobj(f, g)
if not os.path.exists("sars2_tree.nwk"):
    out = subprocess.run(["FastTree", "-nt", "-gtr", "-nosupport", "-quiet", "sars2.fasta.aln"], capture_output=True, text=True).stdout; open("sars2_tree.nwk", "w").write(out)
DATA = {"Ebola": ("ebola_aln.fasta", "ebola_tree.nwk", "ebola_dates.csv"), "SARS-CoV-2": ("sars2.fasta.aln", "sars2_tree.nwk", "sars2_dates.csv")}
def run(ds, tag, extra, outdir=None):
    aln, tree, dates = DATA[ds]; od = outdir or f"tt_runs/{ds.split('-')[0].lower()}_{tag}"
    if os.path.exists(f"{od}/dates.tsv") and os.environ.get("RERUN") != "1": return od                     # reuse an earlier run of the same analysis (RERUN=1 repeats it)
    shutil.rmtree(od, ignore_errors=True)
    subprocess.run(["treetime", "--tree", tree, "--aln", aln, "--dates", dates, "--name-column", "acc", "--date-column", "date", "--outdir", od, "--rng-seed", "1", *extra], capture_output=True, text=True)
    return od
def parse(od):
    mc = open(f"{od}/molecular_clock.txt").read(); rate = float(re.search(r"--rate:\s+([\d.eE+-]+)", mc).group(1)); sd = re.search(r"\+/-\s+([\d.eE+-]+)", mc); r2 = float(re.search(r"--r\^2:\s+([\d.]+)", mc).group(1))
    rows = [l.rstrip("\n").split("\t") for l in open(f"{od}/dates.tsv") if not l.startswith("#")]; nodes = [(r[0], float(r[2]), float(r[3]) if len(r) > 3 else np.nan, float(r[4]) if len(r) > 4 else np.nan) for r in rows if r[0].startswith("NODE_") and r[2] != "--"]
    root = min(nodes, key=lambda x: x[1]); out = [l.split("\t")[0] for l in open(f"{od}/outliers.tsv") if l.strip() and l.split("\t")[0]] if os.path.exists(f"{od}/outliers.tsv") else []        # the header line starts with a tab
    return dict(clock_rate=rate, rate_sd=float(sd.group(1)) if sd else np.nan, r2=r2, root_date=root[1], root_lower=root[2], root_upper=root[3], outliers=len(out), outlier_ids=";".join(out))
rows = []; os.makedirs("tt_runs", exist_ok=True); os.makedirs("figs6", exist_ok=True)
for ds in DATA:
    base = None
    for tag, extra in [("default", []), ("covariation", ["--covariation", "--confidence"]), ("nofilter", ["--covariation", "--confidence", "--clock-filter", "0"])]:
        od = run(ds, tag, extra); r = parse(od); rows.append(dict(dataset=ds, scenario=tag, **r))
        if tag == "default": shutil.copy(f"{od}/divergence_tree.nexus", f"divergence_tree_{ds.split('-')[0].lower()}.nexus")
        if tag == "covariation": base = r["clock_rate"]
    for mult in (0.5, 2.0):
        od = run(ds, f"rate{mult}", ["--clock-rate", f"{mult * base:.6e}", "--clock-std-dev", f"{0.2 * mult * base:.6e}", "--confidence"]); r = parse(od); r["clock_rate"] = mult * base; r["rate_sd"] = 0.2 * mult * base; rows.append(dict(dataset=ds, scenario=f"rate x{mult:g}", **r))
    if ds == "Ebola":
        od = run(ds, "skyline", ["--coalescent", "skyline", "--n-skyline", "8", "--covariation"]); rows.append(dict(dataset=ds, scenario="skyline", **parse(od)))
        sk = pd.read_csv(f"{od}/skyline.tsv", sep="\t", comment="#", header=None, names=["date", "Ne", "lower", "upper"]); sk.to_csv("table_skyline_ebola.csv", index=False)
T = pd.DataFrame(rows); T["rate_per_genome_year"] = T.clock_rate * T.dataset.map(L); T.round(8).to_csv("table_treetime.csv", index=False)
# ---- date randomisation with treetime clock
drt = []
for ds in (DATA if (os.environ.get("RERUN") == "1" or not os.path.exists("table_treetime_drt.csv")) else []):
    aln, tree, dates = DATA[ds]; d = pd.read_csv(dates); rng = np.random.default_rng(7); obs = None
    for k in range(NPERM[ds] + 1):
        dd = d.copy()
        if k > 0: dd["date"] = rng.permutation(d["date"].values)
        with tempfile.TemporaryDirectory() as tmp:
            dp = os.path.join(tmp, "d.csv"); dd.to_csv(dp, index=False); od = os.path.join(tmp, "o")
            subprocess.run(["treetime", "clock", "--tree", tree, "--aln", aln, "--dates", dp, "--name-column", "acc", "--date-column", "date", "--outdir", od, "--rng-seed", "1"], capture_output=True, text=True)
            try: mc = open(f"{od}/molecular_clock.txt").read(); rate = float(re.search(r"--rate:\s+([\d.eE+-]+)", mc).group(1)); r2 = float(re.search(r"--r\^2:\s+([\d.]+)", mc).group(1))
            except Exception: rate, r2 = np.nan, np.nan
        drt.append(dict(dataset=ds, replicate=k, rate=rate, r2=r2))
Dr = pd.DataFrame(drt) if drt else pd.read_csv("table_treetime_drt.csv"); Dr.round(8).to_csv("table_treetime_drt.csv", index=False)
# ---- figure
def todate(x): return 1970 + (pd.Timestamp(year=int(x), month=1, day=1).value / 1e9) / (365.2425 * 86400) if False else x
fig, ax = plt.subplots(2, 2, figsize=(11.6, 7.4))
for a, ds in zip(ax[0], ("Ebola", "SARS-CoV-2")):
    s = T[(T.dataset == ds) & (T.scenario != "skyline")].reset_index(drop=True)
    for i, r in s.iterrows():
        a.plot(r.root_date, -i, "o", color="#2b6cb0")
        if not np.isnan(r.root_lower): a.plot([r.root_lower, r.root_upper], [-i, -i], color="#2b6cb0", lw=2.2)
    a.set(yticks=-np.arange(len(s)), yticklabels=list(s.scenario), xlabel="Date of the common ancestor of the sampled genomes (year)", title=f"{ds}: root date by analysis"); a.invert_yaxis() if False else None
for a, ds, c in zip((ax[1][0],), ("Ebola",), ("#718096",)):
    d = Dr[Dr.dataset == ds]; o = float(d[d.replicate == 0].rate.iloc[0]); p = d[d.replicate > 0].rate.dropna(); a.hist(p * L[ds], bins=15, color=c, alpha=0.85); a.axvline(o * L[ds], color="#c05621", lw=2); a.set(xlabel="Clock rate (substitutions per genome per year)", ylabel="Shuffles", title=f"{ds}: {len(p)} of {len(d) - 1} shuffles gave a rate")
sk = pd.read_csv("table_skyline_ebola.csv"); ax[1][1].fill_between(sk.date, sk.lower, sk.upper, color="#bee3f8"); ax[1][1].plot(sk.date, sk.Ne, "o-", color="#2b6cb0"); ax[1][1].ticklabel_format(axis="x", useOffset=False, style="plain"); ax[1][1].set(yscale="log", xlabel="Year", ylabel="Effective population size (arbitrary scale)", title="Ebola: coalescent skyline from TreeTime")
fig.tight_layout(); fig.savefig("figs6/fig_ch6_timetree.png", dpi=180); print(T.round(5).to_string(index=False)); print(Dr.groupby("dataset").rate.describe().round(6))
