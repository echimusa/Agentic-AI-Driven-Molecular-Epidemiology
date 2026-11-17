"""
Chapter 6, Hands-on 5 - summaries of the BEAST 2 analyses (run_all.sh): posteriors, effective sample sizes, marginal likelihoods and Bayes factors, plus TreeTime on the same genomes.
Reads  xml/mcmc_het_strict_s1.log, _s2.log, xml/mcmc_het_relaxed.log, xml/out_ps_<model>.txt and xml/ps_<model>/step*/likelihood.log (written by BEAST 2),
writes table_beast_posterior.csv, table_beast_marginal.csv, table_beast_bf.csv, table_beast_treetime.csv, ps_step_samples.csv.gz and figs6/fig_ch6_beast.png.
Own implementations: effective sample size (initial positive sequence of autocorrelations), 95% highest posterior density, Gelman-Rubin statistic, path sampling (trapezoid) and stepping-stone estimators from the step samples.
USAGE   python ex6_beast_analysis.py          (needs treetime and the files of ex6_beast.py in the working folder)
"""
import glob, math, os, re, subprocess, shutil, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
os.makedirs("figs6", exist_ok=True)
L = 18961; ALPHA = 0.3; BURN = 0.2; THIN = 1
def ess_ips(x):
    """Geyer's initial positive sequence estimator of the effective sample size."""
    x = np.asarray(x, float); n = len(x); x = x - x.mean(); f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n]; ac = ac / ac[0]; tot = -1.0; k = 0
    while k + 1 < n:
        g = ac[k] + ac[k + 1]
        if g < 0: break
        tot += 2 * g; k += 2
    return n / max(1.0, tot)
def hpd(x, p=0.95):
    x = np.sort(np.asarray(x, float)); m = int(math.floor(p * len(x))); w = x[m:] - x[:len(x) - m]; i = int(np.argmin(w)); return x[i], x[i + m]
def rhat(chains):
    m = len(chains); n = min(len(c) for c in chains); c = np.array([np.asarray(c[:n], float) for c in chains]); w = c.var(axis=1, ddof=1).mean(); b = n * c.mean(axis=1).var(ddof=1); return math.sqrt(((n - 1) / n * w + b / n) / w)
last = max(pd.read_csv("ebola_dates.csv").date.map(pd.Timestamp).map(lambda d: d.year + (d.dayofyear - 0.5) / (366 if d.is_leap_year else 365)))
def readlog(p): d = pd.read_csv(p, sep="\t", comment="#"); return d.iloc[int(len(d) * 0.1):].reset_index(drop=True)
runs = {"het_strict": ["xml/mcmc_het_strict_s1.log", "xml/mcmc_het_strict_s2.log"], "het_relaxed": ["xml/mcmc_het_relaxed.log"]}
rows = []
for model, files in runs.items():
    files = [f for f in files if os.path.exists(f)]
    if not files: continue
    ch = [readlog(f) for f in files]; rate = "clockRate" if model == "het_strict" else "ucldMean"
    for c in ch: c["rootdate"] = last - c["TreeHeight"]; c["rate_genome"] = c[rate] * L; c["doubling_days"] = 365.25 * math.log(2) / c["growthRate"]
    pars = ["posterior", "likelihood", rate, "rate_genome", "TreeHeight", "rootdate", "growthRate", "doubling_days", "popSize", "kappa", "gammaShape"] + (["ucldStdev"] if model == "het_relaxed" else [])
    for p in pars:
        allv = np.concatenate([c[p].to_numpy() for c in ch]); lo, hi = hpd(allv); e = [ess_ips(c[p].to_numpy()) for c in ch]
        rows.append(dict(model=model, parameter=p, mean=allv.mean(), median=float(np.median(allv)), hpd_lower=lo, hpd_upper=hi, chains=len(ch), samples_per_chain=len(ch[0]), ess_min=min(e), ess_total=sum(e), rhat=rhat([c[p].to_numpy() for c in ch]) if len(ch) > 1 else np.nan))
P = pd.DataFrame(rows); P.round(6).to_csv("table_beast_posterior.csv", index=False)
# ---- marginal likelihoods by path sampling and stepping-stone sampling
def betas(k): return np.array([((k - 1 - i) / (k - 1)) ** (1 / ALPHA) for i in range(k)])
def ps_est(b, E): return float(sum((b[i] - b[i + 1]) * (E[i] + E[i + 1]) / 2 for i in range(len(b) - 1)))
def ss_est(b, S):
    tot = 0.0
    for i in range(len(b) - 1):
        x = (b[i] - b[i + 1]) * S[i + 1]; m = x.max(); tot += m + math.log(np.exp(x - m).mean())
    return float(tot)
samples, mrows = [], []
for out in sorted(glob.glob("xml/out_ps_*.txt")):
    model = os.path.basename(out)[7:-4]
    if not os.path.exists(f"xml/done_ps_{model}"): continue                                       # a path-sampling run that has not finished
    dirs = sorted(glob.glob(f"xml/ps_{model}/step*"), key=lambda p: int(p.split("step")[-1]))
    if not dirs: continue
    S = []
    for i, d in enumerate(dirs):
        v = pd.read_csv(os.path.join(d, "likelihood.log"), sep="\t", comment="#").likelihood.to_numpy(); v = v[int(len(v) * BURN):][::THIN]; S.append(v)
        samples += [dict(model=model, step=i, sample=j, loglik=x) for j, x in enumerate(v)]
    k = len(S); b = betas(k); E = [s.mean() for s in S]; txt = open(out).read(); m = re.search(r"marginal L estimate = (-?[\d.]+)", txt)
    ess_steps = [float(r[4]) for r in re.findall(r"^(\d+)\s+(\S+)\s+(-?\d[\d.]*)\s+(-?[\d.]+)\s+([\d.]+)\s*$", txt, flags=re.M)]      # BEAST 2's effective sample size of the likelihood at each step
    mrows.append(dict(model=model, steps=k, ps_own=ps_est(b, E), ss_own=ss_est(b, S), ss_beast=float(m.group(1)) if m else np.nan, mean_loglik_beta1=E[0], samples_per_step=len(S[0]), min_step_ess=min(ess_steps) if ess_steps else np.nan, steps_ess_below_100=sum(x < 100 for x in ess_steps), steps_with_ess=len(ess_steps)))
M = pd.DataFrame(mrows); M.round(4).to_csv("table_beast_marginal.csv", index=False)
if samples: pd.DataFrame(samples).round(4).to_csv("ps_step_samples.csv.gz", index=False)
# ---- Bayes factors (natural log of the marginal likelihood differences, stepping-stone and path sampling)
bf = []
if len(M):
    g = M.set_index("model"); pairs = [("het_strict", "iso", "real dates against no dates (isochronous)")] + [("het_strict", d, f"real dates against randomised dates {d[2:]}") for d in ("dr1", "dr2", "dr3")] + [("het_strict", "het_relaxed", "strict against relaxed clock"), ("het_relaxed", "het_strict", "relaxed against strict clock")]
    for a, b_, lab in pairs:
        if a in g.index and b_ in g.index: bf.append(dict(comparison=lab, model_a=a, model_b=b_, ln_bf_ss=g.ss_own[a] - g.ss_own[b_], ln_bf_ps=g.ps_own[a] - g.ps_own[b_]))
B = pd.DataFrame(bf); B.round(3).to_csv("table_beast_bf.csv", index=False)
# ---- TreeTime on the same genomes (the 60 taxa of the XML), for comparison with BEAST
names = re.findall(r'<sequence [^>]*?taxon="([^"]+)"', open("xml/mcmc_het_strict.xml").read())
if names:
    from Bio import Phylo, AlignIO, SeqIO
    tr = Phylo.read("ebola_tree.nwk", "newick"); keep = set(names)
    for t in [t for t in tr.get_terminals() if t.name not in keep]: tr.prune(t)
    with tempfile.TemporaryDirectory() as tmp:
        Phylo.write(tr, f"{tmp}/t.nwk", "newick"); aln = AlignIO.read("ebola_aln.fasta", "fasta"); SeqIO.write([r for r in aln if r.id in keep], f"{tmp}/a.fasta", "fasta")
        d = pd.read_csv("ebola_dates.csv"); d[d.acc.isin(keep)].to_csv(f"{tmp}/d.csv", index=False)
        subprocess.run(["treetime", "--tree", f"{tmp}/t.nwk", "--aln", f"{tmp}/a.fasta", "--dates", f"{tmp}/d.csv", "--name-column", "acc", "--date-column", "date", "--covariation", "--confidence", "--outdir", f"{tmp}/o", "--rng-seed", "1"], capture_output=True, text=True)
        mc = open(f"{tmp}/o/molecular_clock.txt").read(); rate = float(re.search(r"--rate:\s+([\d.eE+-]+)", mc).group(1)); sd = float(re.search(r"\+/-\s+([\d.eE+-]+)", mc).group(1))
        rws = [l.rstrip("\n").split("\t") for l in open(f"{tmp}/o/dates.tsv") if not l.startswith("#")]; nd = [r for r in rws if r[0].startswith("NODE_") and r[2] != "--"]; root = min(nd, key=lambda r: float(r[2]))
    pd.DataFrame([dict(taxa=len(names), clock_rate=rate, rate_sd=sd, rate_per_genome_year=rate * L, root_date=float(root[2]), root_lower=float(root[3]), root_upper=float(root[4]))]).round(8).to_csv("table_beast_treetime.csv", index=False)
# ---- figure
if len(M):
    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.4)); order = [m for m in ["het_strict", "het_relaxed", "iso", "dr1", "dr2", "dr3"] if m in set(M.model)]; lab = {"het_strict": "Real dates,\nstrict clock", "het_relaxed": "Real dates,\nrelaxed clock", "iso": "No dates", "dr1": "Randomised\ndates 1", "dr2": "Randomised\ndates 2", "dr3": "Randomised\ndates 3"}; g = M.set_index("model")
    ax[0].plot(range(len(order)), [g.ss_own[m] for m in order], "o", color="#2b6cb0", label="stepping stone"); ax[0].plot(range(len(order)), [g.ps_own[m] for m in order], "x", color="#c05621", label="path sampling"); ax[0].set(xticks=range(len(order)), xticklabels=[lab[m] for m in order], ylabel="Log marginal likelihood", title="BETS and clock models: marginal likelihoods"); ax[0].legend(frameon=False)
    if len(P) and "rootdate" in set(P.parameter) and os.path.exists("table_beast_treetime.csv"):
        t = pd.read_csv("table_beast_treetime.csv").iloc[0]; r = P[(P.model == "het_strict") & (P.parameter == "rootdate")].iloc[0]; r2 = P[(P.model == "het_relaxed") & (P.parameter == "rootdate")]
        ax[1].errorbar([r["mean"]], [2], xerr=[[r["mean"] - r.hpd_lower], [r.hpd_upper - r["mean"]]], fmt="o", color="#2b6cb0", capsize=4, label="BEAST 2, strict clock (95% HPD)")
        if len(r2): r2 = r2.iloc[0]; ax[1].errorbar([r2["mean"]], [1], xerr=[[r2["mean"] - r2.hpd_lower], [r2.hpd_upper - r2["mean"]]], fmt="s", color="#1b998b", capsize=4, label="BEAST 2, relaxed clock (95% HPD)")
        ax[1].errorbar([t.root_date], [0], xerr=[[t.root_date - t.root_lower], [t.root_upper - t.root_date]], fmt="D", color="#c05621", capsize=4, label="TreeTime (90% interval)"); ax[1].set(yticks=[], xlabel="Date of the common ancestor (year)", title="Root date of the same genomes"); ax[1].legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout(); fig.savefig("figs6/fig_ch6_beast.png", dpi=180)
print(P.round(5).to_string(index=False) if len(P) else "no posterior logs"); print(M.round(3).to_string(index=False) if len(M) else "no path sampling yet"); print(B.round(2).to_string(index=False) if len(B) else "")
