"""
Chapter 6 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts: its own FASTA reader and sample filters, ordinary least squares and a date-randomisation test on the
shipped table of Hands-on 1, a Python re-implementation of the skyline experiment of Hands-on 2 (200 genealogies, its own random numbers), the closed forms and a fresh simulation of the time-dependent model of Hands-on 4, the reproduction-number formulas, a
least-squares root-to-tip regression on the divergence trees of the TreeTime runs of Hands-on 3, and the relations among the TreeTime results (imposed rates, shuffled dates, skyline). Results are written to validation_results_ch6.csv.
USAGE   python validate_ch6.py [--chapter Chapter_6_....docx] [--doi] [--full] [--sim]
        --chapter  also checks the chapter's tables, sentences, section numbers, figure numbers and references against the results (needs pandoc and python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns Hands-on 1, 2 and 4 in a scratch folder and checks that their result files are reproduced byte for byte (about 1 minute; needs mappy, Rscript with ggplot2)
        --sim      with --full, also reruns the TreeTime analyses of Hands-on 3 from scratch (about 25 minutes; needs treetime and FastTree)
Prints PASS/FAIL/SKIP per check and exits 1 if any check fails.
"""
import argparse, datetime, glob, hashlib, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy import stats, optimize
from Bio import Phylo
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); ap.add_argument("--sim", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data", "../data/beast", "data/beast", "beast", "../data/beast/xml", "data/beast/xml"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
def readfa(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"): names.append(line[1:].split()[0]); seqs.append([])
        elif line: seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]
# ================================================================== A. Hands-on 1: temporal signal in 327 SARS-CoV-2 genomes
nm, sq = readfa(P("sars2.fasta")); seq = dict(zip(nm, sq)); ref = seq["NC_045512.2"]; L = len(ref); meta = pd.read_csv(P("sars2_meta.csv")); tab = pd.read_csv(P("sars2_root_to_tip_table.csv"))
day = meta.date.astype(str).str.match(r"^(\d{4}-\d{2}-\d{2}|\d{2}-[A-Za-z]{3}-\d{4})$"); dt = pd.to_datetime(meta.date.where(day), errors="coerce", format="mixed"); keep = meta[day & dt.notna() & (dt >= "2019-12-01") & (dt <= "2021-03-31") & meta.acc.isin(seq) & (meta.acc != "NC_045512.2")].copy()
keep = keep[[sum(c not in "ACGT" for c in seq[a]) / len(seq[a]) < 0.02 for a in keep.acc]].drop_duplicates("acc")
check("A.1", "sample filters recomputed from the FASTA and the metadata: 431 genomes sampled, 395 with a day-level date, 327 kept; the kept accessions equal those of the result table; reference length 29,903", len(meta) == 431 and int(day.sum()) == 395 and len(keep) == 327 and set(keep.acc) == set(tab.acc) and L == 29903, (len(meta), int(day.sum()), len(keep), L))
t = ((pd.to_datetime(tab.date) - pd.Timestamp("2019-12-26")).dt.days / 365.25).to_numpy(); y = tab.distance_per_site.to_numpy(); n = len(t); sx, sy = t - t.mean(), y - y.mean(); slope = (sx * sy).sum() / (sx ** 2).sum(); icpt = y.mean() - slope * t.mean(); res = y - icpt - slope * t; se = math.sqrt((res ** 2).sum() / (n - 2) / (sx ** 2).sum()); r2 = 1 - (res ** 2).sum() / (sy ** 2).sum()
xint = pd.Timestamp("2019-12-26") + pd.Timedelta(days=-icpt / slope * 365.25)
check("A.2", "ordinary least squares on the table: slope 6.59 × 10⁻⁴ per site per year (SE 1.4 × 10⁻⁵), 19.7 per genome-year, about 1.6 per month, R² 0.87, x-intercept 17 November 2019", abs(slope - 6.59e-4) < 6e-7 and abs(se - 1.4e-5) < 5e-7 and rnd(slope * L, 1) == "19.7" and rnd(slope * L / 12, 1) == "1.6" and rnd(r2, 2) == "0.87" and xint.strftime("%d %b %Y") == "17 Nov 2019", (slope, se, round(slope * L, 2), round(r2, 3), xint.date()))
rng = np.random.default_rng(5); perm = np.array([((t[rng.permutation(n)] - t.mean()) * sy).sum() / (sx ** 2).sum() for _ in range(1000)]); pv = (np.sum(perm >= slope) + 1) / 1001
check("A.3", "own date randomisation (1,000 shuffles): the observed slope exceeds every shuffled slope, which lie within ±1.2 × 10⁻⁴ (p = 0.001)", (perm < slope).all() and abs(perm.min() + 1.2e-4) < 4e-5 and abs(perm.max() - 1.2e-4) < 4e-5 and rnd(pv, 3) == "0.001", (round(perm.min(), 6), round(perm.max(), 6)))
check("A.4", "the Jukes–Cantor distance per site in the table is the substitution count over the genome length to within 5%, genome by genome (a check on the alignment step)", np.all(np.abs(tab.distance_per_site * L - tab.substitutions) <= 0.05 * np.maximum(tab.substitutions, 1) + 1.5), float(np.max(np.abs(tab.distance_per_site * L - tab.substitutions))))
# ================================================================== B. Hands-on 2: skyline experiment, own implementation and random numbers
def sim(n, N0, r, rg):
    t = 0.0; ct = np.empty(n - 1)
    for k in range(n, 1, -1): t = math.log(math.exp(r * t) - r * math.log(rg.random()) / (k * (k - 1) / (2 * N0))) / r; ct[n - k] = t
    return ct
def estimators(ct, n, m=10):
    k = np.arange(n, 1, -1); g = np.diff(np.r_[0, ct]); kk = k * (k - 1) / 2; grp = np.ceil(np.arange(1, len(k) + 1) / m).astype(int); tt, nn = [], []
    for gi in np.unique(grp):
        i = np.where(grp == gi)[0]; tt.append(np.r_[0, ct][[i.min(), i.max() + 1]].mean()); nn.append((g[i] * kk[i]).sum() / len(i))
    sk = -np.polyfit(tt, np.log(nn), 1)[0]; tprev = np.r_[0, ct[:-1]]
    def nll(par): lN0, rr = par[0], math.exp(par[1]); return -(np.sum(np.log(kk) + rr * ct - lN0) - np.sum(kk * (np.exp(rr * ct) - np.exp(rr * tprev)) / (rr * math.exp(lN0))))
    o = optimize.minimize(nll, [math.log(500), math.log(0.3)], method="Nelder-Mead"); return sk, math.exp(o.x[1]), math.exp(o.x[0])
rg = np.random.default_rng(2026); est = np.array([estimators(sim(100, 1000, 0.5, rg), 100)[:2] for _ in range(200)])
check("B.1", "own simulation of 200 genealogies of 100 lineages (N0 = 1,000, r = 0.5): the skyline regression averages 0.46 (SD 0.04) and the exponential maximum-likelihood fit 0.51 (SD 0.04), as in the chapter (means within 0.03, SDs between 0.025 and 0.06)", abs(est[:, 0].mean() - 0.46) < 0.03 and abs(est[:, 1].mean() - 0.51) < 0.03 and 0.025 < est[:, 0].std(ddof=1) < 0.06 and 0.025 < est[:, 1].std(ddof=1) < 0.06 and est[:, 0].mean() < est[:, 1].mean(), (round(est[:, 0].mean(), 3), round(est[:, 0].std(ddof=1), 3), round(est[:, 1].mean(), 3), round(est[:, 1].std(ddof=1), 3)))
out2 = open(P("ex6_skyline_output.txt")).read()
check("B.2", "output of ex6_skyline.R for the one genealogy: rate 0.39 (0.31 to 0.48), maximum likelihood 0.42 with N0 = 686, root at 12.94 years, 200 genealogies 0.46 (0.04) and 0.51 (0.04)", all(x in out2 for x in ("root at 12.94", "0.39 per year (95% CI 0.31 to 0.48)", "N0 = 686, r = 0.42", "skyline regression 0.46 (SD 0.04); exponential ML 0.51 (SD 0.04)")))
# ================================================================== C. Hands-on 4: the time-dependent model and the reproduction numbers
MU_N, MU_D, LAM = 10.0, 40.0, 0.02; Nf = lambda tau: MU_N * tau + MU_D * (1 - math.exp(-LAM * tau)) / LAM
td = pd.read_csv(P("table_time_dependent.csv")); dd = pd.read_csv(P("table_time_dependent_dating.csv")); gr = pd.read_csv(P("table_growth_R.csv"))
check("C.1", "closed-form apparent rate N(τ)/τ of the model recomputed at 1, 10, 100, 1,000 and 10,000 years: 49.6, 46.3, 27.3, 12.0 and 10.2 per genome-year, falling from μn + μd = 50 towards μn = 10", all(abs(Nf(a) / a - c) < 5e-4 for a, c in zip(td.tau_years, td.closed_form_rate)) and [rnd(Nf(a) / a, 1) for a in (1, 10, 100, 1000, 10000)] == ["49.6", "46.3", "27.3", "12.0", "10.2"])
rg2 = np.random.default_rng(99); okm = True; sims = []
for tau in (1, 10, 100):
    tot = np.empty(600)
    for i in range(600): k = rg2.poisson(MU_D * tau); u = rg2.uniform(0, tau, k); life = rg2.exponential(1 / LAM, k); tot[i] = rg2.poisson(MU_N * tau) + int((u + life > tau).sum())
    sims.append(tot.mean() / tau); okm &= abs(tot.mean() / tau - Nf(tau) / tau) < 4 * tot.std() / tau / math.sqrt(600) + 0.01
check("C.2", "a fresh simulation of every mutation (600 lineages at 1, 10 and 100 years, different random numbers) agrees with the closed form within four standard errors", okm, [round(x, 2) for x in sims])
r10 = Nf(10) / 10; check("C.3", "dating with the 10-year rate: a 100-year divergence is dated at 59 years (59%), a 1,000-year divergence at 259 years (26%), a 10,000-year divergence at 2,205 years (22%); the rate at 10 years is 4.6 times neutral", [round(Nf(a) / r10) for a in (100, 1000, 10000)] == [59, 259, 2205] and [round(100 * Nf(a) / r10 / a) for a in (100, 1000, 10000)] == [59, 26, 22] and rnd(r10 / MU_N, 1) == "4.6" and all(abs(Nf(a) / r10 - e) < 5e-3 for a, e in zip(dd.true_age_years, dd.estimated_age_years)))
okg = True
for r_ in gr.itertuples():
    rr = math.log(2) / r_.doubling_time_days; th = r_.mean_generation_interval_days / 2; M = lambda s_: (1 - s_ * th) ** (-2); okg &= abs(1 / M(-rr) - r_.R) < 6e-5 and abs(1 + rr * r_.mean_generation_interval_days - r_.R_exponential_interval) < 6e-5 and abs(rr - r_.growth_rate_per_day) < 6e-5
check("C.4", "reproduction numbers recomputed from R = 1/M(−r) with the gamma moment generating function (shape 2) and from R = 1 + rT for an exponential interval, for three doubling times and three mean intervals; a 14-day doubling gives 1.56, 1.88 and 2.24", okg and [rnd(gr[(gr.doubling_time_days == 14) & (gr.mean_generation_interval_days == m_)].R.iloc[0], 2) for m_ in (10, 15, 20)] == ["1.56", "1.88", "2.24"])
check("C.5", "rate arithmetic of the text: 1 to 2 × 10⁻⁸ per generation at 25 years is 4 to 8 × 10⁻¹⁰ per year; Y-STR pedigree rates 2 to 3 × 10⁻³ against 6.9 × 10⁻⁴ are 2.9 to 4.3 times higher; 10⁻⁷ to 10⁻⁶ per site per year over 4.4 Mb is 0.4 to 4.4 substitutions per genome-year; RNA virus 10⁻³ to 10⁻⁴ against M. tuberculosis 10⁻⁷ to 10⁻⁶ spans two to four orders of magnitude", abs(1e-8 / 25 - 4e-10) < 1e-15 and abs(2e-8 / 25 - 8e-10) < 1e-15 and rnd(2e-3 / 6.9e-4, 1) == "2.9" and rnd(3e-3 / 6.9e-4, 1) == "4.3" and rnd(1e-7 * 4.4e6, 1) == "0.4" and rnd(1e-6 * 4.4e6, 1) == "4.4" and round(math.log10(1e-4 / 1e-6)) == 2 and round(math.log10(1e-3 / 1e-7)) == 4)
# ================================================================== D. Hands-on 3: TreeTime results
T = pd.read_csv(P("table_treetime.csv")); DR = pd.read_csv(P("table_treetime_drt.csv")); SK = pd.read_csv(P("table_skyline_ebola.csv")); LEN = {"Ebola": 18961, "SARS-CoV-2": 29903}
def ols_tree(treefile, datefile, drop=()):
    tr = Phylo.read(P(treefile), "nexus"); dates = dict(pd.read_csv(P(datefile)).pipe(lambda d: zip(d.acc, [pd.Timestamp(x).year + (pd.Timestamp(x).dayofyear - 0.5) / 365.25 for x in d.date]))); xs, ys = [], []
    for tip in tr.get_terminals():
        if tip.name in drop or tip.name not in dates: continue
        xs.append(dates[tip.name]); ys.append(tr.distance(tr.root, tip))
    xs, ys = np.array(xs), np.array(ys); b = np.polyfit(xs, ys, 1); return b[0], float(np.corrcoef(xs, ys)[0, 1] ** 2), len(xs)
okD = True; det = []
for ds, tf, df in (("Ebola", "divergence_tree_ebola.nexus", "ebola_dates.csv"), ("SARS-CoV-2", "divergence_tree_sars.nexus", "sars2_dates.csv")):
    rw = T[(T.dataset == ds) & (T.scenario == "default")].iloc[0]; drop = [x for x in str(rw.outlier_ids).split(";") if x and x != "nan"]; b, r2_, n_ = ols_tree(tf, df, drop); okD &= abs(b - rw.clock_rate) / rw.clock_rate < 0.02 and abs(r2_ - rw.r2) < 0.02; det.append((n_, round(b, 6), round(r2_, 2)))
check("D.1", "own least-squares regression of root-to-tip distance on date, on the rerooted divergence tree of the default TreeTime run with its outlier removed, reproduces the rate (within 2%) and R² (within 0.02) of Ebola and SARS-CoV-2", okD, det)
okS = True
for ds in LEN:
    c = T[(T.dataset == ds) & (T.scenario == "covariation")].iloc[0]; lo, hi = T[(T.dataset == ds) & (T.scenario == "rate x0.5")].iloc[0], T[(T.dataset == ds) & (T.scenario == "rate x2")].iloc[0]
    okS &= abs(lo.clock_rate / c.clock_rate - 0.5) < 1e-3 and abs(hi.clock_rate / c.clock_rate - 2) < 1e-3 and lo.root_date < c.root_date < hi.root_date and c.root_lower < c.root_date < c.root_upper and abs(lo.rate_sd / lo.clock_rate - 0.2) < 1e-3
check("D.2", "relations among the TreeTime runs: the imposed rates are half and twice the covariation rate with a standard deviation of 20%; the root date moves later as the imposed rate rises; each covariation root date lies inside its own 90% interval", okS)
okR = True; det = []
for ds in LEN:
    d = DR[DR.dataset == ds]; o = float(d[d.replicate == 0].rate.iloc[0]); pm = d[d.replicate > 0].rate.dropna(); det.append((len(pm), len(d[d.replicate > 0]), round(o * LEN[ds], 1), round(pm.min() * LEN[ds], 1), round(pm.max() * LEN[ds], 1))); nall = len(d[d.replicate > 0]); okR &= nall == {"Ebola": 50, "SARS-CoV-2": 30}[ds] and (pm < o).all() and 0 < nall - len(pm) <= 6
check("D.3", "date randomisation with TreeTime: 50 shuffles for Ebola and 30 for SARS-CoV-2; every shuffle that could be analysed gives a rate below the observed rate, and the few that could not (no root with a positive rate) count as below it", okR, det)
rsd = lambda d, c: float(T[(T.dataset == d) & (T.scenario == "covariation")].iloc[0][c])
chk = (rnd(T[(T.dataset == "Ebola") & (T.scenario == "default")].iloc[0].clock_rate * 18961, 1), rnd(rsd("Ebola", "clock_rate") * 18961, 1), rnd(rsd("Ebola", "rate_sd") * 18961, 1), rnd(rsd("SARS-CoV-2", "clock_rate") * 29903, 1))
check("D.4", "rates per genome-year: Ebola 21.6 (least squares) and 17.6 ± 3.4 (covariation); SARS-CoV-2 16.6 with covariation", chk == ("21.6", "17.6", "3.4", "16.6"), chk)
def tsv_root(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]; n = [r for r in rows if r[0].startswith("NODE_") and r[2] != "--"]; return min(n, key=lambda r: float(r[2]))
okT = True
for ds, f in (("Ebola", "dates_ebola_covariation.tsv"), ("SARS-CoV-2", "dates_sars_covariation.tsv")):
    r_ = tsv_root(P(f)); c = T[(T.dataset == ds) & (T.scenario == "covariation")].iloc[0]; okT &= abs(float(r_[2]) - c.root_date) < 1e-6 and abs(float(r_[3]) - c.root_lower) < 1e-6 and abs(float(r_[4]) - c.root_upper) < 1e-6 and r_[1] == (datetime.date(int(c.root_date), 1, 1) + datetime.timedelta(days=(c.root_date - int(c.root_date)) * (366 if int(c.root_date) % 4 == 0 else 365))).isoformat()
check("D.5", "the root rows of TreeTime's own dates.tsv (covariation runs) equal the table, and the calendar date TreeTime prints equals the date converted from the decimal year", okT)
fd = lambda f: min(pd.to_datetime(pd.read_csv(P(f)).date)).date(); dB = (fd("ebola_dates.csv") - (datetime.date(2014, 1, 1) + datetime.timedelta(days=rsd("Ebola", "root_date") % 1 * 365))).days; dS = (fd("sars2_dates.csv") - (datetime.date(2019, 1, 1) + datetime.timedelta(days=rsd("SARS-CoV-2", "root_date") % 1 * 365))).days
check("D.6", "first collection dates: Ebola 27 May 2014 and SARS-CoV-2 23 December 2019; the covariation root dates lie 79 and 102 days earlier", str(fd("ebola_dates.csv")) == "2014-05-27" and str(fd("sars2_dates.csv")) == "2019-12-23" and abs(dB - 79) <= 1 and abs(dS - 102) <= 1, (dB, dS))
pk = SK.loc[SK.Ne.idxmax()]
check("D.7", "Ebola skyline: 8 intervals, positive sizes inside their bounds, a rise from 91 to a peak of 1,955 near 2014.2 and a fall to 82 by 2014.7", len(SK) == 8 and (SK.Ne > 0).all() and ((SK.lower <= SK.Ne) & (SK.Ne <= SK.upper)).all() and round(SK.Ne.iloc[0]) == 91 and round(pk.Ne) == 1955 and abs(pk.date - 2014.235) < 0.01 and round(SK.Ne.iloc[-1]) == 82, (round(SK.Ne.iloc[0]), round(pk.Ne), round(SK.Ne.iloc[-1])))
# ================================================================== E. the chapter text (needs --chapter)
if A.chapter:
    from docx import Document
    d = Document(A.chapter); paras = [p.text for p in d.paragraphs]; Tb = d.tables; alltext = "\n".join(paras) + "\n" + "\n".join(c.text for t in Tb for r in t.rows for c in r.cells)
    md = subprocess.run(["pandoc", A.chapter, "-t", "markdown", "--wrap=none"], capture_output=True, text=True).stdout
    heads = [(p.style.name, p.text) for p in d.paragraphs if p.style.name.startswith("Heading 2")]; nums = [int(m.group(1)) for _, h in heads for m in [re.match(r"^(\d+)\. ", h)] if m]; titles = {int(re.match(r"^(\d+)\. ", h).group(1)): h for _, h in heads if re.match(r"^(\d+)\. ", h)}
    h2n = [h for _, h in heads]; order = ["Learning objectives", "At a glance", "Introduction and core concepts", "Clock models", "Calibration", "Bayesian phylodynamic", "Time trees by maximum likelihood", "Case study 1: MRSA", "Case study 2: SARS-CoV-2 Omicron", "Case study 3: TB-MDR", "Case study 4: Zika", "Hands-on 1 (Python)", "Hands-on 2 (R)", "Hands-on 3 (Python)", "Hands-on 4 (Python)", "Hands-on 5 (BEAST 2)", "Agentic AI in phylodynamic", "Common pitfalls", "Links to other chapters", "Sources"]
    pos_ = [next((i for i, h in enumerate(h2n) if h.startswith(w)), -1) for w in order]
    check("E.1", "the 20 main sections run in the template order: objectives, at a glance, five core sections, four case studies, Hands-on 1 to 5, agent, pitfalls, links, sources", all(x_ >= 0 for x_ in pos_) and pos_ == sorted(pos_), [w for w, x_ in zip(order, pos_) if x_ < 0])
    refs = ["The Omicron and Zika case studies", "The MRSA and TB-MDR case studies are ill", "(the MRSA case study)", "as in the MRSA and TB-MDR case studies", "listed in the pitfalls section", "the Zika case study"]
    check("E.2", "the in-text cross-references name the sections and no numbered section reference remains", all(x_ in alltext for x_ in refs) and re.search(r"\bSections? \d", alltext) is None, [x_ for x_ in refs if x_ not in alltext])
    caps = [int(m) for m in re.findall(r"^Figure 6\.(\d+)\. ", alltext, flags=re.M)]; check("E.3", "figure captions are numbered 6.1 to 6.7 in order, one for every picture (6.4 TreeTime, 6.5 time-dependent rates, 6.6 BEAST 2, 6.7 the agent workflow); all 7 pictures have alt text", caps == list(range(1, 8)) and len(re.findall(r'<wp:docPr [^>]*descr="[^"]{10,}"', __import__("zipfile").ZipFile(A.chapter).read("word/document.xml").decode("utf8", "ignore"))) == 7, caps)
    objs = (lambda a, b: [p for p in paras[a + 1:b] if len(p) > 60])(paras.index("Learning objectives"), paras.index("About the worked examples")); check("E.4", "ten learning objectives; no mention of uploaded files or drafts; the Soares and Rito review is cited and listed", len(objs) == 10 and re.search(r"uploaded (draft|file|files|script|document)", alltext.lower()) is None and "Soares and Rito, 2026" in alltext, len(objs))
    # ---- Hands-on 3 table, rebuilt from the result table with dates converted independently of the chapter patch
    def dstr(x): y = int(x); dd = datetime.date(y, 1, 1) + datetime.timedelta(days=(x - y) * (366 if y % 4 == 0 else 365)); return f"{dd.day} {dd.strftime('%b')} {dd.year}"
    def dint(lo, hi): a, b = dstr(lo), dstr(hi); return (a.rsplit(" ", 1)[0] if a.split()[-1] == b.split()[-1] else a) + " to " + b
    def cell(r, ds): g = r.clock_rate * LEN[ds]; return f"{g:.1f}" + (f" ± {r.rate_sd * LEN[ds]:.1f}" if not np.isnan(r.rate_sd) else "")
    rt = lambda r: dstr(r.root_date) + ("" if np.isnan(r.root_lower) else f" ({dint(r.root_lower, r.root_upper)})")
    tt = [t for t in Tb if t.rows[0].cells[0].text == "Analysis"][0]; bad = []
    for i, sc in enumerate(("default", "covariation", "nofilter", "rate x0.5", "rate x2")):
        cs = [c.text.strip() for c in tt.rows[i + 1].cells]; e, s_ = T[(T.dataset == "Ebola") & (T.scenario == sc)].iloc[0], T[(T.dataset == "SARS-CoV-2") & (T.scenario == sc)].iloc[0]
        exp = [cell(e, "Ebola"), f"{e.r2:.2f}", rt(e), cell(s_, "SARS-CoV-2"), f"{s_.r2:.2f}", rt(s_)]
        if cs[1:] != exp: bad.append((sc, cs[1:], exp))
    check("E.5", "the Hands-on 3 table (rate per genome-year, R² and root date with its 90% interval, five analyses, two data sets) equals table_treetime.csv", not bad and len(tt.rows) == 6, bad[:1])
    def days_before(x, first): y = int(x); dd = datetime.date(y, 1, 1) + datetime.timedelta(days=(x - y) * (366 if y % 4 == 0 else 365)); return (first - dd).days
    ce, cs_ = T[(T.dataset == "Ebola") & (T.scenario == "covariation")].iloc[0], T[(T.dataset == "SARS-CoV-2") & (T.scenario == "covariation")].iloc[0]; g = lambda ds, sc: T[(T.dataset == ds) & (T.scenario == sc)].iloc[0]
    sent = [f"{days_before(ce.root_date, datetime.date(2014, 5, 27))} days before the first sample", f"{days_before(cs_.root_date, datetime.date(2019, 12, 23))} days before the first sample", f"a spread of {abs(g('Ebola', 'rate x2').root_date - g('Ebola', 'rate x0.5').root_date) * 365:.0f} days", f"{abs(g('SARS-CoV-2', 'rate x2').root_date - g('SARS-CoV-2', 'rate x0.5').root_date) * 365:.0f} days.", "the genome whose private substitutions cluster in one window in Chapter 5", "KR105312.1, the genome that Chapter 5's outlier rule flagged", f"falls from {g('Ebola', 'covariation').r2:.2f} to {g('Ebola', 'nofilter').r2:.2f}", "1.25 × 10⁻³ (23.7 per genome-year, R² 0.58)", f"rises from {SK.Ne.iloc[0]:.0f} in late 2013 to a peak of {SK.Ne.max():,.0f} in early 2014 and falls to {SK.Ne.iloc[-1]:.0f} by late 2014", "one-sided p = 0.020", "p = 0.032", "5 of 50 had no root with a positive rate"]
    miss = [x for x in sent if x not in alltext.replace("\u2019", "'")]; c5 = pd.read_csv(P("ch5_table_root_to_tip.csv")); b5 = c5[(c5.dataset.str.startswith("Ebola")) & (c5.rooting == "best")].iloc[0]
    check("E.6", "sentences of Hands-on 3 that quote results (days before the first sample, the spread under imposed rates, outlier identities, R², the skyline, the shuffle p-values) match the tables; Chapter 5's best-root Ebola regression is 1.25 × 10⁻³ (23.7, R² 0.58)", not miss and rnd(b5.slope_per_site_per_year * 1000, 2) == "1.25" and rnd(b5.substitutions_per_genome_year, 1) == "23.7" and rnd(b5.r2, 2) == "0.58", miss[:2])
    ok4 = []; tab_ = lambda head: [t for t in Tb if t.rows[0].cells[0].text == head][0]
    t1 = tab_("Age of the lineage (years)"); ok4.append(all([c.text.strip() for c in t1.rows[i + 1].cells] == [f"{int(r_.tau_years):,}", f"{r_.closed_form_rate:.1f}", (f"{r_.simulated_rate:.1f} ({r_.simulated_sd:.1f})" if not np.isnan(r_.simulated_rate) else "not simulated"), f"{r_.ratio_to_neutral:.2f}"] for i, r_ in enumerate(td.itertuples())))
    t2 = tab_("True age (years)"); ok4.append(all([c.text.strip() for c in t2.rows[i + 1].cells] == [f"{int(r_.true_age_years):,}", f"{r_.expected_substitutions:,.0f}", f"{r_.estimated_age_years:,.0f}", f"{100 * r_.ratio_estimated_to_true:.0f}%"] for i, r_ in enumerate(dd.itertuples())))
    t3 = tab_("Doubling time (days)"); ok4.append(all([c.text.strip() for c in t3.rows[i + 1].cells] == [str(int(r_.doubling_time_days)), f"{r_.growth_rate_per_day:.3f}", str(int(r_.mean_generation_interval_days)), f"{r_.R:.2f}", f"{r_.R_exponential_interval:.2f}"] for i, r_ in enumerate(gr.itertuples())))
    check("E.7", "the Hands-on 4 tables (apparent rate by age, dating with the 10-year rate, reproduction numbers by doubling time and generation interval) equal the result tables", all(ok4), ok4)
    d_in, d_ref = (lambda t_: set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", t_)))(md.split("## Sources")[0]), (lambda t_: set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", t_)))(md.split("## Sources")[1])
    check("E.8", "every DOI cited in the text is in the reference list, which holds 45 DOIs", d_in <= d_ref and len(d_ref) == 45, (len(d_in), len(d_ref), sorted(d_in - d_ref)[:3]))
    if A.doi:
        import urllib.request, urllib.parse, time
        badd = []
        for dd_ in sorted(d_ref | d_in):
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(dd_, safe="/()"), headers={"User-Agent": "validate-ch6"}), timeout=30).read()
            except Exception as ex: badd.append((dd_, str(ex)[:30]))
            time.sleep(0.15)
        check("DOI", f"all {len(d_ref | d_in)} DOIs of the chapter resolve in Crossref", not badd, badd[:4])
else: skip("E.*", "chapter text checks", "no --chapter given")
# ================================================================== F. reruns (needs --full)
if A.full:
    tmp = tempfile.mkdtemp(); src = os.path.dirname(os.path.abspath(P("sars2.fasta"))); code = os.path.dirname(os.path.abspath(__file__)); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) + glob.glob(os.path.join(code, "*.py")) + glob.glob(os.path.join(code, "*.R")) if os.path.isfile(f)]
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); names = ["sars2_root_to_tip_table.csv", "table_time_dependent.csv", "table_time_dependent_dating.csv", "table_growth_R.csv"]; before = {n_: md5(os.path.join(tmp, n_)) for n_ in names}; [os.remove(os.path.join(tmp, n_)) for n_ in names]
    subprocess.run([sys.executable, "ex6_root_to_tip.py"], cwd=tmp, capture_output=True); subprocess.run([sys.executable, "ex6_time_dependent.py"], cwd=tmp, capture_output=True); r2o = subprocess.run(["Rscript", "ex6_skyline.R"], cwd=tmp, capture_output=True, text=True).stdout
    for n_ in names: check("REP", f"{n_} reproduced byte for byte", os.path.exists(os.path.join(tmp, n_)) and md5(os.path.join(tmp, n_)) == before[n_])
    check("REP", "ex6_skyline.R prints the numbers of ex6_skyline_output.txt", r2o.strip() == open(os.path.join(tmp, "ex6_skyline_output.txt")).read().strip())
    if A.sim:
        env = dict(os.environ, RERUN="1"); shutil.copy(os.path.join(tmp, "table_treetime.csv"), os.path.join(tmp, "old_table_treetime.csv")); subprocess.run([sys.executable, "ex6_treetime.py"], cwd=tmp, capture_output=True, env=env)
        a_, b_ = pd.read_csv(os.path.join(tmp, "old_table_treetime.csv")), pd.read_csv(os.path.join(tmp, "table_treetime.csv")); num = ["clock_rate", "r2", "root_date", "root_lower", "root_upper", "outliers"]
        check("REP", "the TreeTime runs rerun from scratch (RERUN=1) reproduce table_treetime.csv to six figures", len(a_) == len(b_) and all(np.allclose(a_[c].fillna(0), b_[c].fillna(0), rtol=1e-5, atol=1e-6) for c in num))
    shutil.rmtree(tmp, ignore_errors=True)
# ================================================================== G. Hands-on 5: BEAST 2 (BETS, Bayes factors, posteriors)
bp = pd.read_csv(P("table_beast_posterior.csv")); bm = pd.read_csv(P("table_beast_marginal.csv")).set_index("model"); bb = pd.read_csv(P("table_beast_bf.csv")); btt = pd.read_csv(P("table_beast_treetime.csv")).iloc[0]
decy = lambda d: pd.Timestamp(d).year + (pd.Timestamp(d).dayofyear - 0.5) / (366 if pd.Timestamp(d).is_leap_year else 365)
ed = pd.read_csv(P("ebola_dates.csv")); ed["dec"] = ed.date.map(decy); srt = ed.sort_values(["date", "acc"]).reset_index(drop=True); sub60 = [srt.acc[int(round(k))] for k in np.linspace(0, len(srt) - 1, 60)]; real60 = dict(zip(ed.acc, ed.dec))
def xml_info(f):
    t = open(P("xml/" + f), encoding="utf8").read(); seqs = re.findall(r'<sequence [^>]*?taxon="([^"]+)"', t); tr = re.search(r'traitname="date" value="([^"]+)"', t).group(1); dts = {a.split("=")[0]: float(a.split("=")[1]) for a in tr.split(",")}
    prior = {k: re.search(k + r'" name="distr"([^/]*)/>', t).group(1).strip() for k in ("LogNormalClock", "LogNormalKappa", "LogNormalPop", "LaplaceGrowth", "ExpShape")}; return t, seqs, dts, prior
models = ["het_strict", "het_relaxed", "iso", "dr1", "dr2", "dr3"]; okx, det = True, []
infos = {(m, k): xml_info(f"{k}_{m}.xml") for m in models for k in ("mcmc", "ps")}
base_prior = infos[("het_strict", "mcmc")][3]
for (m, k), (t, seqs, dts, prior) in infos.items():
    ok = len(seqs) == 60 and sorted(seqs) == sorted(sub60) and prior == base_prior
    if m.startswith("het"): ok &= all(abs(dts[a] - real60[a]) < 2e-6 for a in seqs)
    elif m == "iso": ok &= len(set(round(v, 6) for v in dts.values())) == 1 and abs(max(dts.values()) - max(real60[a] for a in seqs)) < 2e-6
    else: ok &= np.allclose(sorted(dts.values()), sorted(real60[a] for a in seqs), atol=2e-6) and any(abs(dts[a] - real60[a]) > 1e-4 for a in seqs)
    ok &= (f'chainLength="6000000"' in t) if k == "mcmc" else (all(x in t for x in ('nrOfSteps="20"', 'chainLength="250000"', 'alpha="0.3"', 'preBurnin="100000"', 'burnInPercentage="20"')))
    okx &= ok; det.append((m, k, ok)) if not ok else None
check("G.1", "the 12 BEAST 2 files: the 60 genomes are every third of the 195 by collection date (recomputed from the dates file); the real-date models carry the real dates, the isochronous model one date, each randomised model a permutation of the real dates; the priors are identical in all 12 files; chain settings are 6,000,000 steps and 20 path steps of 250,000", okx, det[:3])
def ess_sokal(x):
    x = np.asarray(x, float); n = len(x); x = x - x.mean(); f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n]; ac = ac / ac[0]; tau = 1.0
    for m in range(1, n):
        tau = 1 + 2 * ac[1:m + 1].sum()
        if m >= 5 * tau: break
    return n / max(tau, 1.0)
lastd = max(real60[a] for a in sub60); logs = {"het_strict": [P("mcmc_het_strict_s1.log"), P("mcmc_het_strict_s2.log")], "het_relaxed": [P("mcmc_het_relaxed.log")]}; okp, okr, oke, detp = True, True, True, []
for m, fs in logs.items():
    ch = [pd.read_csv(f, sep="\t", comment="#") for f in fs]; ch = [c.iloc[int(len(c) * 0.1):] for c in ch]; rate = "clockRate" if m == "het_strict" else "ucldMean"
    okr &= len(ch) < 2 or bool((ch[0]["posterior"].to_numpy()[:200] != ch[1]["posterior"].to_numpy()[:200]).any())      # two different runs, not two copies of one log
    for c in ch: c["rootdate"] = lastd - c["TreeHeight"]; c["rate_genome"] = c[rate] * 18961; c["doubling_days"] = 365.25 * math.log(2) / c["growthRate"]
    for par in ("rate_genome", "rootdate", "growthRate", "doubling_days", "kappa", "TreeHeight"):
        v = np.concatenate([c[par].to_numpy() for c in ch]); row = bp[(bp.model == m) & (bp.parameter == par)].iloc[0]; xs = np.sort(v); k = int(0.95 * len(xs)); w = xs[k:] - xs[:len(xs) - k]; i0 = int(np.argmin(w))
        okp &= abs(v.mean() - row["mean"]) < 1e-5 * max(1, abs(row["mean"])) and abs(xs[i0] - row.hpd_lower) < 1e-3 * max(1, abs(row.hpd_lower) * 1e-3) + 1e-4 and abs(xs[i0 + k] - row.hpd_upper) < 1e-3 * max(1, abs(row.hpd_upper) * 1e-3) + 1e-4
        e_own = [ess_sokal(c[par].to_numpy()) for c in ch]; oke &= all(0.6 < a / b < 1.7 for a, b in zip(e_own, [row.ess_total / len(ch)] * len(ch))) and min(e_own) >= 80
        if len(ch) == 2: cc = np.array([c[par].to_numpy()[:min(len(ch[0]), len(ch[1]))] for c in ch]); nn = cc.shape[1]; W = cc.var(axis=1, ddof=1).mean(); Bv = nn * cc.mean(axis=1).var(ddof=1); rh = math.sqrt(((nn - 1) / nn * W + Bv / nn) / W); okr &= rh < 1.05 and abs(rh - row.rhat) < 5e-3
check("G.2", "posterior means and 95% highest posterior density intervals of rate, root date, growth rate, doubling time, kappa and tree height recomputed from the three BEAST 2 logs (10% burn-in) equal table_beast_posterior.csv; the root date is the last tip date minus the tree height", okp)
check("G.3", "effective sample sizes recomputed with a different estimator (Sokal's automatic windowing) agree with the table within a factor of 1.7 and are at least 80 for those parameters in every chain", oke)
check("G.4", "the two strict-clock logs are different runs (not copies of one chain) and agree: the Gelman–Rubin statistic recomputed from the logs is below 1.05 for every parameter and equals the table", okr)
st = pd.read_csv(P("ps_step_samples.csv.gz")); okm = True; detm = []
for m in models:
    d = st[st.model == m]; k = d.step.nunique(); b = np.array([((k - 1 - i) / (k - 1)) ** (1 / 0.3) for i in range(k)]); S = [d[d.step == i].loglik.to_numpy() for i in range(k)]
    ps = sum((b[i] - b[i + 1]) * (S[i].mean() + S[i + 1].mean()) / 2 for i in range(k - 1)); ss = 0.0
    for i in range(k - 1): x = (b[i] - b[i + 1]) * S[i + 1]; mx = x.max(); ss += mx + math.log(np.mean(np.exp(x - mx)))
    r = bm.loc[m]; okm &= k == 20 and abs(ps - r.ps_own) < 0.01 and abs(ss - r.ss_own) < 0.01 and abs(r.ss_beast - r.ss_own) < 0.01 and abs(r.ss_own - r.ps_own) < 5; detm.append((m, round(ps, 1), round(ss, 1), r.ss_beast))
okess = True
for m in models:
    tx = open(P(f"xml/out_ps_{m}.txt")).read(); es = [float(r[4]) for r in re.findall(r"^(\d+)\s+(\S+)\s+(-?\d[\d.]*)\s+(-?[\d.]+)\s+([\d.]+)\s*$", tx, flags=re.M)]; okess &= len(es) == 20 and abs(min(es) - bm.loc[m].min_step_ess) < 1e-6 and sum(x < 100 for x in es) == bm.loc[m].steps_ess_below_100
check("G.5b", "the effective sample size of the likelihood at each of the 20 steps, read from BEAST 2's own output of every path-sampling run, gives the smallest value and the number of steps below 100 in table_beast_marginal.csv", okess)
check("G.5", "log marginal likelihoods recomputed from the step samples of the 6 path-sampling runs (20 steps, Beta(0.3, 1) heating values, trapezoid rule for path sampling and the stepping-stone formula) equal table_beast_marginal.csv; the estimate BEAST 2 itself prints at the end of each run equals the recomputed stepping-stone estimate", okm, detm[:2])
okb = True
for r in bb.itertuples(): okb &= abs(r.ln_bf_ss - (bm.loc[r.model_a].ss_own - bm.loc[r.model_b].ss_own)) < 0.01 and abs(r.ln_bf_ps - (bm.loc[r.model_a].ps_own - bm.loc[r.model_b].ps_own)) < 0.01
check("G.6", "the log Bayes factors of table_beast_bf.csv are the differences of the log marginal likelihoods; the real-date model is better than the isochronous and each randomised model by more than 5 (stepping stone and path sampling)", okb and all(bb[bb.model_b.isin(["iso", "dr1", "dr2", "dr3"])].ln_bf_ss > 5) and all(bb[bb.model_b.isin(["iso", "dr1", "dr2", "dr3"])].ln_bf_ps > 5), bb.round(1).ln_bf_ss.tolist())
rs_ = lambda m_: dstr_x(bp[(bp.model == m_) & (bp.parameter == "rootdate")].iloc[0]["mean"])
def dstr_x(x): y = int(x); dd = datetime.date(y, 1, 1) + datetime.timedelta(days=(x - y) * (366 if y % 4 == 0 else 365)); return f"{dd.day} {dd.strftime(chr(37) + chr(98))} {dd.year}"
check("G.7", "the posterior mean root dates under the strict and the relaxed clock are within two days of each other, and the BEAST 2 clock rate (27.8 per genome-year) is higher than TreeTime's on all 195 genomes (17.6) and on the same 60 (10.8)", abs(bp[(bp.model == "het_strict") & (bp.parameter == "rootdate")].iloc[0]["mean"] - bp[(bp.model == "het_relaxed") & (bp.parameter == "rootdate")].iloc[0]["mean"]) * 365 <= 2 and bp[(bp.model == "het_strict") & (bp.parameter == "rate_genome")].iloc[0]["mean"] > T[(T.dataset == "Ebola") & (T.scenario == "covariation")].iloc[0].clock_rate * 18961 > btt.rate_per_genome_year)
if A.full:
    from Bio import Phylo, SeqIO
    tmpd = tempfile.mkdtemp(); tr = Phylo.read(P("ebola_tree.nwk"), "newick")
    for t_ in [t_ for t_ in tr.get_terminals() if t_.name not in set(sub60)]: tr.prune(t_)
    Phylo.write(tr, f"{tmpd}/t.nwk", "newick"); SeqIO.write([r for r in SeqIO.parse(P("ebola_aln.fasta"), "fasta") if r.id in set(sub60)], f"{tmpd}/a.fasta", "fasta"); ed[ed.acc.isin(sub60)][["acc", "date"]].to_csv(f"{tmpd}/d.csv", index=False)
    subprocess.run(["treetime", "--tree", f"{tmpd}/t.nwk", "--aln", f"{tmpd}/a.fasta", "--dates", f"{tmpd}/d.csv", "--name-column", "acc", "--date-column", "date", "--covariation", "--confidence", "--outdir", f"{tmpd}/o", "--rng-seed", "1"], capture_output=True, text=True)
    mc = open(f"{tmpd}/o/molecular_clock.txt").read(); rt_ = float(re.search(r"--rate:\s+([\d.eE+-]+)", mc).group(1)); nd_ = [l.rstrip("\n").split("\t") for l in open(f"{tmpd}/o/dates.tsv") if not l.startswith("#")]; rt0 = min([r_ for r_ in nd_ if r_[0].startswith("NODE_") and r_[2] != "--"], key=lambda r_: float(r_[2]))
    check("REP", "TreeTime on the same 60 genomes (pruned tree, covariation, 90% intervals) reproduces table_beast_treetime.csv", abs(rt_ - btt.clock_rate) / btt.clock_rate < 1e-3 and abs(float(rt0[2]) - btt.root_date) < 2e-3 and abs(float(rt0[3]) - btt.root_lower) < 2e-3 and abs(float(rt0[4]) - btt.root_upper) < 2e-3); shutil.rmtree(tmpd, ignore_errors=True)
if A.chapter:
    def dstr_(x): y = int(x); dd = datetime.date(y, 1, 1) + datetime.timedelta(days=(x - y) * (366 if y % 4 == 0 else 365)); return f"{dd.day} {dd.strftime('%b')} {dd.year}"
    def dint_(lo, hi): a, b = dstr_(lo), dstr_(hi); return (a.rsplit(" ", 1)[0] if a.split()[-1] == b.split()[-1] else a) + " to " + b
    pr_ = lambda m, par: bp[(bp.model == m) & (bp.parameter == par)].iloc[0]; hs, hr = "het_strict", "het_relaxed"; rows5 = []
    for lab, par, f in [("Clock rate (substitutions per genome-year)", "rate_genome", 1), ("Date of the common ancestor", "rootdate", 0), ("Growth rate (per year)", "growthRate", 1), ("Doubling time (days)", "doubling_days", 1), ("HKY kappa", "kappa", 1)]:
        a_, b_ = pr_(hs, par), pr_(hr, par)
        if par == "rootdate": rows5.append([lab, f"{dstr_(a_['mean'])} ({dint_(a_.hpd_lower, a_.hpd_upper)})", f"{a_.ess_min:.0f}", f"{a_.rhat:.3f}", f"{dstr_(b_['mean'])} ({dint_(b_.hpd_lower, b_.hpd_upper)})", f"{b_.ess_min:.0f}"])
        else: rows5.append([lab, f"{a_['mean']:.{f}f} ({a_.hpd_lower:.{f}f} to {a_.hpd_upper:.{f}f})", f"{a_.ess_min:.0f}", f"{a_.rhat:.3f}", f"{b_['mean']:.{f}f} ({b_.hpd_lower:.{f}f} to {b_.hpd_upper:.{f}f})", f"{b_.ess_min:.0f}"])
    sd_ = pr_(hr, "ucldStdev"); rows5.append(["Relaxed-clock standard deviation", "not in the model", "", "", f"{sd_['mean']:.2f} ({sd_.hpd_lower:.2f} to {sd_.hpd_upper:.2f})", f"{sd_.ess_min:.0f}"])
    t5 = [t for t in Tb if t.rows[0].cells[0].text == "Quantity"][0]; got = [[c.text.strip() for c in r.cells] for r in t5.rows[1:]]
    check("E.9", "the Hands-on 5 posterior table (rate, root date, growth rate, doubling time, kappa, relaxed-clock standard deviation; means, 95% HPD, smallest ESS and R̂) equals table_beast_posterior.csv", got == rows5 and len(t5.rows) == 7, [g_ for g_, e_ in zip(got, rows5) if g_ != e_][:1])
    order5 = [("het_strict", "Real dates, strict clock"), ("het_relaxed", "Real dates, relaxed clock"), ("iso", "No dates (isochronous), strict clock"), ("dr1", "Randomised dates 1, strict clock"), ("dr2", "Randomised dates 2, strict clock"), ("dr3", "Randomised dates 3, strict clock")]
    exp_ml = [[lab, f"{bm.loc[m].ss_own:.1f}", f"{bm.loc[m].ps_own:.1f}", "reference" if m == "het_strict" else f"{(bm.loc['het_strict'].ss_own - bm.loc[m].ss_own):+.1f}"] for m, lab in order5]
    t6 = [t for t in Tb if t.rows[0].cells[0].text == "Model"][0]; got6 = [[c.text.strip().replace("−", "-") for c in r.cells] for r in t6.rows[1:]]
    check("E.10", "the Hands-on 5 marginal-likelihood table (stepping stone, path sampling, log Bayes factor against the real-date strict-clock model) equals table_beast_marginal.csv", got6 == exp_ml, [g_ for g_, e_ in zip(got6, exp_ml) if g_ != e_][:1])
    bfx = lambda mb: bb[bb.model_b == mb].iloc[0]; drs = [bfx(f"dr{k}").ln_bf_ss for k in (1, 2, 3)]; drps = [bfx(f"dr{k}").ln_bf_ps for k in (1, 2, 3)]; r_b = pr_(hs, "rate_genome"); gi_ = lambda mean, r_yr: (1 + (r_yr / 365.25) * mean / 2) ** 2; rg = pr_(hs, "growthRate")["mean"]
    logs_t = open(P("run_all.log")).read(); mins = round(sum(int(v) for v in dict(re.findall(r"end (\S+) exit 0 (\d+) s", logs_t)).values()) / 60); tt195 = T[(T.dataset == "Ebola") & (T.scenario == "covariation")].iloc[0]
    labx = {"posterior": "the log posterior", "likelihood": "the log likelihood", "gammaShape": "the gamma shape", "popSize": "the population size", "growthRate": "the growth rate", "doubling_days": "the doubling time", "ucldStdev": "the relaxed-clock standard deviation", "kappa": "kappa", "TreeHeight": "the tree height", "rootdate": "the root date", "rate_genome": "the clock rate", "clockRate": "the clock rate", "ucldMean": "the clock rate"}
    worst_ = lambda m_: (lambda d_: f"{d_.ess_min.min():.0f} ({labx[d_.loc[d_.ess_min.idxmin()].parameter]})")(bp[bp.model == m_])
    allbf = bb[bb.model_b.isin(['iso', 'dr1', 'dr2', 'dr3'])][['ln_bf_ss', 'ln_bf_ps']]
    sent5 = [f"below 100 at {int(bm.steps_ess_below_100.min())} to {int(bm.steps_ess_below_100.max())} of the 20 steps", f"down to {bm.min_step_ess.min():.0f} at the worst step", f"by a log Bayes factor of {bfx('iso').ln_bf_ss:.1f} (path sampling {bfx('iso').ln_bf_ps:.1f})", f"by {min(drs):.1f} to {max(drs):.1f} (path sampling {min(drps):.1f} to {max(drps):.1f})", f"is {bfx('het_relaxed').ln_bf_ss:.1f} by stepping stone and {bfx('het_relaxed').ln_bf_ps:.1f} by path sampling", f"The smallest effective sample size is {worst_(hs)} for the strict clock and {worst_(hr)} for the relaxed clock", f"is {bp[(bp.model == hs) & bp.rhat.notna()].rhat.max():.3f}", f"(largest difference {abs(bm.ss_beast - bm.ss_own).max():.3f})", f"by up to {abs(bm.ss_own - bm.ps_own).max():.1f} log units", f"BEAST 2 gives {r_b['mean']:.1f} substitutions per genome-year (95% HPD {r_b.hpd_lower:.1f} to {r_b.hpd_upper:.1f})", f"TreeTime gives {btt.rate_per_genome_year:.1f} ± {btt.rate_sd * 18961:.1f} and {dstr_(btt.root_date)} ({dint_(btt.root_lower, btt.root_upper)})", f"on all 195 genomes TreeTime gave {tt195.clock_rate * 18961:.1f} ± {tt195.rate_sd * 18961:.1f} and {dstr_(tt195.root_date)}", f"a reproduction number of {gi_(15, rg):.2f}", f"{gi_(10, rg):.2f} to {gi_(20, rg):.2f} for mean intervals of 10 to 20 days", f"about {mins} minutes", f"The Bayes factors of {allbf.min().min():.0f} to {allbf.max().max():.0f} for the real dates", f"The root date is {rs_('het_strict')} under the strict clock and {rs_('het_relaxed')} under the relaxed clock"]
    miss5 = [x for x in sent5 if x not in alltext.replace("\u2019", "'")]
    check("E.11", "the sentences of Hands-on 5 that quote results (Bayes factors, step and chain effective sample sizes, R̂, estimator agreement, BEAST 2 against TreeTime, reproduction numbers, running time) match the tables and the run log", not miss5, miss5[:2])
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch6.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=100)); npass = (Rd.result == "PASS").sum(); nskip = (Rd.result == "SKIP").sum(); print(f"\n{npass} of {len(Rd) - nskip} checks passed" + (f" ({nskip} skipped)" if nskip else "")); sys.exit(0 if (Rd.result == "FAIL").sum() == 0 else 1)
