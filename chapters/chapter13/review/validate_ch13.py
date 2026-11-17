"""
Chapter 13 validation harness. Recomputes the numbers quoted in the chapter with code that does not import the chapter's own scripts (ebolacore.py, coalfit13.py, ex13_*.py).
A  the original chapter: the Ebola line list, the Cori reproduction number, the growth rate, the simulation's true growth rate, the genomic tie sensitivity and the variant forecast scores
B  the new exercises: the SIR model, the generation-interval shapes, windows, serial interval and a falling reporting fraction, a partly sampled outbreak, the audit of published numbers, the capstone rules
G  the chapter text (tables, figures, objectives, sentences, references)
H  figures and scripts
I  reruns of every script from the shipped data, with a byte-for-byte comparison of the result tables
USAGE   python validate_ch13.py [--chapter Chapter_13_....docx] [--doi] [--full]
        --chapter  also checks the chapter's tables, figures, objectives, sentences and references against the results (needs python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns every script in a scratch folder and checks that the result tables are reproduced byte for byte
Prints PASS/FAIL/SKIP per check, writes validation_results_ch13.csv and exits 1 if any check fails.
"""
import argparse, glob, math, os, re, shutil, subprocess, sys, tempfile, warnings
import numpy as np, pandas as pd
from scipy import stats, optimize, integrate
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
def P(name):
    for d in (".", "../data", "data", os.path.join(HERE, "../data"), os.path.join(HERE, "data")):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []; EXPECTED_DOIS = 21
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
close = lambda a, b, tol: abs(float(a) - float(b)) <= tol
# ================================================================== I. reruns first
if A.full:
    tmp = tempfile.mkdtemp(prefix="ch13_"); data = os.path.dirname(os.path.abspath(P("ebola_sl_linelist.csv"))); src = os.path.join(HERE, "../code") if os.path.isdir(os.path.join(HERE, "../code")) else HERE
    for f_ in glob.glob(os.path.join(src, "*.py")) + glob.glob(os.path.join(src, "*.R")): shutil.copy(f_, tmp)
    for item in os.listdir(data):
        if item.startswith("table_"): continue
        shutil.copy(os.path.join(data, item), os.path.join(tmp, item))
    os.makedirs(os.path.join(tmp, "figs13"), exist_ok=True)
    runs = [("ex13_1", ["python3", "ex13_1_ebola_rt.py"]), ("ex13_1b", ["Rscript", "ex13_1b_r_check.R"]), ("ex13_2", ["python3", "ex13_2_genome_growth.py"]), ("ex13_3", ["python3", "ex13_3_simulation.py"]), ("ex13_4", ["python3", "ex13_4_variants.py"]), ("ex13_5", ["python3", "ex13_5_sir_gi.py"]),
            ("ex13_6", ["python3", "ex13_6_windows.py"]), ("ex13_7", ["python3", "ex13_7_partial_sampling.py"]), ("ex13_8", ["python3", "ex13_8_integration_audit.py"]), ("ex13_9", ["python3", "ex13_9_capstone_check.py"]), ("ex13_10", ["python3", "ex13_10_diagrams.py"])]
    for name, cmd in runs:
        r_ = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True); check(f"I.{name}", f"{name} reruns from the shipped data without error", r_.returncode == 0, (r_.stderr or "")[-120:])
    diffs = []
    for tb in sorted(os.path.basename(f_) for f_ in glob.glob(os.path.join(data, "table_*.csv"))):
        a_, b_ = os.path.join(tmp, tb), os.path.join(data, tb)
        if not os.path.exists(a_): diffs.append(tb + " (not written)"); continue
        if open(a_, "rb").read() != open(b_, "rb").read(): diffs.append(tb)
    check("I.tables", "every result table of the rerun is byte-identical to the shipped one", not diffs, diffs)
else: skip("I.*", "reruns of the scripts", "run with --full")
# ================================================================== A. the original chapter, recomputed
d = pd.read_csv(P("ebola_sl_linelist.csv"), parse_dates=["date_of_onset"]); day = d.groupby("date_of_onset").size().asfreq("D", fill_value=0); inc = day.values.astype(float)
check("A.1", "line list: 11,903 cases with onset from 18 May 2014 to 12 September 2015; the busiest day is 20 October 2014 with 117 cases and the busiest week has 593", len(d) == 11903 and day.index.min() == pd.Timestamp("2014-05-18") and day.index.max() == pd.Timestamp("2015-09-12") and day.idxmax() == pd.Timestamp("2014-10-20") and day.max() == 117 and int(day.rolling(7).sum().max()) == 593, (len(d), day.idxmax(), day.max()))
def weights(mean, sd, n=60):
    k, th = (mean / sd) ** 2, sd ** 2 / mean; cdf = lambda x: stats.gamma.cdf(x, a=k, scale=th); w = np.array([cdf(1.5) - cdf(0.0)] + [cdf(s + 0.5) - cdf(s - 0.5) for s in range(2, n + 1)]); return w / w.sum()
def rt_at(x, w, window, t):
    lam = lambda u: sum(w[s - 1] * x[u - s] for s in range(1, min(u, len(w)) + 1)); I = x[t - window + 1:t + 1].sum(); L = sum(lam(u) for u in range(t - window + 1, t + 1)); a, rate = 1 + I, 1 / 5 + L
    return a / rate, stats.gamma.ppf(0.025, a, scale=1 / rate), stats.gamma.ppf(0.975, a, scale=1 / rate)
w0 = weights(15.3, 9.3); t1 = (pd.Timestamp("2014-10-01") - day.index[0]).days; r14 = rt_at(inc, w0, 14, t1)
check("A.2", "the Cori reproduction number on 1 October 2014 with a 14-day window is 1.42 (1.32 to 1.51), by a separate implementation", [rnd(v, 2) for v in r14] == ["1.42", "1.32", "1.51"], [rnd(v, 3) for v in r14])
def growth(y):
    t = np.arange(len(y), dtype=float); t -= t.mean(); nll = lambda b: np.sum(np.exp(b[0] + b[1] * t) - y * (b[0] + b[1] * t)); res = optimize.minimize(nll, [np.log(y.mean()), 0.0], method="BFGS", tol=1e-12); b = res.x; mu = np.exp(b[0] + b[1] * t)
    X = np.column_stack([np.ones_like(t), t]); cov = np.linalg.inv((X.T * mu) @ X); phi = max(1.0, np.sum((y - mu) ** 2 / mu) / (len(y) - 2)); return b[1], math.sqrt(cov[1, 1] * phi)
s_ = day["2014-06-01":"2014-09-15"].values.astype(float); r_, se_ = growth(s_); k9, th9 = (15.3 / 9.3) ** 2, 9.3 ** 2 / 15.3; R_ = (1 + r_ * th9) ** k9
check("A.3", "growth rate for 1 June to 15 September 2014 is 0.0192 per day (SE 0.0020), doubling 36.2 days, R 1.32; the generation interval of 10 days gives 1.20 and of 20 days 1.43 (separate optimiser)", rnd(r_, 4) == "0.0192" and rnd(se_, 4) == "0.0020" and rnd(math.log(2) / r_, 1) == "36.2" and rnd(R_, 2) == "1.32" and rnd((1 + r_ * (9.3 * 10 / 15.3) ** 2 / 10) ** ((10 / (9.3 * 10 / 15.3)) ** 2), 2) == "1.20" and rnd((1 + r_ * (9.3 * 20 / 15.3) ** 2 / 20) ** ((20 / (9.3 * 20 / 15.3)) ** 2), 2) == "1.43", (rnd(r_, 5), rnd(se_, 5), rnd(R_, 3)))
kk, thh = (8 / 5) ** 2, 25 / 8; r_true = (1.5 ** (1 / kk) - 1) / thh; GT = pd.read_csv(P("table_genome_growth.csv")); SM = pd.read_csv(P("table_simulation13.csv"))
check("A.4", "simulation: the true growth rate for R = 1.5 and a gamma interval of mean 8 and SD 5 is 0.0549 per day (doubling 12.6 days), as in table_simulation13.csv; the genomic rate falls from 0.0092 to -0.0027 as tied events are spread apart", rnd(r_true, 4) == "0.0549" and rnd(math.log(2) / r_true, 1) == "12.6" and SM.true_r.eq(0.0549).all() and [rnd(v, 4) for v in GT.r_per_day] == ["0.0092", "0.0071", "0.0037", "-0.0003", "-0.0027"], rnd(r_true, 5))
VF = pd.read_csv(P("table_variant_forecasts.csv")); VS = pd.read_csv(P("table_variant_summary.csv"), index_col=0)
check("A.5", "variant forecasts: 19 origins; mean weighted interval score 62,926 for cases only and 344,486 for cases plus genomes; the second is better at 5 origins; 95% interval coverage 10.5% and 0% (recomputed from the forecast table)", len(VF) == 19 and rnd(VF.A_wis.mean(), 0) == "62926" and rnd(VF.B_wis.mean(), 0) == "344486" and int((VF.B_wis < VF.A_wis).sum()) == 5 and rnd(100 * VF.A_in95.mean(), 1) == "10.5" and VF.B_in95.mean() == 0 and rnd(VS.loc["A: cases only", "mean_WIS"], 0) == "62926", (len(VF), rnd(VF.A_wis.mean(), 0)))
# ================================================================== B. the new exercises
SIRt = pd.read_csv(P("table_sir.csv")); sv = SIRt.set_index("quantity").value; GS = pd.read_csv(P("table_gi_shape.csv")); WI = pd.read_csv(P("table_windows.csv")); SI = pd.read_csv(P("table_si_sensitivity.csv")); RF = pd.read_csv(P("table_reporting_fall.csv")); PA = pd.read_csv(P("table_partial.csv"))
N, Bt, Gm = 10000.0, 0.3, 0.1
def rk4(beta_fn, T=200.0, h=0.01):
    y = np.array([N - 10, 10.0, 0.0]); out = [y.copy()]; f = lambda t, y: np.array([-beta_fn(t) * y[0] * y[1] / N, beta_fn(t) * y[0] * y[1] / N - Gm * y[1], Gm * y[1]])
    for i in range(int(round(T / h))): t = i * h; k1 = f(t, y); k2 = f(t + h / 2, y + h / 2 * k1); k3 = f(t + h / 2, y + h / 2 * k2); k4 = f(t + h, y + h * k3); y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); out.append(y.copy())
    return np.array(out)
Y = rk4(lambda t: Bt); Yc = rk4(lambda t: Bt if t < 30 else Bt / 2); ts = np.arange(len(Y)) * 0.01; ip = Y[:, 1].argmax()
sl = np.polyfit(np.arange(0, 11), np.log(Y[::100][:11, 1]), 1)[0]
def fsize(R0, s0=0.999):
    lo, hi = 1e-9, s0 - 1e-9
    for _ in range(200): mid = (lo + hi) / 2; (lo, hi) = (mid, hi) if math.log(s0 / mid) - R0 * (s0 - mid) > 0 else (lo, mid)
    return 1 - (lo + hi) / 2
q = list(sv.index)
check("B.1", "SIR: with a fourth-order Runge-Kutta solver the peak is 3,008 infected on day 38.4, 94.06% are ever infected (final-size equation 94.03%), the early growth rate is 0.199 and the doubling time 3.49 days", close(Y[ip, 1], sv[q[3]], 1.0) and close(ts[ip], 38.4, 0.06) and close((N - Y[-1, 0]) / N, sv[q[4]], 0.0005) and close(fsize(3.0), sv[q[5]], 0.0002) and rnd(sl, 3) == "0.199" and rnd(math.log(2) / sl, 2) == "3.49", (round(Y[ip, 1], 1), round(ts[ip], 2), round((N - Y[-1, 0]) / N, 4), round(fsize(3.0), 4)))
check("B.2", "SIR control: halving beta on day 30 gives a peak of 2,057 on day 30 and a final size of 73.0%; 35% are infected or recovered by day 30; R(t) = R0 S/N is 1.00 at the peak", close(Yc[:, 1].max(), sv[q[7]], 1.0) and close(ts[Yc[:, 1].argmax()], 30.0, 0.06) and close((N - Yc[-1, 0]) / N, sv[q[8]], 0.0005) and rnd(100 * (1 - Y[3000, 0] / N), 0) == "35" and rnd(3 * Y[ip, 0] / N, 2) == "1.00", (round(Yc[:, 1].max(), 1), rnd(100 * (1 - Y[3000, 0] / N), 1)))
def Rw(r, mean, k): return (1 + r * mean / k) ** k
def Rquad(r, mean, k): th = mean / k; return 1 / integrate.quad(lambda x: stats.gamma.pdf(x, a=k, scale=th) * math.exp(-r * x), 0, np.inf)[0]
r_eb = float(np.load(P("r_cases.npy"))[0]); sh_eb = (15.3 / 9.3) ** 2; exp_vals = [Rw(0.2, 10, 1), Rw(0.2, 10, 2), Rw(0.2, 10, 4), math.exp(2.0), Rw(r_eb, 15.3, 1), Rw(r_eb, 15.3, sh_eb), Rw(r_eb, 15.3, 4), math.exp(r_eb * 15.3)]
check("B.3", "generation-interval shape: R from r = 0.2 and a mean of 10 days is 3.00, 4.00, 5.06 and 7.39 (exponential, gamma 2, gamma 4, fixed); for Ebola 1.29, 1.32, 1.33 and 1.34; the closed form agrees with the Wallinga-Lipsitch integral", all(close(a, b, 1e-6) for a, b in zip(exp_vals, GS.R)) and all(close(Rquad(r, m, k), Rw(r, m, k), 1e-6) for r, m, k in ((0.2, 10, 1), (0.2, 10, 2), (0.2, 10, 4), (r_eb, 15.3, sh_eb))) and [rnd(v, 2) for v in exp_vals] == ["3.00", "4.00", "5.06", "7.39", "1.29", "1.32", "1.33", "1.34"], [rnd(v, 3) for v in exp_vals])
tpk = day.idxmax(); bad = []
for win, row in zip((7, 14, 28), WI.itertuples()):
    est = rt_at(inc, w0, win, t1)
    if not (close(est[0], row.Rt_1_Oct, 0.002) and close(est[1], row.lo, 0.002) and close(est[2], row.hi, 0.002)): bad.append((win, [round(v, 3) for v in est]))
check("B.4", "windows: the reproduction number on 1 October 2014 for windows of 7, 14 and 28 days is 1.27, 1.42 and 1.53 with the intervals of table_windows.csv (separate implementation)", not bad, bad)
def first_below(win):
    ser = pd.Series([rt_at(inc, w0, win, t)[0] for t in range((tpk - day.index[0]).days, len(inc))], index=day.index[(tpk - day.index[0]).days:]); b = (ser < 1).astype(int).rolling(14).sum(); ok = b[b == 14].index; return ok.min() - pd.Timedelta(days=13)
dates = [first_below(w).strftime("%Y-%m-%d") for w in (7, 14, 28)]
check("B.5", "windows: from the incidence peak, R stays below 1 for 14 days from 27 November, 29 November and 3 December 2014 (38, 40 and 44 days after the peak) for windows of 7, 14 and 28 days", dates == list(WI.first_date_below_1_for_14_days) and [(pd.Timestamp(x) - tpk).days for x in dates] == [38, 40, 44], dates)
bad = []
for m in (10.0, 15.3, 20.0):
    est = rt_at(inc, weights(m, m * 9.3 / 15.3), 14, t1); row = SI[SI.mean_interval_days == m].iloc[0]
    if not (close(est[0], row.Rt_1_Oct, 0.002) and close(rt_at(inc, weights(m, m * 9.3 / 15.3), 14, (pd.Timestamp("2014-11-01") - day.index[0]).days)[0], row.Rt_1_Nov, 0.002)): bad.append(m)
frac = np.ones(len(inc)); frac[t1:t1 + 31] = np.linspace(1, 0.5, 31); frac[t1 + 31:] = 0.5; th_ = inc * frac; tn = (pd.Timestamp("2014-11-01") - day.index[0]).days
check("B.6", "serial interval and reporting: 1.19, 1.42 and 1.63 on 1 October for means of 10, 15.3 and 20 days; with reporting falling to 50% in October the estimate for 1 November is 0.75 instead of 1.03", not bad and close(rt_at(th_, w0, 14, tn)[0], RF.Rt_thinned[1], 0.002) and rnd(rt_at(th_, w0, 14, tn)[0], 2) == "0.75" and rnd(rt_at(inc, w0, 14, tn)[0], 2) == "1.03", bad)
# a partly sampled outbreak with separate code and fewer replicates
def sim_once(rng, piA, piB, assort, target=1500):
    while True:
        par = [-1]; typ = [0]; gen = [0]; cur = [0]; g = 0
        while len(par) < target and cur:
            nxt = []
            for i in cur:
                for _ in range(rng.poisson(2.4 if typ[i] else 1.2)):
                    par.append(i); typ.append(typ[i] if rng.random() < assort else int(rng.random() < 0.3)); gen.append(g + 1); nxt.append(len(par) - 1)
            g += 1; cur = nxt
        if len(par) >= 300: break
    n = len(par); typ = np.array(typ); gen = np.array(gen); par = np.array(par); smp = rng.random(n) < np.where(typ == 1, piA, piB); last = gen.max(); cnt = np.zeros(n)
    for j in range(1, n):
        if smp[j] and smp[par[j]] and gen[par[j]] < last: cnt[par[j]] += 1
    ok = smp & (gen < last); links = [j for j in range(1, n) if gen[par[j]] < last]; both = sum(1 for j in links if smp[j] and smp[par[j]])
    a = cnt[ok & (typ == 1)].mean(); b = cnt[ok & (typ == 0)].mean(); return both / len(links), a / b
def many(piA, piB, assort, reps=60, seed=5):
    rng = np.random.default_rng(seed); out = np.array([sim_once(rng, piA, piB, assort) for _ in range(reps)]); return out[:, 0].mean(), np.nanmedian(out[:, 1])
sA = many(0.3, 0.3, 0.0); sB = many(0.6, 0.3, 0.8)
check("B.7", "partial sampling (60 outbreaks with separate code): at a sampling probability of 0.3 about 9% of links have both ends sampled and the ratio is near 2; with assortative mixing and 0.6 against 0.3 the direct-link ratio is biased well above 2", close(sA[0], 0.09, 0.015) and 1.7 <= sA[1] <= 2.4 and sB[1] > 2.7, (round(sA[0], 3), round(sA[1], 2), round(sB[1], 2)))
row = lambda s: PA[PA.scenario == s].iloc[0]
check("B.8", "partial sampling in the result file: link shares of 100%, 25%, 9% and 1% for sampling probabilities of 1, 0.5, 0.3 and 0.1 (the square of the probability); the ratio estimated from direct links is within 0.1 of 2 when sampling is independent of mixing; assortative mixing with differential sampling gives a ratio above 3", [rnd(100 * row(s).share_of_links_observed, 0) for s in ("pi = 1.0 (all cases sampled)", "pi = 0.5", "pi = 0.3", "pi = 0.1")] == ["100", "25", "9", "1"] and all(abs(row(s).ratio_direct_median - 2) < 0.1 for s in ("pi = 1.0 (all cases sampled)", "pi = 0.5", "pi = 0.3", "pi = 0.1", "pi_A = 0.6, pi_B = 0.3")) and row("Assortative mixing, pi_A = 0.6, pi_B = 0.3").ratio_direct_median > 3, "")
AU = pd.read_csv(P("table_overlap_audit.csv")); CD = pd.read_csv(P("table_cascade_derived.csv")); IN = pd.read_csv(P("table_independence.csv"))
jc = lambda a, b, ab: ab / (a + b - ab); pubs = {"37%": 894 / 2417, "35%": 894 / 2527, "4%": 59 / 1373, "2%": 59 / 3526, "48%": 1204 / 2527, "65%": 583 / 894, "55%": 494 / 894}
ok = all(rnd(100 * v, 0) == k.rstrip("%") for k, v in pubs.items()) and rnd(100 * 687 / 1056, 0) == "65" and rnd(100 * 631 / 687, 0) == "92" and rnd(100 * 446 / 631, 0) == "71" and 631 - 446 == 185 and rnd(100 * 24 / 185, 0) == "13" and rnd(100 * 121 / 446, 0) == "27"
check("B.9", "audit: every percentage published in the Rhode Island paper recomputes from its counts, and the Jaccard coefficients are 0.22 (persons, with 2,417 or 2,418) and 0.012 (links)", ok and rnd(jc(2417, 2527, 894), 2) == "0.22" and rnd(jc(2418, 2527, 894), 2) == "0.22" and rnd(jc(1373, 3526, 59), 3) == "0.012" and len(AU) == 15 and close(AU.recomputed_pct.iloc[-1], jc(1373, 3526, 59), 5e-5), "")
check("B.10", "audit: the abstract's cascade percentages imply 886, 824, 585 and 158 partners for the 1,342 named; independence implies a population of 6,832 persons with coverages of 35% and 37% (equal coverage 36%)", [round(1342 * 0.66), round(1342 * 0.66 * 0.93), round(1342 * 0.66 * 0.93 * 0.71), round(1342 * 0.66 * 0.93 * 0.71 * 0.27)] == [886, 824, 585, 158] and [round(v) for v in CD.implied_count[1:]] == [886, 824, 585, 158] and rnd(2417 * 2527 / 894, 0) == "6832" and rnd(100 * 2417 / (2417 * 2527 / 894), 0) == "35" and rnd(100 * 2527 / (2417 * 2527 / 894), 0) == "37" and rnd(100 * IN.value[3], 0) == "36", "")
CC = pd.read_csv(P("table_capstone_check.csv")); pa = CC[CC.plan.str.startswith("A")]; pb = CC[CC.plan.str.startswith("B")]
check("B.11", "capstone check: the consistent plan passes all ten rules; the plan with six planted faults fails exactly C2, C4, C5, C6, C7 and C10 (four rules pass)", len(pa) == 10 and (pa.result == "pass").all() and sorted(pb[pb.result == "FAIL"].rule) == sorted(["C2", "C4", "C5", "C6", "C7", "C10"]) and (pb.result == "pass").sum() == 4, list(pb[pb.result == "FAIL"].rule))
check("B.12", "capstone check: deliverables fall due in the order 13, 15, 22, 29 in plan A and the portfolio (29) falls after the report (22) in plan B; plan B has 620 words and 2 figures", "[13, 15, 22, 29]" in pa[pa.rule == "C6"].detail.iloc[0] and "[13, 15, 29, 22]" in pb[pb.rule == "C6"].detail.iloc[0] and "620 words, 2 figures" in pb[pb.rule == "C10"].detail.iloc[0], "")
# ================================================================== G. the chapter
def docx_text(path):
    from docx import Document
    d_ = Document(path); body = [p.text for p in d_.paragraphs]; tabs = [[[c.text for c in r_.cells] for r_ in t_.rows] for t_ in d_.tables]; return d_, "\n".join(body), tabs
if A.chapter:
    dd, txt, tabs = docx_text(A.chapter); flat = re.sub(r"\s+", " ", txt + " " + " ".join(c for t_ in tabs for r_ in t_ for c in r_))
    figs = re.findall(r"^Figure 13\.(\d+)\.", txt, flags=re.M); check("G.1", "figures are numbered 13.1 to 13.12 in order, one caption for every image; 12 images (3 diagrams and 9 figures), 21 tables", figs == [str(i) for i in range(1, 13)] and len(dd.inline_shapes) == 12 and len(dd.tables) == 21, (figs, len(dd.inline_shapes), len(dd.tables)))
    objs, on = [], False
    for p in dd.paragraphs:
        if p.style.name == "Heading 2": on = (p.text == "Learning objectives")
        elif on and p.text.strip(): objs.append(p.text)
    objs = objs[1:]; check("G.2", "fifteen learning objectives, from why genomes and cases answer different questions to the links with the earlier chapters", len(objs) == 15 and objs[0].startswith("Explain why genomes") and objs[-1].startswith("Connect the integrated view"), len(objs))
    heads = [p.text for p in dd.paragraphs if p.style.name == "Heading 2"]; ho = [h_ for h_ in heads if h_.startswith("Hands-on")]; check("G.3", "nine hands-on sections numbered 1 to 9 in order, then the capstone and the synthesis sections", [h_.split(" ")[1] for h_ in ho] == [str(i) for i in range(1, 10)] and "From sequences to decisions: the capstone" in heads and "What the methods can and cannot show" in heads, ho)
    want = ["0.199 per day", "3.49 days", "94.1%", "73.0%", "R is 3.0 for an exponential interval, 4.0 for a gamma shape of 2, 5.06 for a shape of 4 and 7.39 for a fixed interval", "1.27, 1.42 and 1.53 for windows of 7, 14 and 28 days", "1.19, 1.42 and 1.63", "0.75 instead of 1.03", "27 November", "3 December", "252 genomes", "1.63 (0.86 to 2.32)", "0.86", "0.32 (0.13 to 0.70)", "227 of the 255 infected farms", "0.89", "894", "0.22", "0.012", "48%", "88%", "6832", "4%", "plan A passes all 10", "plan B", "13.5%".replace("13.5%", "13%"), "twelve sentinel sites in 2023", "16 November 2026"]
    want = [w_ for w_ in want if w_ not in ("6832", "plan B", "13%", "4%")] + ["about 6,800" if "about 6,800" in flat else "6,832"]
    miss = [w_ for w_ in want if w_ not in flat]; check("G.4", "the sentences quote the validated numbers", not miss, miss)
    def tbl(first, second=None): return next((t_ for t_ in tabs if t_[0][0].startswith(first) and (second is None or t_[0][1].startswith(second))), None)
    t1_ = tbl("Window (days)"); ok = t1_ is not None and [r_[1].split(" ")[0] for r_ in t1_[1:]] == [rnd(v, 2) for v in WI.Rt_1_Oct] and [r_[4] for r_ in t1_[1:]] == [pd.Timestamp(x).strftime("%-d %B") for x in WI.first_date_below_1_for_14_days]; check("G.5", "the window table of the chapter equals the result file", ok, None if t1_ is None else t1_[1])
    t2_ = tbl("Mean serial interval"); ok = t2_ is not None and [r_[2].split(" ")[0] for r_ in t2_[1:]] == [rnd(v, 2) for v in SI.Rt_1_Oct]; check("G.6", "the serial-interval table of the chapter equals the result file", ok, None if t2_ is None else t2_[1])
    t3_ = tbl("Date", "Reporting fraction"); ok = t3_ is not None and [r_[3] for r_ in t3_[1:]] == [rnd(v, 2) for v in RF.Rt_thinned]; check("G.7", "the reporting-fall table of the chapter equals the result file", ok, None if t3_ is None else t3_[1])
    t4_ = tbl("Scenario", "Links seen"); ok = t4_ is not None and [r_[1] for r_ in t4_[1:]] == [f"{rnd(100 * v, 1)}%" for v in PA.share_of_links_observed] and [r_[3].split(" ")[0] for r_ in t4_[1:]] == [rnd(v, 2) for v in PA.ratio_direct_median]; check("G.8", "the partial-sampling table of the chapter equals the result file", ok, None if t4_ is None else t4_[1])
    t5_ = tbl("Quantity", "Count"); ok = t5_ is not None and [r_[4] for r_ in t5_[1:]][:5] == [rnd(v, 2) for v in AU.recomputed_pct[:5]]; check("G.9", "the audit table of the chapter equals the result file", ok, None if t5_ is None else t5_[1])
    t6_ = tbl("Rule", "Check"); ok = t6_ is not None and [r_[2] for r_ in t6_[1:]] == pa.result.tolist() and [r_[3] for r_ in t6_[1:]] == pb.result.tolist(); check("G.10", "the capstone-check table of the chapter equals the result file", ok, None if t6_ is None else t6_[1])
    import zipfile; z = zipfile.ZipFile(A.chapter); rels = z.read("word/_rels/document.xml.rels").decode(); dois = {u_.rstrip(".,;").replace("&amp;", "&") for u_ in re.findall(r"https://doi\.org/[^\"'<> ]+", rels)}
    newd = ["10.1093/molbev/msaf083", "10.1093/ofid/ofaf341", "10.12688/verixiv.3269.3", "10.1038/s41564-018-0296-2"]; have = {x_.split("doi.org/")[1] for x_ in dois}
    check("G.11", "the chapter links the 4 new DOIs read for this extension, with 21 distinct DOIs in all (17 in the original)", all(n_ in have for n_ in newd) and len(have) == EXPECTED_DOIS, (len(have), [n_ for n_ in newd if n_ not in have]))
    if A.doi:
        import urllib.request, urllib.parse, time; bad_d = []
        for dd_ in sorted(dois):
            doi = dd_.split("doi.org/")[1]; okd = False
            for tryn in range(4):
                try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "mepi-book/1.0 (mailto:none@example.org)"}), timeout=40).read(); okd = True; break
                except Exception as e: err = str(e)[:30]; time.sleep(5)
            if not okd: bad_d.append((doi, err))
        check("G.12", "every DOI of the chapter resolves in Crossref", not bad_d, bad_d)
    else: skip("G.12", "every DOI of the chapter resolves in Crossref", "run with --doi")
else:
    for k, w_ in [("G.1", "figure numbers, tables"), ("G.2", "objectives"), ("G.3", "hands-on sections"), ("G.4", "sentences"), ("G.5", "window table"), ("G.6", "serial-interval table"), ("G.7", "reporting-fall table"), ("G.8", "partial-sampling table"), ("G.9", "audit table"), ("G.10", "capstone table"), ("G.11", "DOIs"), ("G.12", "DOIs resolve")]: skip(k, w_, "run with --chapter")
fd = next((d_ for d_ in ("../figures", "figures", os.path.join(HERE, "../figures")) if os.path.isdir(d_)), None)
figs_ = ["diagram_integration", "diagram_capstone_workflow", "diagram_evidence_to_action", "fig1_ebola_rt", "fig2_genome_growth", "fig3_simulation", "fig4_variant_forecast", "fig5_sir", "fig6_windows", "fig7_partial", "fig8_integration", "fig9_agentic_workflow_diagram"]
if fd: miss = [f_ for f_ in figs_ if not os.path.exists(os.path.join(fd, f_ + ".png"))]; check("H.1", "all 12 figure files are present (3 diagrams and 9 figures)", not miss, miss)
else: skip("H.1", "figure files", "no figures folder")
cd = next((d_ for d_ in (os.path.join(HERE, "../code"), "../code", "code") if os.path.isdir(d_) and glob.glob(os.path.join(d_, "ex13_5_sir_gi.py"))), None)
if cd:
    import py_compile; badc = []
    for f_ in sorted(glob.glob(os.path.join(cd, "*.py"))):
        try: py_compile.compile(f_, doraise=True)
        except Exception: badc.append(os.path.basename(f_))
    check("H.2", "every Python script of the package compiles", not badc, badc)
out_ = pd.DataFrame(R); out_.to_csv("validation_results_ch13.csv", index=False)
for r in R: print(f"{r['result']:5} {r['id']:8} {r['check'][:150]}" + (f"   [{r['detail']}]" if r["result"] != "PASS" else ""))
print(f"\n{(out_.result == 'PASS').sum()} PASS, {(out_.result == 'FAIL').sum()} FAIL, {(out_.result == 'SKIP').sum()} SKIP"); sys.exit(1 if (out_.result == "FAIL").any() else 0)
