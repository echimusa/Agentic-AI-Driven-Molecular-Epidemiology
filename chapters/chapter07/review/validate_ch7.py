"""
Chapter 7 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts:
Hagelloch (Exercises 1 and 4): serial interval, offspring distribution with its own negative-binomial fit, the closed-form R_t of Cori et al., the timing-only ranking probabilities, generations, family, class and distance patterns, recorded fatality with a score
interval, the effect of the onset definition, and tracing recall and workload for every window and contact definition; Ebola (Exercise 2): the effective numbers of candidates from the ranking file; TransPhylo (Hands-on 5): every metric recomputed
from the raw posterior pair probabilities and the true infectors; and the chapter text (tables, sentences, figure numbers, DOIs).
USAGE   python validate_ch7.py [--chapter Chapter_7_....docx] [--doi] [--full]
        --chapter  also checks the chapter's tables, sentences, figure numbers and references against the results (needs pandoc and python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns Exercises 1, 2 and 4 in a scratch folder and checks that their result files are reproduced (needs mappy, pyreadr, biopython, networkx)
Prints PASS/FAIL/SKIP per check, writes validation_results_ch7.csv and exits 1 if any check fails.
"""
import argparse, glob, hashlib, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd, pyreadr
from scipy import stats, optimize
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); ap.add_argument("--sim", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
# ================================================================== A. Hagelloch: Exercises 1 and 4
d = pyreadr.read_r(P("measles_hagelloch_1861.RData"))["measles_hagelloch_1861"].copy()
d["onset"] = pd.to_datetime(d.date_of_prodrome); d["t"] = (d.onset - d.onset.min()).dt.days; d["infector"] = pd.to_numeric(d.infector); d.index = d.case_ID.astype(int)
n = len(d); inf = d.infector.dropna().astype(int); has = d.infector.notna()
check("A.1", "188 cases, 184 with a documented infector, prodrome from 30 October 1861 to 24 January 1862 (86 days)", n == 188 and int(has.sum()) == 184 and str(d.onset.min().date()) == "1861-10-30" and str(d.onset.max().date()) == "1862-01-24" and int(d.t.max()) == 86)
si = (d.t[has.values].values - d.t.reindex(inf.values).values).astype(float)
check("A.2", "serial interval of the 184 pairs: mean 10.4 d, median 10 d, SD 1.7 d, range 8 to 16 d (no pair below 1 day)", rnd(si.mean(), 1) == "10.4" and np.median(si) == 10 and rnd(si.std(ddof=1), 1) == "1.7" and si.min() >= 1 and int(si.max()) == 16, (si.mean(), np.median(si), si.std(ddof=1), si.min(), si.max()))
cnt = np.bincount(inf.values, minlength=d.index.max() + 1)[d.index.values]; m_ = cnt.mean()
def nll(lk): k = math.exp(lk); return -sum(stats.nbinom.logpmf(c, k, k / (k + m_)) for c in cnt)
kh = math.exp(optimize.minimize_scalar(nll, bounds=(-4, 4), method="bounded").x); top = np.sort(cnt)[::-1]
check("A.3", "offspring distribution recomputed: mean 0.98, maximum 30, negative-binomial dispersion k = 0.25, the top 10% of cases cause 61% of infections, 67% infect nobody", rnd(m_, 2) == "0.98" and cnt.max() == 30 and rnd(kh, 2) == "0.25" and rnd(top[:int(round(0.1 * n))].sum() / top.sum(), 2) == "0.61" and rnd((cnt == 0).mean(), 2) == "0.67", (m_, cnt.max(), kh))
shape, _, scale = stats.gamma.fit(si, floc=0); s_ = np.arange(0, 61); w = stats.gamma.cdf(s_ + 0.5, shape, scale=scale) - stats.gamma.cdf(s_ - 0.5, shape, scale=scale); w[0] = 0; w /= w.sum()
I = np.bincount(d.t.values, minlength=int(d.t.max()) + 1).astype(float); Lam = np.convolve(I, w)[:len(I)]; best = None; rts = []
for t in range(7, len(I)):
    a = 1 + I[t - 6:t + 1].sum(); rate = 0.2 + Lam[t - 6:t + 1].sum()
    if I[t - 6:t + 1].sum() >= 5: rts.append((t, a / rate, *stats.gamma.ppf([.025, .975], a, scale=1 / rate)))
pk = max(rts, key=lambda x: x[1]); below = min(r_[0] for r_ in rts if r_[1] < 1 and r_[0] > pk[0])
check("A.4", "the closed-form R_t of Cori et al. (7-day window, Gamma(1, scale 5) prior, at least 5 cases in the window) peaks at 8.1 (95% interval 6.0 to 10.4) on day 27 and falls below 1 from day 41", rnd(pk[1], 1) == "8.1" and rnd(pk[2], 1) == "6.0" and rnd(pk[3], 1) == "10.4" and pk[0] == 27 and below == 41, (pk, below))
tt = d.t.values; ids = d.index.values; pos = {c: i for i, c in enumerate(ids)}; pt = []; npair = []
for j in range(n):
    if not has.iloc[j]: continue
    dt = tt[j] - tt; ok = (dt >= 1) & (dt < len(w)); wt = np.zeros(n); wt[ok] = w[dt[ok]]; ptrue = wt[pos[int(d.infector.iloc[j])]] / wt.sum(); pt.append(ptrue); npair.append(ok.sum())
check("A.5", "timing-only ranking: the mean probability on the documented infector is 10.5% (0.105), and 1 over the mean number of earlier candidates is below it", rnd(np.mean(pt), 3) == "0.105" and np.mean(pt) > 1 / np.mean(npair), (np.mean(pt), 1 / np.mean(npair)))
# ---- Exercise 4: descriptive epidemiology (own implementation: iterative propagation, vectorised pair matrices)
gen = np.full(n, -1); root = ~has.values
gen[root] = 0
for _ in range(20):
    for j in range(n):
        if gen[j] < 0 and has.iloc[j]:
            g = gen[pos[int(d.infector.iloc[j])]]
            if g >= 0: gen[j] = g + 1
Gn = pd.DataFrame({"g": gen, "t": tt}); G = Gn.groupby("g").agg(cases=("t", "size"), med=("t", "median")); gt = pd.read_csv(P("table_generations.csv")).set_index("generation")
check("A.6", "generations by propagation along the documented infectors: 6 generations of 4, 8, 50, 102, 22 and 2 cases with the median onset days of the table", list(G.cases) == [4, 8, 50, 102, 22, 2] and np.allclose(G.med.values, gt.median_day.values) and (gen >= 0).all(), list(G.cases))
dead = d.date_of_death.notna().values; age = d.age.values.astype(float); grp = np.where(age <= 4, "0-4", np.where(age <= 9, "5-9", "10+")); pp = pd.read_csv(P("table_person.csv")).set_index("group")
def score_ci(k, nn, z=1.959964): ph = k / nn; lo = (ph + z * z / (2 * nn) - z * math.sqrt(ph * (1 - ph) / nn + z * z / (4 * nn * nn))) / (1 + z * z / nn); hi = (ph + z * z / (2 * nn) + z * math.sqrt(ph * (1 - ph) / nn + z * z / (4 * nn * nn))) / (1 + z * z / nn); return lo, hi
okp = int(dead.sum()) == 12 and pp.loc["all", "deaths"] == 12 and np.median(age) == 7 and np.quantile(age, .25) == 3 and np.quantile(age, .75) == 11
for g_ in ("0-4", "5-9", "10+"):
    k, nn = int(dead[grp == g_].sum()), int((grp == g_).sum()); lo, hi = score_ci(k, nn); r = pp.loc[f"age {g_}"]; okp &= r.deaths == k and r.n == nn and abs(r.cfr_lo - lo) < 2e-4 and abs(r.cfr_hi - hi) < 2e-4
check("A.7", "person: median age 7 (IQR 3 to 11); 12 recorded deaths of 188; fatality by age group (0-4 10.1%, 5-9 4.5%, 10+ 4.0%) with score intervals recomputed to 0.0002; the complications column is 'yes' for every case", okp and (d.complications.astype(str) == "yes").all() and rnd(dead[grp == "0-4"].mean(), 3) == "0.101", d.complications.value_counts().to_dict())
fam = d.family_ID.values; cls = d["class"].astype(str).to_numpy(); X = d.x_loc.values; Y = d.y_loc.values; later = tt[None, :] > tt[:, None]                         # [i, j]: i earlier than j
sf = (fam[:, None] == fam[None, :]); sc = (cls[:, None] == cls[None, :]); dist = np.hypot(X[:, None] - X[None, :], Y[:, None] - Y[None, :]); pl = pd.read_csv(P("table_place.csv")).set_index("measure").value
ii = np.array([pos[int(x)] for x in d.infector.dropna().astype(int)]); jj = np.where(has.values)[0]; fc = pd.Series(fam).value_counts()
okpl = int(len(fc)) == 69 and int((fc >= 2).sum()) == 53 and abs(fc[fc >= 2].sum() / n - pl["share of cases in those families"]) < 1e-4 and abs(sf[ii, jj].mean() - pl["documented pairs in the same family"]) < 1e-4 and abs(sf[later].mean() - pl["all earlier-later pairs in the same family"]) < 1e-4
okpl &= abs(sc[ii, jj].mean() - pl["documented pairs in the same class"]) < 1e-4 and abs(sc[later].mean() - pl["all earlier-later pairs in the same class"]) < 1e-4 and abs(np.median(dist[ii, jj]) - pl["median distance, documented pairs (map units)"]) < 1e-3 and abs(np.median(dist[later]) - pl["median distance, all earlier-later pairs (map units)"]) < 1e-3
check("A.8", "place: 69 families, 53 with two or more cases (91.5% of cases); the documented infector shares the family in 41.9% of pairs (1.3% of all earlier-later pairs) and the class in 59.2% (37.3%); median distance 10.2 against 107.7 map units", okpl, (sf[ii, jj].mean(), sf[later].mean()))
rash = pd.to_datetime(d.date_of_rash); lag = (rash - d.onset).dt.days.values; sir = (rash[has.values].values - rash.reindex(inf.values).values) / np.timedelta64(1, "D"); rt_ = ((rash - d.onset.min()).dt.days).value_counts()
sm = pd.read_csv(P("table_descriptive_summary.csv")).set_index("key").value
check("A.9", "onset definition: the rash follows the prodrome by a median of 4 days (IQR 3 to 5); dating onset by the rash moves the peak from day 32 to day 38 and the serial interval from 10.4 (SD 1.7) to 10.5 (SD 2.5) days", np.median(lag) == 4 and np.quantile(lag, .25) == 3 and np.quantile(lag, .75) == 5 and int(pd.Series(tt).value_counts().idxmax()) == 32 and int(rt_.idxmax()) == 38 and rnd(si.mean(), 1) == "10.4" and rnd(sir.mean(), 1) == "10.5" and rnd(sir.std(ddof=1), 1) == "2.5" and abs(float(sm["si_mean_rash"]) - sir.mean()) < 1e-6, (sir.mean(), sir.std(ddof=1)))
tab = pd.read_csv(P("table_tracing.csv")); okt = True; dfun = {"anyone": np.ones((n, n), bool), "same family": sf, "same class": sc, "family or class": sf | sc, "within 25 map units": dist <= 25, "within 100 map units": dist <= 100}; dtm = tt[None, :] - tt[:, None]                    # [i, j]: onset j minus onset i
for W in (7, 10, 14, 21):
    for name, M in dfun.items():
        cand = (dtm > 0) & (dtm <= W) & M; rec = np.mean([cand[pos[int(d.infector.iloc[j])], j] for j in jj]) ; r = tab[(tab.window_days == W) & (tab.contact_definition == name)].iloc[0]
        okt &= abs(rec - r.recall) < 1e-4 and abs(cand.sum(axis=0).mean() - r.candidates_per_case) < 1e-3
check("A.10", "tracing: the recall of the documented infector and the candidates per case, for 4 windows and 6 contact definitions, recomputed from boolean pair matrices equal table_tracing.csv", okt)
q = np.percentile(si, [50, 90, 95, 99]); r14 = tab[tab.window_days == 14].set_index("contact_definition")
check("A.11", "serial interval quantiles 10, 12.7, 13 and 15 days (maximum 16); a 14-day window finds 98.4% of infectors with 73.0 candidates per case when anyone counts, 89.1% with 28.8 for family or class, 40.8% with 1.2 for family alone; windows of 7, 10 and 21 days find 2.2%, 57.1% and 100.0%", rnd(q[0], 0) == "10" and rnd(q[1], 1) == "12.7" and rnd(q[2], 0) == "13" and rnd(q[3], 0) == "15" and [rnd(100 * r14.loc[x, "recall"], 1) for x in ("anyone", "family or class", "same family")] == ["98.4", "89.1", "40.8"] and [rnd(r14.loc[x, "candidates_per_case"], 1) for x in ("anyone", "family or class", "same family")] == ["73.0", "28.8", "1.2"] and [rnd(100 * tab[(tab.window_days == w_) & (tab.contact_definition == "anyone")].recall.iloc[0], 1) for w_ in (7, 10, 21)] == ["2.2", "57.1", "100.0"], q)
# ================================================================== B. Ebola: Exercise 2 (from the ranking file)
rk = pd.read_csv(P("ebola_candidate_ranking.csv"))
check("B.1", "Ebola ranking file: 61 early cases; the median effective number of candidate infectors is 21.0 with timing only and 12.5 with timing plus genetics; genetics lowers it for most cases", len(rk) == 61 and rnd(rk.neff_time.median(), 1) == "21.0" and rnd(rk.neff_gen.median(), 1) == "12.5" and (rk.neff_gen <= rk.neff_time + 1e-9).mean() > 0.9, ((rk.neff_gen <= rk.neff_time + 1e-9).mean()))
# ================================================================== C. Hands-on 5: TransPhylo scored from the raw pair probabilities
runs = pd.read_csv(P("transphylo_runs.csv")); truth = pd.read_csv(P("transphylo_truth.csv")); pairs = pd.read_csv(P("transphylo_pairs.csv")); sruns = pd.read_csv(P("transphylo_short_runs.csv")); spairs = pd.read_csv(P("transphylo_short_pairs.csv")); strue = pd.read_csv(P("transphylo_short_truth.csv"))
tsum = pd.read_csv(P("table_transphylo_summary.csv")).set_index("pi_true"); ssum = pd.read_csv(P("table_transphylo_short_summary.csv")).set_index("pi_true")
okd = len(runs) == 12 and sorted(runs.pi_true.unique()) == [0.25, 0.75] and (runs.groupby("pi_true").size() == 6).all() and runs.sampled.between(15, 35).all() and (runs.iterations == 100000).all() and (runs.thinning == 10).all() and (runs.samples == 10000).all() and (runs.burnin_samples == 5000).all() and (sruns.iterations == 2000).all() and (sruns.sampled == runs.sampled).all() and (sruns.hosts == runs.hosts).all()
for r_ in runs.itertuples():
    t_ = truth[truth.run == r_.run].sort_values("host"); ns = int(r_.sampled); okd &= len(t_) == r_.hosts and t_.sampling_time.iloc[:ns].notna().all() and t_.sampling_time.iloc[ns:].isna().all() and t_.infector.between(0, len(t_)).all()
check("C.1", "12 simulated outbreaks, 6 at each sampling probability (0.25 and 0.75), 15 to 35 sampled hosts each; the sampled hosts are rows 1 to n of the truth table; the long chains have 100,000 iterations thinned to 10,000 samples with 5,000 discarded; the short chains are the same outbreaks at 2,000 iterations", okd)
def metrics(pr, tr, rn):
    out = []
    for k in rn.run:
        ns = int(rn[rn.run == k].sampled.iloc[0]); t_ = tr[tr.run == k].set_index("host"); pp = pr[pr.run == k]; prob = {(int(i), int(j)): v for i, j, v in zip(pp.i, pp.j, pp.prob)}
        true_inf = {j: int(t_.infector[j]) for j in range(1, ns + 1) if 1 <= int(t_.infector[j]) <= ns}
        links = len(true_inf); corr = 0
        for j, ti in true_inf.items():
            col = [prob.get((i, j), 0.0) for i in range(1, ns + 1)]
            if sum(col) > 0 and (int(np.argmax(col)) + 1) == ti: corr += 1
        called = [(i, j) for (i, j), v in prob.items() if v > 0.5]; tc = sum(1 for (i, j) in called if true_inf.get(j) == i)
        pn_u = np.mean([1 - sum(prob.get((i, j), 0.0) for i in range(1, ns + 1)) for j in range(1, ns + 1) if j not in true_inf]) if ns > links else np.nan; pn_s = np.mean([1 - sum(prob.get((i, j), 0.0) for i in range(1, ns + 1)) for j in true_inf]) if links else np.nan
        out.append(dict(pi=float(rn[rn.run == k].pi_true.iloc[0]), ns=ns, links=links, top1=corr, called=len(called), true_called=tc, pnu=pn_u, pns=pn_s))
    return pd.DataFrame(out)
M = metrics(pairs, truth, runs); agg = M.groupby("pi").apply(lambda g: pd.Series(dict(links=g.links.sum(), top1=g.top1.sum(), called=g.called.sum(), tc=g.true_called.sum(), hosts=g.ns.sum(), pnu=g.pnu.mean(), pns=g.pns.mean())), include_groups=False)
okm = True
for pv in (0.25, 0.75):
    a, t = agg.loc[pv], tsum.loc[pv]; okm &= a.links == t.links and a.top1 == t.top1_n and a.called == t.pairs_called and a.tc == t.true_called and abs(a.tc / a.called - t.precision) < 1e-3 and abs(a.tc / a.links - t.recall) < 1e-3 and abs(a.top1 / a.links - t.pooled_top1) < 1e-3 and abs(a.pnu - t.p_none_if_unsampled) < 1e-3 and abs(a.pns - t.p_none_if_sampled) < 1e-3
check("C.2", "TransPhylo scored from the raw pair probabilities and true infectors with separate code: true links, top-1 correct, links called above 0.5, true among them, and the mean probability of no sampled infector, per sampling probability, equal table_transphylo_summary.csv", okm, agg.round(3).to_dict())
check("C.3", "pooled results: at a sampling probability of 0.25, 23 of 110 hosts have a sampled infector, all 23 are top-ranked correctly, 24 links are called of which 18 are true (precision 75%, recall 78%), and the probability of no sampled infector is 0.84 when the truth is unsampled against 0.30 when sampled; at 0.75, 74 of 100 hosts, top-1 70%, 49 called and 42 true (precision 86%, recall 57%)", [int(agg.loc[0.25, c]) for c in ("links", "hosts", "top1", "called", "tc")] == [23, 110, 23, 24, 18] and [int(agg.loc[0.75, c]) for c in ("links", "hosts", "top1", "called", "tc")] == [74, 100, 52, 49, 42] and rnd(agg.loc[0.25, "pnu"], 2) == "0.84" and rnd(agg.loc[0.25, "pns"], 2) == "0.30" and rnd(agg.loc[0.75, "tc"] / agg.loc[0.75, "called"], 2) == "0.86" and rnd(agg.loc[0.75, "tc"] / agg.loc[0.75, "links"], 2) == "0.57")
Ms = metrics(spairs, strue, sruns); aggs = Ms.groupby("pi").apply(lambda g: pd.Series(dict(called=g.called.sum(), tc=g.true_called.sum())), include_groups=False)
check("C.4", "effective sample sizes: every 100,000-iteration chain has at least 100 (smallest 106), every 2,000-iteration chain at most 20 (smallest 4.3); the precision of the called links at 0.25 falls from 75% to 59% with the short chains", runs.ess_log_pTTree.min() >= 100 and rnd(runs.ess_log_pTTree.min(), 0) == "106" and sruns.ess_log_pTTree.max() <= 20 and rnd(sruns.ess_log_pTTree.min(), 1) == "4.3" and rnd(runs.ess_log_pTTree.min(), 0) == rnd(tsum.ess_min.min(), 0) and rnd(aggs.loc[0.25, "tc"] / aggs.loc[0.25, "called"], 2) == "0.59", (runs.ess_log_pTTree.min(), sruns.ess_log_pTTree.min(), sruns.ess_log_pTTree.max()))
pm = runs.groupby("pi_true")[["pi_post", "off_r_post", "neg_post"]].mean()
check("C.5", "posterior parameters: the sampling probability is 0.33 for a true 0.25 and 0.68 for 0.75, the offspring parameter 5.6 and 4.9 (true 5), and the within-host scale 0.55 and 0.68 (true 0.27, a factor of two or more too high)", [rnd(pm.loc[0.25, "pi_post"], 2), rnd(pm.loc[0.75, "pi_post"], 2), rnd(pm.loc[0.25, "off_r_post"], 1), rnd(pm.loc[0.75, "off_r_post"], 1), rnd(pm.loc[0.25, "neg_post"], 2), rnd(pm.loc[0.75, "neg_post"], 2)] == ["0.33", "0.68", "5.6", "4.9", "0.55", "0.68"] and bool((pm.neg_post / (100 / 365) >= 2).all()))
tabc = []; allp = []
for k in runs.run:
    ns = int(runs[runs.run == k].sampled.iloc[0]); t_ = truth[truth.run == k].set_index("host"); pp = pairs[pairs.run == k]
    for i, j, v in zip(pp.i, pp.j, pp.prob): allp.append((float(runs[runs.run == k].pi_true.iloc[0]), v, 1 if int(t_.infector[j]) == i else 0))
al = pd.DataFrame(allp, columns=["pi", "p", "true"]); al = al[al.p >= 0.01]; al["bin"] = pd.cut(al.p, [0.0099, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0001], include_lowest=True)
g = al.groupby("bin", observed=True).agg(n=("p", "size"), pred=("p", "mean"), obs=("true", "mean")); gap = float((g.n * (g.pred - g.obs).abs()).sum() / g.n.sum()); gp = {pv: al[al.pi == pv].groupby("bin", observed=True).agg(n=("p", "size"), pred=("p", "mean"), obs=("true", "mean")) for pv in (0.25, 0.75)}
gapp = {pv: float((x.n * (x.pred - x.obs).abs()).sum() / x.n.sum()) for pv, x in gp.items()}
check("C.6", "calibration recomputed from every pair with probability of at least 0.01: 390 pairs, an average gap between predicted and observed share of 0.02 overall (0.06 at 0.25 and 0.05 at 0.75), and 97 true links in all", len(al) == 390 and rnd(gap, 2) == "0.02" and rnd(gapp[0.25], 2) == "0.06" and rnd(gapp[0.75], 2) == "0.05" and int(tsum.links.sum()) == 97, (len(al), gap, gapp))
# ================================================================== D. the chapter text (needs --chapter)
if A.chapter:
    from docx import Document
    dx = Document(A.chapter); paras = [p.text for p in dx.paragraphs]; Tb = dx.tables; alltext = "\n".join(paras) + "\n" + "\n".join(c.text for t in Tb for r in t.rows for c in r.cells); alltext_n = alltext.replace("\u2019", "'").replace("\u201c", '"').replace("\u201d", '"')
    md = subprocess.run(["pandoc", A.chapter, "-t", "markdown", "--wrap=none"], capture_output=True, text=True).stdout
    h2 = [p.text for p in dx.paragraphs if p.style.name.startswith("Heading 2")]
    want = ["Learning objectives", "At a glance", "Core concepts", "Genetic distance and clustering", "Inference methods", "Reproduction number and superspreading", "Outbreaks, cases and the line list", "Tracing and containment zones", "Hands-on 1 (Python)", "Hands-on 2 (Python)", "Hands-on 3 (Python)", "Hands-on 4 (R)", "Hands-on 5 (R)", "Outbreak investigation workflow", "Naming a source", "Case studies revisited", "Agentic AI in outbreak investigation", "Common pitfalls", "Links to other chapters", "Sources"]
    pos_ = [next((i for i, h in enumerate(h2) if h.startswith(w)), -1) for w in want]
    check("D.1", "the 20 main sections are present in the intended order (concepts, outbreaks and the line list, tracing, Hands-on 1 to 5, workflow, naming a source, case studies, agent, pitfalls, links, sources)", all(x >= 0 for x in pos_) and pos_ == sorted(pos_), [w for w, x in zip(want, pos_) if x < 0])
    caps = [int(m) for m in re.findall(r"^Figure 7\.(\d+)\. ", alltext, flags=re.M)]; zx = __import__("zipfile").ZipFile(A.chapter).read("word/document.xml").decode("utf8", "ignore"); descr = re.findall(r'<wp:docPr [^>]*descr="[^"]{10,}"', zx)
    check("D.2", "figure captions run from 7.1 to 7.16 in order, one for every picture, and all 16 pictures have alt text", caps == list(range(1, 17)) and len(descr) == 16 and zx.count("<wp:docPr ") == 16, (caps, len(descr)))
    nobj = len(re.findall(r"^\d+\.\s+", md.split("## Learning objectives")[1].split("## At a glance")[0], flags=re.M))
    check("D.3", "the front matter is clean: no 'Nov 8, 2026' export line, no uploaded-file or draft mention, the title is set; the learning objectives list 10 items", "Nov 8, 2026" not in alltext and re.search(r"uploaded (draft|file|files|script|document)", alltext.lower()) is None and dx.core_properties.title.startswith("Chapter 7") and nobj == 10, nobj)
    cell = lambda t: [[c.text.strip() for c in r.cells] for r in t.rows]; pct1 = lambda x: f"{100 * x:.1f}%"; pct0 = lambda x: f"{100 * x:.0f}%"
    gt = pd.read_csv(P("table_generations.csv")); tab_g = [t for t in Tb if t.rows[0].cells[0].text == "Generation"][0]
    exp = [[str(int(r.generation)), str(int(r.cases)), f"{r.median_day:g}", f"{int(r.first_day)} to {int(r.last_day)}", "" if np.isnan(r.gap_to_previous_median) else f"{r.gap_to_previous_median:.1f}"] for r in gt.itertuples()]
    check("D.4", "the generation table of Hands-on 3 (cases, median onset day, onset span and gap between generations) equals table_generations.csv", cell(tab_g)[1:] == exp, cell(tab_g)[1:2])
    pe = pd.read_csv(P("table_person.csv")); tab_p = [t for t in Tb if t.rows[0].cells[0].text == "Group" and "Recorded" in t.rows[0].cells[2].text][0]
    exp = [[r.group, str(int(r.n)), str(int(r.deaths)), f"{pct1(r.cfr)} ({pct1(r.cfr_lo)} to {pct1(r.cfr_hi)})"] for r in pe.itertuples()]
    check("D.5", "the person table (cases, recorded deaths, recorded fatality with its 95% score interval, for all cases and by age, sex and class) equals table_person.csv, with no negative interval", cell(tab_p)[1:] == exp and "-0.0%" not in alltext and (pe.cfr_lo >= 0).all(), cell(tab_p)[1:2])
    pl = pd.read_csv(P("table_place.csv")); tab_pl = [t for t in Tb if t.rows[0].cells[0].text == "Measure" and t.rows[0].cells[1].text == "Value"][0]
    fmt = lambda m, v: pct1(v) if (m.startswith("share") or "same" in m) else (f"{v:.1f}" if "distance" in m else str(int(v)))
    check("D.6", "the place table (families, shares of cases and pairs in the same family or class, median distances) equals table_place.csv", cell(tab_pl)[1:] == [[r.measure, fmt(r.measure, r.value)] for r in pl.itertuples()], cell(tab_pl)[1:2])
    trc = pd.read_csv(P("table_tracing.csv")); tab_t = [t for t in Tb if t.rows[0].cells[0].text == "Contact definition"][0]; expt = []
    for d_ in ("anyone", "same family", "same class", "family or class", "within 25 map units", "within 100 map units"):
        expt.append([d_] + [f"{100 * trc[(trc.window_days == w) & (trc.contact_definition == d_)].recall.iloc[0]:.0f}% ({trc[(trc.window_days == w) & (trc.contact_definition == d_)].candidates_per_case.iloc[0]:.1f})" for w in (7, 10, 14, 21)])
    check("D.7", "the tracing table (share of documented infectors found, with the candidates per case in brackets, for four windows and six contact definitions) equals table_tracing.csv", cell(tab_t)[1:] == expt, cell(tab_t)[1:2])
    sm = pd.read_csv(P("table_descriptive_summary.csv")).set_index("key").value; gi = gt.set_index("generation"); f_ = lambda k: float(sm[k])
    tr_ = lambda w, d_, c: float(trc[(trc.window_days == w) & (trc.contact_definition == d_)][c].iloc[0])
    sent = [f"from {sm['onset_first']} to {sm['onset_last']}, with a peak of {int(f_('peak_cases'))} cases on day {int(f_('peak_day'))}", f"have median onsets on days {gi.loc[2].median_day:g} and {gi.loc[3].median_day:g}, {gi.loc[3].median_day - gi.loc[2].median_day:.0f} days apart", f"The median age is {f_('median_age'):.0f} years (interquartile range {f_('age_q1'):.0f} to {f_('age_q3'):.0f})", f"for {int(f_('recorded_deaths'))} of the {int(f_('cases'))} cases", f"{100 * pl.set_index('measure').value['share of cases in those families']:.0f}% of cases are in families with two or more cases", f"shares the infectee's family in {pct1(pl.set_index('measure').value['documented pairs in the same family'])} of pairs against {pct1(pl.set_index('measure').value['all earlier-later pairs in the same family'])}", f"the school class in {pct0(pl.set_index('measure').value['documented pairs in the same class'])} against {pct0(pl.set_index('measure').value['all earlier-later pairs in the same class'])}", f"a median of {pl.set_index('measure').value['median distance, documented pairs (map units)']:.0f} map units apart, against {pl.set_index('measure').value['median distance, all earlier-later pairs (map units)']:.0f} for all pairs",
            f"by a median of {f_('rash_lag_median'):.0f} days (interquartile range {f_('rash_lag_q1'):.0f} to {f_('rash_lag_q3'):.0f})", f"moves the peak from day {int(f_('peak_day'))} to day {int(f_('peak_day_rash'))}", f"({f_('si_mean_prodrome'):.1f} days by prodrome, {f_('si_mean_rash'):.1f} by rash)", f"from {f_('si_sd_prodrome'):.1f} to {f_('si_sd_rash'):.1f} days", f"has a median of {f_('si_q50'):.0f} days, a 90th and 95th percentile of {f_('si_q95'):.0f} and a 99th of {f_('si_q99'):.0f}, and a maximum of {int(f_('si_max'))}" if rnd(f_('si_q90'), 0) == rnd(f_('si_q95'), 0) else "90th and 95th percentiles differ",
            f"A 7-day window finds only {pct0(tr_(7, 'anyone', 'recall'))} of documented infectors even when anyone is a candidate, a 10-day window {pct0(tr_(10, 'anyone', 'recall'))} and a 14-day window {pct0(tr_(14, 'anyone', 'recall'))}", f"anyone finds {pct0(tr_(14, 'anyone', 'recall'))} of infectors at {tr_(14, 'anyone', 'candidates_per_case'):.0f} candidates per case, family or class finds {pct0(tr_(14, 'family or class', 'recall'))} at {tr_(14, 'family or class', 'candidates_per_case'):.0f}, and family alone {pct0(tr_(14, 'same family', 'recall'))} at {tr_(14, 'same family', 'candidates_per_case'):.1f}", f"Within 25 map units finds {pct0(tr_(14, 'within 25 map units', 'recall'))} at {tr_(14, 'within 25 map units', 'candidates_per_case'):.0f} candidates, and within 100 map units {pct0(tr_(14, 'within 100 map units', 'recall'))} at {tr_(14, 'within 100 map units', 'candidates_per_case'):.0f}"]
    miss = [x for x in sent if x not in alltext_n]; check("D.8", "the sentences of Hands-on 3 that quote results (dates, peak, generations, age, recorded deaths, families, class, distances, onset definition, serial-interval quantiles, tracing windows and workload) match the tables and the summary", not miss, miss[:2])
    hh = 1 - 1 / np.array([2, 3.5, 8.1]); zz = [float(__import__("scipy.optimize", fromlist=["brentq"]).brentq(lambda z, r=r: z - 1 + math.exp(-r * z), 1e-6, 1)) for r in (2, 3.5)]
    check("D.9", "herd-immunity arithmetic: 1 − 1/R0 is 50%, 71% and 88% for R0 of 2, 3.5 and 8.1, as the text says (the 88% from one early estimate, flagged as arithmetic)", [rnd(100 * x, 0) for x in hh] == ["50", "71", "88"] and "needs 50% immunity, 3.5 needs 71%, and the early $R_t$ peak of 8.1" in md.replace("\\", "") or "50% immunity, 3.5 needs 71%" in alltext_n, [round(float(x), 3) for x in hh] + [round(z, 3) for z in zz])
    S = pd.read_csv(P("table_transphylo_summary.csv")).set_index("pi_true"); Sh = pd.read_csv(P("table_transphylo_short_summary.csv")).set_index("pi_true"); CBl = pd.read_csv(P("table_transphylo_calibration.csv"))
    t1 = [t for t in Tb if t.rows[0].cells[0].text == "Sampling probability" and "Recall" in " ".join(c.text for c in t.rows[0].cells)][0]; t2 = [t for t in Tb if t.rows[0].cells[0].text == "Sampling probability" and "ESS, 2,000" in t.rows[0].cells[1].text][0]
    r1 = lambda p_: [f"{p_:g}", f"{S.loc[p_].sampled_mean:.1f}", pct0(S.loc[p_].share_sampled_infector), str(int(S.loc[p_].links)), pct0(S.loc[p_].top1), f"{S.loc[p_].mean_prob_true:.2f}", pct0(S.loc[p_].recall_50), pct0(S.loc[p_].precision_50), f"{S.loc[p_].p_no_inf_when_sampled:.2f} / {S.loc[p_].p_no_inf_when_unsampled:.2f}", f"{S.loc[p_].ess_min:.0f} / {S.loc[p_].ess_median:.0f}", f"{S.loc[p_].pi_post:.2f}", f"{S.loc[p_].off_r_post:.1f}"]
    r2 = lambda p_: [f"{p_:g}", f"{Sh.loc[p_].ess_min:.1f} / {Sh.loc[p_].ess_median:.1f}", f"{S.loc[p_].ess_min:.0f} / {S.loc[p_].ess_median:.0f}", f"{Sh.loc[p_].precision_50:.2f}", f"{S.loc[p_].precision_50:.2f}", f"{Sh.loc[p_].recall_50:.2f}", f"{S.loc[p_].recall_50:.2f}", f"{Sh.loc[p_].pi_post:.2f}", f"{S.loc[p_].pi_post:.2f}"]
    check("D.10", "the two TransPhylo tables (accuracy and posterior by sampling probability, pooled over the six outbreaks; effective sample sizes and precision, recall and sampling probability for 2,000 against 100,000 iterations) equal the summary tables", cell(t1)[1:] == [r1(0.25), r1(0.75)] and cell(t2)[1:] == [r2(0.25), r2(0.75)], (cell(t1)[1:2], cell(t2)[1:2]))
    npairs = CBl.groupby("pi_true").n.sum(); runs_ = pd.read_csv(P("transphylo_runs.csv")); rs_ = pd.read_csv(P("transphylo_short_runs.csv"))
    sent2 = [f"only {pct0(S.loc[0.25].share_sampled_infector)} of sampled cases have a sampled true infector, against {pct0(S.loc[0.75].share_sampled_infector)} at 0.75", f"a mean probability of {S.loc[0.25].p_no_inf_when_unsampled:.2f} to \"no sampled infector\" at 0.25, against {S.loc[0.25].p_no_inf_when_sampled:.2f} for cases whose infector was sampled", f"Top-1 is {pct0(S.loc[0.25].top1)} at 0.25 (over {int(S.loc[0.25].links)} links) and {pct0(S.loc[0.75].top1)} at 0.75 (over {int(S.loc[0.75].links)}); calls above 0.5 have a precision of {pct0(S.loc[0.25].precision_50)} and {pct0(S.loc[0.75].precision_50)}", f"({int(npairs[0.25])} and {int(npairs[0.75])} pairs with probability at least 0.01)", f"Only {pct0(S.loc[0.25].recall_50)} and {pct0(S.loc[0.75].recall_50)} of true links reach 0.5", f"The posterior mean sampling probability is {S.loc[0.25].pi_post:.2f} against a true 0.25 and {S.loc[0.75].pi_post:.2f} against 0.75", f"effective sample size of the log posterior is {rs_.ess_log_pTTree.min():.1f} at worst and a median of {rs_.ess_log_pTTree.median():.1f} across the 12 outbreaks, against {runs_.ess_log_pTTree.min():.0f} at worst for 100,000", f"at 0.25 precision falls from {S.loc[0.25].precision_50:.2f} to {Sh.loc[0.25].precision_50:.2f}"]
    miss2 = [x for x in sent2 if x not in alltext_n]; check("D.11", "the sentences of the TransPhylo companion that quote results (shares with a sampled infector, probabilities of no sampled infector, top-1, precision, recall, calibration pairs, posterior sampling probability, effective sample sizes of long and short chains) match the tables", not miss2, miss2[:2])
    d_in = set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", md.split("## Sources")[0])); d_ref = set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", md.split("## Sources")[1]))
    check("D.12", "every DOI cited in the text is in the reference list, which holds 21 DOIs", d_in <= d_ref and len(d_ref) == 21, (len(d_in), len(d_ref), sorted(d_in - d_ref)[:3]))
    if A.doi:
        import urllib.request, urllib.parse, time
        badd = []
        for dd_ in sorted(d_ref | d_in):
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(dd_, safe="/()"), headers={"User-Agent": "validate-ch7"}), timeout=30).read()
            except Exception as ex: badd.append((dd_, str(ex)[:30]))
            time.sleep(0.15)
        check("DOI", f"all {len(d_ref | d_in)} DOIs of the chapter resolve in Crossref", not badd, badd[:4])
else: skip("D.*", "chapter text checks", "no --chapter given")
# ================================================================== E. reruns (needs --full)
if A.full:
    tmp = tempfile.mkdtemp(); src = os.path.dirname(os.path.abspath(P("measles_hagelloch_1861.RData"))); code = os.path.dirname(os.path.abspath(__file__))
    for f in glob.glob(os.path.join(src, "*")) + glob.glob(os.path.join(code, "*.py")) + glob.glob(os.path.join(code, "*.R")):
        if os.path.isfile(f): shutil.copy(f, tmp)
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); names = ["table_generations.csv", "table_person.csv", "table_place.csv", "table_tracing.csv", "table_descriptive_summary.csv", "table_transphylo_summary.csv", "table_transphylo_short_summary.csv", "table_transphylo_runs.csv", "table_transphylo_short_runs.csv", "table_transphylo_calibration.csv", "table_transphylo_short_calibration.csv"]; before = {n_: md5(os.path.join(tmp, n_)) for n_ in names}; [os.remove(os.path.join(tmp, n_)) for n_ in names]
    subprocess.run([sys.executable, "ex4_descriptive.py"], cwd=tmp, capture_output=True); subprocess.run([sys.executable, "ex5_transphylo_analysis.py"], cwd=tmp, capture_output=True)
    for n_ in names: check("REP", f"{n_} reproduced byte for byte", os.path.exists(os.path.join(tmp, n_)) and md5(os.path.join(tmp, n_)) == before[n_])
    if os.path.exists(os.path.join(tmp, "ebola.fasta")):
        old = md5(os.path.join(tmp, "ebola_candidate_ranking.csv")); subprocess.run([sys.executable, "ex2_ebola.py"], cwd=tmp, capture_output=True); check("REP", "Exercise 2 reproduces ebola_candidate_ranking.csv (as ranking.csv) byte for byte", os.path.exists(os.path.join(tmp, "ranking.csv")) and md5(os.path.join(tmp, "ranking.csv")) == old)
    r1 = subprocess.run([sys.executable, "ex1_hagelloch.py"], cwd=tmp, capture_output=True, text=True).stdout; check("REP", "Exercise 1 prints its documented numbers (serial interval 10.4 d, k = 0.25, R_t 8.1 on day 27, ranking 0.147 and 0.408)", all(x in r1 for x in ("mean 10.4 d", "k = 0.25", "Rt peaks at 8.1", "0.147", "0.408")))
    if A.sim:
        shutil.copy(os.path.join(tmp, "transphylo_runs.csv"), os.path.join(tmp, "old_transphylo_runs.csv")); [os.remove(f) for f in glob.glob(os.path.join(tmp, "transphylo*_run*.rds"))]
        subprocess.run(["Rscript", "ex5_transphylo.R"], cwd=tmp, capture_output=True); subprocess.run(["Rscript", "ex5_transphylo.R"], cwd=tmp, capture_output=True, env=dict(os.environ, TP_ITER="2000", TP_THIN="1", TP_PREFIX="transphylo_short"))
        a_, b_ = pd.read_csv(os.path.join(tmp, "old_transphylo_runs.csv")), pd.read_csv(os.path.join(tmp, "transphylo_runs.csv")); check("REP", "ex5_transphylo.R rerun from scratch (12 outbreaks, 100,000 iterations; about 40 minutes) reproduces the effective sample sizes and the posterior means of transphylo_runs.csv to six figures", len(a_) == len(b_) and all(np.allclose(a_[c], b_[c], rtol=1e-5) for c in ("ess_log_pTTree", "pi_post", "off_r_post", "neg_post", "sampled", "hosts")))
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch7.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=100)); npass = (Rd.result == "PASS").sum(); nskip = (Rd.result == "SKIP").sum(); print(f"\n{npass} of {len(Rd) - nskip} checks passed" + (f" ({nskip} skipped)" if nskip else "")); sys.exit(0 if (Rd.result == "FAIL").sum() == 0 else 1)
