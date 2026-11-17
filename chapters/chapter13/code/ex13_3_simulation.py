"""
Chapter 13, Hands-on 3 - A known truth: what do genomes add to case data?
A branching-process epidemic (R = 1.5, gamma generation interval, mean 8 d, SD 5 d) is simulated with its transmission tree. Cases are reported with a probability that may fall over time,
and genomes are sampled from infected people. The growth rate is estimated from reported cases (quasi-Poisson regression), from the genealogy of the genomes (exponential-growth coalescent),
and by combining the two with inverse-variance weights. 'Resolution' rounds coalescent times to multiples of tau days to mimic genomes with too few mutations to resolve the tree.
Output: table_simulation13.csv and Figure 3.
"""
import numpy as np, pandas as pd
from scipy import stats, optimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ebolacore import gamma_si, poisson_growth
from coalfit13 import fit_r, se_r
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
GI_MEAN, GI_SD, R0, T, I0, REPS = 8.0, 5.0, 1.5, 90, 10, 120; W = gamma_si(GI_MEAN, GI_SD, 40)
k_, th_ = (GI_MEAN / GI_SD) ** 2, GI_SD ** 2 / GI_MEAN; R_TRUE_R = (R0 ** (1 / k_) - 1) / th_                      # the growth rate that goes with R0 and the generation interval
rng = np.random.default_rng(13)

def simulate():
    inf_day = [0] * I0; parent = [-1] * I0; by_day = {0: list(range(I0))}
    for t in range(1, T + 1):
        lam = R0 * sum(W[s - 1] * len(by_day.get(t - s, ())) for s in range(1, 41)); n = rng.poisson(lam)
        if n == 0: by_day[t] = []; continue
        pw = np.array([W[s - 1] * len(by_day.get(t - s, ())) for s in range(1, 41)]); pw = pw / pw.sum(); lags = rng.multinomial(n, pw); new = []
        for s, c in enumerate(lags, 1):
            src = by_day.get(t - s, [])
            for _ in range(c): new.append(len(inf_day)); inf_day.append(t); parent.append(src[rng.integers(len(src))])
        by_day[t] = new
        if len(inf_day) > 6000: break
    return np.array(inf_day), np.array(parent)

def genealogy(inf_day, parent, tips, tau):
    """Coalescent events of the sampled genomes from the transmission tree; times rounded to multiples of tau days to remove resolution."""
    tipset = {int(i): int(inf_day[i]) + 5 for i in tips}; arrivals = {}                                           # a tip is sampled 5 days after infection
    for i in tips: arrivals.setdefault(int(i), []).append(tipset[int(i)])
    order = np.argsort(-inf_day, kind="stable"); events = []
    for x in order:
        a = arrivals.get(int(x))
        if not a: continue
        a.sort(reverse=True); events += a[1:]                                                                      # all but the most recent arrival coalesce at their own times
        if parent[x] >= 0: arrivals.setdefault(int(parent[x]), []).append(int(inf_day[x]))
    ev_c = np.array(events, float)
    if tau > 0: ev_c = np.round(ev_c / tau) * tau
    tmax = max(tipset.values()); ev = [(tmax - t, +1) for t in tipset.values()] + [(tmax - t, -1) for t in ev_c]; ev.sort(key=lambda e: (e[0], -e[1])); return ev

def one(report, n_gen, tau):
    inf_day, parent = simulate(); N = len(inf_day)
    if N < 200: return None
    rho = report(inf_day.max()); daily = np.bincount(inf_day, minlength=T + 1)[: inf_day.max() + 1]; probs = np.array([report(t) for t in range(len(daily))]); obs = rng.binomial(daily, probs)
    lo = 10; rc, sec = poisson_growth(obs[lo:])
    tips = rng.choice(N, size=min(n_gen, N), replace=False); ev = genealogy(inf_day, parent, tips, tau)
    try: rg, lnN = fit_r(ev); sg = se_r(ev, rg, lnN)
    except Exception: return None
    w1, w2 = 1 / sec ** 2, 1 / sg ** 2 if np.isfinite(sg) and sg > 0 else 0.0; rcomb = (rc * w1 + rg * w2) / (w1 + w2); scomb = 1 / np.sqrt(w1 + w2)
    return dict(r_cases=rc, se_cases=sec, r_gen=rg, se_gen=sg, r_comb=rcomb, se_comb=scomb, n_cases=N)

SC = [("1 Full reporting, resolved genomes", lambda t: 1.0, 60, 0), ("2 Reporting falls 80% to 30%", lambda t: 0.8 - 0.5 * min(t, 90) / 90, 60, 0), ("3 As 2, genomes unresolved (14 d)", lambda t: 0.8 - 0.5 * min(t, 90) / 90, 60, 14),
      ("4 Constant 30%, only 15 genomes", lambda t: 0.3, 15, 0), ("5 Full reporting, 60 genomes, unresolved (14 d)", lambda t: 1.0, 60, 14)]
rows, store = [], {}
for name, rep, ng, tau in SC:
    res = [x for x in (one(rep, ng, tau) for _ in range(REPS)) if x]; D = pd.DataFrame(res); store[name] = D
    for src, c, s in [("Cases", "r_cases", "se_cases"), ("Genomes", "r_gen", "se_gen"), ("Combined", "r_comb", "se_comb")]:
        cov = (np.abs(D[c] - R_TRUE_R) < 1.96 * D[s]).mean(); rows.append(dict(scenario=name, source=src, replicates=len(D), mean_cases=D.n_cases.mean(), true_r=R_TRUE_R, mean_r=D[c].mean(), bias=D[c].mean() - R_TRUE_R, sd=D[c].std(), coverage=cov))
T_ = pd.DataFrame(rows); T_.round(4).to_csv("table_simulation13.csv", index=False); print(f"true growth rate r = {R_TRUE_R:.4f} per day (R = {R0}); doubling {np.log(2) / R_TRUE_R:.1f} d"); print(T_.round(3).to_string(index=False))
fl = []
for name, D in store.items(): z_ = (D.r_cases - D.r_gen) / np.hypot(D.se_cases, D.se_gen); fl.append(dict(scenario=name, disagreement_flagged_pct=(np.abs(z_) > 1.96).mean() * 100, mean_se_cases=D.se_cases.mean(), sd_r_cases=D.r_cases.std(), mean_se_gen=D.se_gen.mean(), sd_r_gen=D.r_gen.std()))
F = pd.DataFrame(fl); F.round(4).to_csv("table_simulation13_flags.csv", index=False); print(F.round(4).to_string(index=False))
fig, ax = plt.subplots(1, 5, figsize=(14.2, 4.4), sharey=True); cols = {"Cases": "#dd6b20", "Genomes": "#2b6cb0", "Combined": "#2f855a"}
for a, (name, *_ ) in zip(ax, SC):
    D = store[name]; bp = a.boxplot([D.r_cases.dropna(), D.r_gen.dropna(), D.r_comb.dropna()], tick_labels=["Cases", "Genomes", "Both"], patch_artist=True, showfliers=False, widths=0.6)
    for b, c in zip(bp["boxes"], cols.values()): b.set(facecolor=c, alpha=0.55)
    for m_ in bp["medians"]: m_.set(color="black")
    a.axhline(R_TRUE_R, color="black", ls="--", lw=1); a.set_title(name, fontsize=8.6); a.tick_params(axis="x", labelsize=8.5)
ax[0].set_ylabel("Estimated growth rate (per day; true value dashed)"); fig.tight_layout(); fig.savefig("figs13/fig3_simulation.png", dpi=200); plt.close(fig); print("figure written")
