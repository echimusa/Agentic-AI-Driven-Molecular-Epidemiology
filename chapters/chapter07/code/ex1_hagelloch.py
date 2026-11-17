"""
Chapter 7, Exercise 1 - Who infected whom in the Hagelloch measles outbreak (1861).

Public data: 188 cases in a German village with the documented infector of 184 of
them (Pfeilsticker 1863; Oesterle 1992), distributed in the CRAN package `outbreaks`.
Everything below runs in plain Python: no R needed.
"""
import tarfile, tempfile, urllib.request, os
import numpy as np
import pandas as pd
import pyreadr
from scipy import stats, optimize
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter

URL = "https://cran.r-project.org/src/contrib/outbreaks_1.9.0.tar.gz"
MEMBER = "outbreaks/data/measles_hagelloch_1861.RData"
OUT = "figs"; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)

# ---------------------------------------------------------------- 1. load
def load_hagelloch(local="measles_hagelloch_1861.RData"):
    if os.path.exists(local):                         # supplied in the data/ folder of the chapter zip
        return pyreadr.read_r(local)["measles_hagelloch_1861"]
    tmp = tempfile.mkdtemp()
    tgz = os.path.join(tmp, "outbreaks.tar.gz")
    for url in (URL, URL.replace("/src/contrib/", "/src/contrib/Archive/outbreaks/")):   # CRAN moves old versions to Archive/
        try: urllib.request.urlretrieve(url, tgz); break
        except Exception: continue
    with tarfile.open(tgz) as tar:
        tar.extract(MEMBER, tmp)
    return pyreadr.read_r(os.path.join(tmp, MEMBER))["measles_hagelloch_1861"]

d = load_hagelloch()
d["t"] = (pd.to_datetime(d.date_of_prodrome) - pd.to_datetime(d.date_of_prodrome).min()).dt.days
d["infector"] = pd.to_numeric(d.infector)
d = d.set_index("case_ID", drop=False)
print(f"{len(d)} cases, onset day 0-{d.t.max()}, documented infector for {d.infector.notna().sum()}")

# ------------------------------------------- 2. serial interval from known pairs
pairs = d.dropna(subset=["infector"]).copy()
pairs["infector"] = pairs.infector.astype(int)
pairs["si"] = pairs.t - d.t.loc[pairs.infector].values
n_bad = int((pairs.si < 1).sum())
si = pairs.si[pairs.si >= 1]
shape, _, scale = stats.gamma.fit(si, floc=0)
print(f"serial interval (prodrome to prodrome): n={len(si)} pairs, mean {si.mean():.1f} d, "
      f"median {si.median():.0f} d, gamma shape {shape:.1f}, scale {scale:.2f}; "
      f"{n_bad} pairs with interval < 1 d set aside")

def si_pmf(shape, scale, smax=60):
    s = np.arange(0, smax + 1)
    p = stats.gamma.cdf(s + 0.5, shape, scale=scale) - stats.gamma.cdf(s - 0.5, shape, scale=scale)
    p[0] = 0.0
    return p / p.sum()

w = si_pmf(shape, scale)

# ------------------------------------------------ 3. offspring distribution
counts = pairs.infector.value_counts().reindex(d.case_ID, fill_value=0)
m = counts.mean()
def nll(logk):
    k = np.exp(logk); p = k / (k + m)
    return -stats.nbinom.logpmf(counts, k, p).sum()
k_hat = float(np.exp(optimize.minimize_scalar(nll, bounds=(-4, 4), method="bounded").x))
top = np.sort(counts.values)[::-1]
share10 = top[: int(round(0.10 * len(top)))].sum() / top.sum()
print(f"offspring: mean {m:.2f}, max {counts.max()} (case {counts.idxmax()}), "
      f"NB dispersion k = {k_hat:.2f}, top 10% of cases cause {share10:.0%} of infections, "
      f"{(counts == 0).mean():.0%} infect nobody")

# --------------------------------------------- 4. Rt (Cori et al. 2013, closed form)
I = np.bincount(d.t.values, minlength=d.t.max() + 1).astype(float)
Lam = np.convolve(I, w)[: len(I)]
tau, a0, b0 = 7, 1.0, 5.0                      # 7-day window, Gamma(1, scale 5) prior
rt = []
for t in range(tau, len(I)):
    sl = slice(t - tau + 1, t + 1)
    a, rate = a0 + I[sl].sum(), 1 / b0 + Lam[sl].sum()
    rt.append((t, a / rate, *stats.gamma.ppf([0.025, 0.975], a, scale=1 / rate)))
rt = pd.DataFrame(rt, columns=["t", "mean", "lo", "hi"])
rt = rt[[I[max(0, t - tau + 1): t + 1].sum() >= 5 for t in rt.t]]   # need >= 5 cases in window
peak = rt.loc[rt["mean"].idxmax()]
print(f"Rt peaks at {peak['mean']:.1f} (95% CrI {peak.lo:.1f}-{peak.hi:.1f}) on day {int(peak.t)}; "
      f"Rt < 1 from day {int(rt[(rt['mean'] < 1) & (rt.t > peak.t)].t.min())}")

# ------------------------------- 5. rank candidate infectors (5-fold cross-validated)
xy = d[["x_loc", "y_loc"]].to_numpy(); ids = d.case_ID.to_numpy(); tt = d.t.to_numpy()
idx = {c: i for i, c in enumerate(ids)}
D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])

def cand_weights(j, w, lam):
    dt = tt[j] - tt                                  # infectee minus candidate infector
    ok = (dt >= 1) & (dt < len(w))
    wt = np.zeros(len(tt)); wt[ok] = w[dt[ok]]
    if np.isfinite(lam): wt *= np.exp(-D[j] / lam)
    return wt

def loglik(js, w, lam):
    ll = 0.0
    for j in js:
        wt = cand_weights(j, w, lam); true = idx[int(d.infector.iloc[j])]
        if wt[true] > 0 and wt.sum() > 0: ll += np.log(wt[true] / wt.sum())
    return ll

known = np.array([j for j in range(len(d)) if not np.isnan(d.infector.iloc[j])])
folds = np.array_split(rng.permutation(known), 5)
grid = [2, 5, 10, 20, 50, 100, 200, np.inf]
res = {"time only": [], "time + distance": []}; base = []; lams = []
for f in folds:
    train = np.setdiff1d(known, f)
    s_tr = pairs.si.loc[ids[train]]; s_tr = s_tr[s_tr >= 1]
    sh, _, sc = stats.gamma.fit(s_tr, floc=0); w_tr = si_pmf(sh, sc)
    lam = max(grid, key=lambda L: loglik(train, w_tr, L))
    lams.append(lam)
    for name, L in [("time only", np.inf), ("time + distance", lam)]:
        for j in f:
            wt = cand_weights(j, w_tr, L); true = idx[int(d.infector.iloc[j])]
            p = wt / wt.sum() if wt.sum() > 0 else np.full(len(wt), 1 / len(wt))
            res[name].append((p.argmax() == true, p[true]))
    for j in f:
        base.append(1 / max(1, int(((tt[j] - tt) >= 1).sum())))
summary = pd.DataFrame({k: {"top-1 correct": np.mean([r[0] for r in v]),
                            "mean probability on true infector": np.mean([r[1] for r in v])}
                        for k, v in res.items()}).T
summary.loc["random among earlier cases"] = [np.mean(base), np.mean(base)]
print(summary.round(3)); print("distance scale chosen per fold:", lams)

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
C0, C1, C2 = "#2b6cb0", "#c05621", "#718096"

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].bar(np.arange(len(I)), I, color=C0, width=0.9)
ax[0].set(xlabel="Day since first prodrome (30 Oct 1861)", ylabel="New cases per day",
          title="Epidemic curve, 188 cases")
ax[0].yaxis.set_major_locator(MaxNLocator(integer=True))
bins = np.arange(0.5, 32.5, 1)
ax[1].hist(si, bins=bins, color=C0, alpha=0.8, density=True)
xs = np.linspace(1, 32, 200)
ax[1].plot(xs, stats.gamma.pdf(xs, shape, scale=scale), color=C1, lw=2)
ax[1].axvline(si.mean(), color=C2, ls="--")
ax[1].set(xlabel="Serial interval (days)", ylabel="Density",
          title=f"Serial interval, {len(si)} known pairs")
ax[1].text(si.mean() + 0.5, ax[1].get_ylim()[1] * 0.9, f"mean {si.mean():.1f} d", color=C2)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_epicurve_si.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
vals = np.arange(0, counts.max() + 1)
ax[0].bar(vals, [(counts == v).mean() for v in vals], color=C0, alpha=0.8, label="observed")
kk = k_hat; ax[0].plot(vals, stats.nbinom.pmf(vals, kk, kk / (kk + m)), "o-", color=C1, label=f"negative binomial, k = {k_hat:.2f}")
ax[0].set(xlabel="Secondary cases per case", ylabel="Proportion of cases", title="Offspring distribution")
ax[0].set_yscale("log"); ax[0].legend(frameon=False)
cum = np.cumsum(top) / top.sum(); frac = np.arange(1, len(top) + 1) / len(top)
ax[1].plot(np.r_[0, frac], np.r_[0, cum], color=C0, lw=2)
ax[1].plot([0, 1], [0, 1], color=C2, ls="--")
ax[1].set(xlabel="Fraction of cases, most infectious first", ylabel="Fraction of all transmissions",
          title=f"Top 10% of cases cause {share10:.0%} of infections")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_offspring.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4))
ax.fill_between(rt.t, rt.lo, rt.hi, color=C0, alpha=0.2)
ax.plot(rt.t, rt["mean"], color=C0, lw=2); ax.axhline(1, color=C2, ls="--")
ax.set(xlabel="Day since first prodrome", ylabel="Rt (7-day window)", yscale="log",
       title="Time-varying reproduction number (mean and 95% credible interval)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_rt.png", dpi=200); plt.close(fig)

# transmission forest: x = onset day, y = depth-first leaf order, colour = school class
children = pairs.groupby("infector").case_ID.apply(list).to_dict()
roots = d.index[d.infector.isna()].tolist()
ypos, y = {}, [0]
def place(c, seen):
    seen.add(c); kids = [k for k in children.get(c, []) if k not in seen]
    if not kids: ypos[c] = y[0]; y[0] += 1; return
    for k in kids: place(k, seen)
    ypos[c] = np.mean([ypos[k] for k in kids])
seen = set()
for r in roots: place(r, seen)
for c in d.index:
    if c not in ypos: place(c, seen)
cols = {"0": "#718096", "1": "#c05621", "2": "#2b6cb0"}
fig, ax = plt.subplots(figsize=(11, 6.5))
for c, inf in pairs.infector.items():
    ax.plot([d.t[inf], d.t[c]], [ypos[inf], ypos[c]], color="#cbd5e0", lw=0.7, zorder=1)
ax.scatter(d.t, [ypos[c] for c in d.index], s=14 + 14 * counts.loc[d.index].values,
           c=[cols[str(k)] for k in d["class"]], zorder=2, edgecolor="white", linewidth=0.4)
big = counts.idxmax(); ax.annotate(f"case {big}: {counts.max()} secondary cases", (d.t[big], ypos[big]),
            xytext=(d.t[big] - 3, ypos[big] - 14), ha="right", arrowprops=dict(arrowstyle="-", color=C2))
ax.set(xlabel="Day since first prodrome", yticks=[], xticks=range(0, 90, 20), xlim=(-26, 90),
       title="Documented transmission forest, Hagelloch 1861 (dot size = secondary cases)")
ax.spines["left"].set_visible(False)
for k, lab in [("0", "pre-school"), ("1", "first class"), ("2", "second class")]:
    ax.scatter([], [], c=cols[k], label=lab, s=30)
ax.legend(frameon=False, title="School group", loc="center right")
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_transmission_forest.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4))
names = list(summary.index); xs = np.arange(len(names))
ax.bar(xs - 0.2, summary["top-1 correct"], 0.4, color=C0, label="top-ranked candidate is the infector")
ax.bar(xs + 0.2, summary["mean probability on true infector"], 0.4, color=C1, label="mean probability on the infector")
ax.set_xticks(xs); ax.set_xticklabels(names); ax.set_ylim(0, 0.6); ax.yaxis.set_major_formatter(PercentFormatter(1.0))
for i, n in enumerate(names):
    ax.text(i - 0.2, summary["top-1 correct"].iloc[i] + 0.02, f"{summary['top-1 correct'].iloc[i]:.0%}", ha="center")
    ax.text(i + 0.2, summary["mean probability on true infector"].iloc[i] + 0.02,
            f"{summary['mean probability on true infector'].iloc[i]:.0%}", ha="center")
ax.set(title="Five-fold cross-validation against the documented infectors"); ax.legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_ranking.png", dpi=200); plt.close(fig)
print("figures written to", OUT)
