"""
Chapter 11, Hands-on 3 - A known truth: how weak instruments, sample overlap and pleiotropy bias Mendelian randomisation estimates.
Summary-level simulation: 100 variants, a true causal effect of 0.3, 500 replicates per scenario, three estimators (IVW random effects, MR-Egger, weighted median).
Output: table_simulation.csv and Figure 3.   Needs numpy pandas scipy matplotlib and mr11.py.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.0, "axes.spines.top": False, "axes.spines.right": False})
BETA, M, N, REPS, NBOOT = 0.3, 100, 100_000, 500, 100

def simulate(rng, F, alpha_mu=0.0, alpha_sd=0.0, overlap=0.0, invalid_frac=0.0):
    maf = rng.uniform(0.05, 0.5, M); sx = 1 / np.sqrt(2 * N * maf * (1 - maf)); sy = sx.copy()                      # standard errors from the sample size and frequency
    gamma = np.sqrt(F) * sx * rng.normal(size=M)                                                                  # variant effects on the exposure; mean F is the target
    alpha = rng.normal(alpha_mu, alpha_sd, M) * np.sign(gamma)                                                    # pleiotropic effects, expressed relative to the allele that raises the exposure
    if invalid_frac > 0:                                                                                           # a few invalid variants that are also the strongest: InSIDE fails
        bad = rng.random(M) < invalid_frac; gamma = np.where(bad, 2 * gamma, gamma); alpha = np.where(bad, 0.03 * np.sign(gamma), alpha)
    ex = rng.normal(size=M); ey = overlap * ex + np.sqrt(1 - overlap ** 2) * rng.normal(size=M)                    # correlated errors when the two samples overlap
    bx = gamma + sx * ex; by = BETA * gamma + alpha + sy * ey
    return bx, by, sy, sx

SC = [("1 Valid, strong (F about 60)", dict(F=60)), ("2 Weak (F about 4), no overlap", dict(F=4)), ("3 Weak, full sample overlap", dict(F=4, overlap=0.5)),
      ("4 Balanced pleiotropy", dict(F=60, alpha_sd=0.01)), ("5 Directional pleiotropy", dict(F=60, alpha_mu=0.01, alpha_sd=0.01)), ("6 Invalid and strongest (15%)", dict(F=60, invalid_frac=0.15))]
rng = np.random.default_rng(11); rows, store = [], {}
for name, kw in SC:
    est = {"IVW": [], "MR-Egger": [], "Weighted median": []}; cov = {k: [] for k in est}; Fs = []
    for _ in range(REPS):
        bx, by, sy, sx = simulate(rng, **kw); s = np.sign(bx); Fs.append(np.mean(mr.fstat(bx, sx)))
        I = mr.ivw(bx, by, sy); E = mr.egger(bx, by, sy); W = mr.weighted_median(by / bx, sy / np.abs(bx), nboot=NBOOT, seed=int(rng.integers(1e9)))
        for k, (b, se) in {"IVW": (I["beta"], I["se_random"]), "MR-Egger": (E["slope"], E["se_slope"]), "Weighted median": (W["beta"], W["se"])}.items():
            est[k].append(b); cov[k].append(abs(b - BETA) < 1.96 * se)
    store[name] = est
    for k in est: e = np.array(est[k]); rows.append(dict(scenario=name, method=k, mean_F=round(float(np.mean(Fs)), 1), mean_estimate=round(float(e.mean()), 3), bias=round(float(e.mean() - BETA), 3), sd=round(float(e.std(ddof=1)), 3), coverage=round(float(np.mean(cov[k])), 3)))
T = pd.DataFrame(rows); T.to_csv("table_simulation.csv", index=False); print(T.to_string(index=False))

fig, ax = plt.subplots(2, 3, figsize=(12.2, 6.6), sharey=True); cols = {"IVW": "#2b6cb0", "MR-Egger": "#c53030", "Weighted median": "#2f855a"}
for a, (name, _) in zip(ax.ravel(), SC):
    bp = a.boxplot([store[name][k] for k in cols], tick_labels=["IVW", "Egger", "W. median"], patch_artist=True, showfliers=False, widths=0.6)
    for b, c in zip(bp["boxes"], cols.values()): b.set(facecolor=c, alpha=0.55)
    for m in bp["medians"]: m.set(color="black"); 
    a.axhline(BETA, color="black", ls="--", lw=1); a.set_title(name, fontsize=9.5)
ax[0, 0].set_ylabel("Estimate (true effect 0.3, dashed)"); ax[1, 0].set_ylabel("Estimate (true effect 0.3, dashed)"); ax[0, 0].set_ylim(-0.35, 0.95); fig.tight_layout(); fig.savefig("figs11/fig3_simulation.png", dpi=200); plt.close(fig); print("figure written")
