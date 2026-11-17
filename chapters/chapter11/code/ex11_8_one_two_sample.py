"""
Chapter 11, Hands-on 8 - weak instruments, sample overlap and the winner's curse, one sample against two (individual-level simulation with a known truth).
Yao et al. (2025) conclude that with weak instruments a one-sample estimate is biased toward the confounded observational estimate and a two-sample estimate toward zero. This script checks that and adds the winner's curse.
A: 20 independent variants with equal effects on X, scaled to a target mean F; U confounds X and Y (0.6 each); the true effect of X on Y is 0.3; the observational slope is about 0.58 (confounded).
   Estimators: ordinary least squares (observational); one-sample 2SLS in 10,000 people; two-sample IVW with the exposure effects from one independent half of 10,000 and the outcome effects from the other; the same IVW with the two
   effects taken from the SAME sample (complete overlap).  Mean F of about 1.5, 5, 10, 30 and 100; 400 replicates.
B: winner's curse. 600 variants, each with a small effect on X drawn from a normal distribution (mean F about 5), 150 replicates; the instruments are the variants with P < 5e-8 for X.  Selected and estimated in the same sample, or selected in one sample and estimated in a second, independent one.
Output: Figure 7, table_one_two_sample.csv, table_winner.csv.   Needs numpy, pandas, scipy, matplotlib and mr11.py only.
"""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
rng = np.random.default_rng(12); BETA, N, REPS = 0.3, 10_000, 400
def gen(n, M, gam, maf):
    G = rng.binomial(2, maf, (n, M)).astype(float); G = (G - 2 * maf) / np.sqrt(2 * maf * (1 - maf)); U = rng.normal(size=n); X = G @ gam + 0.6 * U + rng.normal(size=n); Y = BETA * X + 0.6 * U + rng.normal(size=n); return G, X, Y
def gw(G, v): n = len(v); b = G.T @ (v - v.mean()) / n; se = np.sqrt((v.var() - b ** 2) / n); return b, se
def tsls(G, X, Y):
    Xh = G @ np.linalg.lstsq(G, X - X.mean(), rcond=None)[0]; return float(Xh @ (Y - Y.mean()) / (Xh @ Xh))
rows = []
for F in (1.5, 5, 10, 30, 100):
    M = 20; res = {k: [] for k in ("ols", "onesample", "twosample", "overlap")}; Fs = []
    for _ in range(REPS):
        maf = rng.uniform(0.1, 0.5, M); gam = np.sqrt(F / N) * np.ones(M) * rng.choice([-1, 1], M)
        G1, X1, Y1 = gen(N, M, gam, maf); G2, X2, Y2 = gen(N, M, gam, maf)
        res["ols"].append(np.cov(X1, Y1)[0, 1] / X1.var()); res["onesample"].append(tsls(G1, X1, Y1))
        bx, sx = gw(G1, X1); by, sy = gw(G2, Y2); s = np.sign(bx); res["twosample"].append(mr.ivw(bx * s, by * s, sy)["beta"]); Fs.append(np.mean(mr.fstat(bx, sx)))
        bx, sx = gw(G1, X1); by, sy = gw(G1, Y1); s = np.sign(bx); res["overlap"].append(mr.ivw(bx * s, by * s, sy)["beta"])
    rows.append(dict(target_F=F, mean_F=np.mean(Fs), truth=BETA, **{k: float(np.mean(v)) for k, v in res.items()}, **{f"{k}_sd": float(np.std(v)) for k, v in res.items()}, reps=REPS))
A = pd.DataFrame(rows); A.to_csv("table_one_two_sample.csv", index=False); print(A.round(3).to_string(index=False))
# ---------------------------------------------------------------- B: winner's curse
M, REPS_B, rows = 600, 150, []; thr = stats.norm.isf(5e-8 / 2); inflation = {"same sample": [], "independent sample": []}; nsel = []; true_sel = []
for _ in range(REPS_B):
    maf = rng.uniform(0.1, 0.5, M); gam = rng.normal(0, np.sqrt(5 / N), M)
    G1, X1, Y1 = gen(N, M, gam, maf); G2, X2, Y2 = gen(N, M, gam, maf); b1, s1 = gw(G1, X1); b2, s2 = gw(G2, X2); sel = np.abs(b1 / s1) > thr; nsel.append(sel.sum())
    if sel.sum() >= 3:
        true_sel.append(np.mean(np.abs(gam[sel]))); inflation["same sample"].append(np.mean(np.abs(b1[sel])) / np.mean(np.abs(gam[sel]))); inflation["independent sample"].append(np.mean(np.abs(b2[sel])) / np.mean(np.abs(gam[sel])))
B = pd.DataFrame([dict(estimated_in=k, variants_selected_mean=float(np.mean(nsel)), mean_abs_estimate_over_mean_abs_truth=float(np.mean(v)), replicates=len(v)) for k, v in inflation.items()]); B.to_csv("table_winner.csv", index=False); print(B.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3)); x = np.arange(len(A)); cols = {"ols": "#a0aec0", "onesample": "#E4572E", "twosample": "#1B998B", "overlap": "#F3A712"}; lab = {"ols": "Observational (OLS)", "onesample": "One-sample 2SLS", "twosample": "Two-sample IVW", "overlap": "IVW, same sample twice"}
for k in cols: ax[0].plot(x, A[k], "o-", color=cols[k], lw=2, label=lab[k])
ax[0].axhline(BETA, color="#2E294E", lw=1.2, ls="--"); ax[0].text(len(A) - 1.02, BETA + 0.015, "truth 0.3", ha="right", fontsize=9, color="#2E294E"); ax[0].set_xticks(x); ax[0].set_xticklabels([f"{f:.0f}" for f in A.mean_F]); ax[0].set(xlabel="Mean F-statistic of the instruments", ylabel="Mean estimate of the effect of X on Y (400 replicates)", title="Weak instruments pull in two directions"); ax[0].legend(frameon=False, fontsize=8.5, loc="lower right")
ax[1].bar([0, 1], B.mean_abs_estimate_over_mean_abs_truth, color=["#E4572E", "#1B998B"], width=0.55); ax[1].axhline(1, color="#2E294E", lw=1.2, ls="--"); ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["Selected and estimated\nin the same sample", "Selected in one sample,\nestimated in another"]); ax[1].set(ylabel="Estimated / true effect size (selected variants)", title="Winner's curse"); 
for i, v in enumerate(B.mean_abs_estimate_over_mean_abs_truth): ax[1].text(i, v + 0.03, f"{v:.2f}", ha="center")
fig.tight_layout(); fig.savefig("figs11/fig7_one_two_sample.png", dpi=200); print("figure written")
