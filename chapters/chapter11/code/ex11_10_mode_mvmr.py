"""
Chapter 11, Hands-on 10 - invalid instruments revisited: the weighted mode and the plurality rule, and multivariable MR.
A  Summary-level simulation (100 variants, 100,000 people in each of two studies, 300 replicates, true effect 0.3, mean F about 60) with three patterns of invalid instruments:
     P1  all instruments valid (control)
     P2  plurality valid: 35% of the variants are valid and the other 65% have Wald ratios shifted by an amount that is scattered (mean 0.15, standard deviation 0.4), so the valid ones form the largest cluster of ratios but not a majority
     P3  plurality invalid: 40% valid, and the other 60% share the same shift of 0.2 in the ratio (a second causal pathway), so the largest cluster of ratios is the wrong one
   Estimators: IVW, MR-Egger, weighted median, weighted mode.  Reported: mean estimate, standard deviation and root mean squared error.
B  The same weighted mode on the 893 body mass index instruments of Hands-on 1 (the harmonised set).
C  Multivariable MR. Two exposures whose effects across variants are correlated (0.7); the outcome depends on X1 with coefficient 0.0 and on X2 with 0.4 (case C1), or on 0.3 and 0.4 (case C2). 150 variants, 500 replicates.
   Univariable IVW for X1 against multivariable IVW for X1 and X2 together.
Output: Figure 9, table_mode.csv, table_mode_bmi.csv, table_mvmr.csv.   Needs numpy, pandas, scipy, matplotlib and mr11.py; the cached BMI inputs for part B.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
rng = np.random.default_rng(14); N = 100_000; s = 1 / np.sqrt(N); BETA = 0.3
# ---------------------------------------------------------------- A
def pattern(name, M=100, F=60):
    gam = np.sqrt(F) * s * rng.normal(size=M); alpha = np.zeros(M)
    if name == "P2": k = int(0.35 * M); alpha[k:] = rng.normal(0.15, 0.4, M - k) * gam[k:]      # pleiotropic effect = shift in the ratio estimate times the exposure effect, so shifts are relative to the exposure-raising allele
    if name == "P3": k = int(0.40 * M); alpha[k:] = 0.2 * gam[k:]
    bx = gam + s * rng.normal(size=M); by = BETA * gam + alpha + s * rng.normal(size=M); sg = np.sign(bx); return bx * sg, by * sg, np.full(M, s)
rows = []
for name, label in (("P1", "P1 All valid"), ("P2", "P2 Plurality valid (35%), invalid scattered"), ("P3", "P3 Plurality invalid (60% share one pleiotropic effect)")):
    est = {"IVW": [], "MR-Egger": [], "Weighted median": [], "Weighted mode": []}
    for _ in range(300):
        bx, by, sy = pattern(name); bw, sw = mr.wald(bx, by, sy); est["IVW"].append(mr.ivw(bx, by, sy)["beta"]); est["MR-Egger"].append(mr.egger(bx, by, sy)["slope"]); est["Weighted median"].append(mr.weighted_median(bw, sw, nboot=0 + 1)["beta"]); est["Weighted mode"].append(mr.weighted_mode(bw, sw)["beta"])
    for k, v in est.items(): v = np.array(v); rows.append(dict(pattern=label, estimator=k, mean=v.mean(), sd=v.std(), rmse=float(np.sqrt(np.mean((v - BETA) ** 2))), truth=BETA))
T = pd.DataFrame(rows); T.to_csv("table_mode.csv", index=False); print(T.round(3).to_string(index=False))
# ---------------------------------------------------------------- B: the mode on the body mass index instruments
X = pd.read_csv("bmi_exposure_snps.csv").rename(columns={"Allele1": "ea", "Allele2": "oa", "Freq1": "eaf_x", "Effect": "bx", "StdErr": "sx", "P-value": "px"}); X[["ea", "oa"]] = X[["ea", "oa"]].apply(lambda c: c.str.upper()); Y = pd.read_csv("hgi_for_bmi.csv")
D = X.drop_duplicates("rsid").merge(Y.drop_duplicates("rsid"), on="rsid"); keep, sign = mr.harmonise(D.ea.values, D.oa.values, D.eaf_x.values, D.oa1.values, D.oa2.values, D.eaf_y.values); D = D[keep].copy(); D["by"] = D["by"].values * sign[keep]
bx, by, sy = D.bx.values, D.by.values, D.sy.values; bw, sw = mr.wald(bx, by, sy); I = mr.ivw(bx, by, sy); Wm = mr.weighted_median(bw, sw, nboot=1000); Mo = mr.weighted_mode(bw, sw, nboot=500); Mo5 = mr.weighted_mode(bw, sw, phi=0.5)
Bm = pd.DataFrame([dict(method="IVW, random effects", beta=I["beta"], se=I["se_random"]), dict(method="Weighted median", beta=Wm["beta"], se=Wm["se"]), dict(method="Weighted mode (phi = 1)", beta=Mo["beta"], se=Mo["se"]), dict(method="Weighted mode (phi = 0.5), point estimate only", beta=Mo5["beta"], se=np.nan)])
Bm["odds_ratio"] = np.exp(Bm.beta); Bm["ci_low"] = np.exp(Bm.beta - 1.96 * Bm.se); Bm["ci_high"] = np.exp(Bm.beta + 1.96 * Bm.se); Bm["instruments"] = len(D); Bm.to_csv("table_mode_bmi.csv", index=False); print(Bm.round(3).to_string(index=False))
# ---------------------------------------------------------------- C: multivariable MR
rows = []; M = 150
for case, (b1, b2) in (("C1 X1 has no effect, X2 does", (0.0, 0.4)), ("C2 both have effects", (0.3, 0.4))):
    u1, mv1, mv2, c1, c2 = [], [], [], [], []
    for _ in range(500):
        z1, z2 = rng.normal(size=M), rng.normal(size=M); g1 = np.sqrt(50) * s * z1; g2 = np.sqrt(50) * s * (0.7 * z1 + np.sqrt(1 - 0.49) * z2)
        bx1, bx2 = g1 + s * rng.normal(size=M), g2 + s * rng.normal(size=M); by = b1 * g1 + b2 * g2 + s * rng.normal(size=M); sg = np.sign(bx1); sy = np.full(M, s)
        u1.append(mr.ivw(bx1 * sg, by * sg, sy)["beta"]); V = mr.mvmr_ivw(np.column_stack([bx1 * sg, bx2 * sg]), by * sg, sy); mv1.append(V["beta"][0]); mv2.append(V["beta"][1]); c1.append(abs(V["beta"][0] - b1) < 1.96 * V["se"][0])
    rows.append(dict(case=case, true_x1=b1, true_x2=b2, univariable_x1=np.mean(u1), multivariable_x1=np.mean(mv1), multivariable_x2=np.mean(mv2), multivariable_x1_covers_truth=float(np.mean(c1)), reps=500))
Cm = pd.DataFrame(rows); Cm.to_csv("table_mvmr.csv", index=False); print(Cm.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.2)); pats = list(dict.fromkeys(T.pattern)); est = ["IVW", "MR-Egger", "Weighted median", "Weighted mode"]; cols = ["#a0aec0", "#1B998B", "#F3A712", "#E4572E"]; w = 0.2
for i, e in enumerate(est): sub = T[T.estimator == e]; ax[0].bar(np.arange(3) + (i - 1.5) * w, sub["mean"], w, yerr=sub.sd, color=cols[i], label=e, error_kw=dict(ecolor="#2E294E", elinewidth=0.9, capsize=2))
ax[0].axhline(BETA, color="#2E294E", lw=1.2, ls="--"); ax[0].set_xticks(range(3)); ax[0].set_xticklabels(["P1", "P2", "P3"]); ax[0].set(ylabel="Mean estimate ± SD (truth 0.3)", title="Mode and median under invalid instruments"); ax[0].legend(frameon=False, fontsize=8, loc="upper left")
k = np.argsort(bw); ax[1].hist(bw, bins=np.linspace(-3, 3, 60), color="#3D7EA6", alpha=0.85); ax[1].axvline(Mo["beta"], color="#E4572E", lw=2, label=f"weighted mode {Mo['beta']:.2f}"); ax[1].axvline(I["beta"], color="#2E294E", lw=1.6, ls="--", label=f"IVW {I['beta']:.2f}"); ax[1].set(xlabel="Wald ratio, BMI on hospitalised COVID-19 (log odds ratio)", ylabel="Instruments", title="893 body mass index instruments"); ax[1].legend(frameon=False, fontsize=8.5)
x = np.arange(2); ax[2].bar(x - 0.27, Cm.univariable_x1, 0.27, color="#a0aec0", label="Univariable IVW, X1"); ax[2].bar(x, Cm.multivariable_x1, 0.27, color="#E4572E", label="Multivariable IVW, X1"); ax[2].bar(x + 0.27, Cm.multivariable_x2, 0.27, color="#1B998B", label="Multivariable IVW, X2")
for i, (a, b) in enumerate(zip(Cm.true_x1, Cm.true_x2)): ax[2].plot([i - 0.4, i + 0.4], [a, a], color="#2E294E", lw=1.2, ls="--"); ax[2].plot([i - 0.4, i + 0.4], [b, b], color="#2E294E", lw=1.2, ls=":")
ax[2].set_xticks(x); ax[2].set_xticklabels(["C1", "C2"]); ax[2].set(ylabel="Mean estimate (dashed: true X1, dotted: true X2)", title="Multivariable MR separates exposures"); ax[2].legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig("figs11/fig9_mode_mvmr.png", dpi=200); print("figure written")
