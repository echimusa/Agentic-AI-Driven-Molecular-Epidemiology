"""
Chapter 11, Hands-on 7 - reconciling an observational association with a null Mendelian randomisation result (the "think critically" question of Session 11), with a known truth.
The slide offers three explanations: confounding in the observational data, reverse causation, and weak instruments. Five scenarios are simulated, 200 replicates each, with two independent samples of 50,000 people
(one gives the exposure effects, the other the outcome effects and the observational association), 30 independent variants, and a standardised genotype:
  S1  confounding only: a common cause U raises X and Y; X has no effect on Y
  S2  reverse causation: Y raises X (0.4 per unit); X has no effect on Y; the variants act on X directly
  S3  reverse causation with the wrong instruments: the variants act on Y and Y raises X (0.5 per unit), so they are associated with X only through Y; X has no effect on Y
  S4  a true effect (0.3) with weak instruments (mean F about 3)
  S5  a true effect (0.3) with strong instruments (mean F about 44)
For each: the observational slope, the inverse-variance weighted MR estimate, the share of replicates whose 95% interval excludes zero, and the Steiger test of direction (the variance each variant explains in X and in Y).
Output: Figure 6, table_reconcile.csv.   Needs numpy, pandas, scipy, matplotlib and mr11.py only (no downloads).
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
N, M, REPS = 50_000, 30, 200; rng = np.random.default_rng(11)
def sample(n, maf, mode, b, cx, cy, r, gam, delta):
    G = rng.binomial(2, maf, (n, M)).astype(float); G = (G - 2 * maf) / np.sqrt(2 * maf * (1 - maf)); U = rng.normal(size=n)
    if mode == "forward": X = G @ gam + cx * U + rng.normal(size=n); Y = b * X + cy * U + G @ delta + rng.normal(size=n)
    else: Y = cy * U + G @ delta + rng.normal(size=n); X = G @ gam + cx * U + r * Y + rng.normal(size=n)
    return G, X, Y
def gwas(G, v): n = len(v); b = G.T @ (v - v.mean()) / n; se = np.sqrt((v.var() - b ** 2) / n); return b, se
SC = [("S1 Confounding only", dict(mode="forward", b=0.0, cx=0.6, cy=0.6, r=0, g=np.sqrt(60 / N), d=0.0), 0.0), ("S2 Reverse causation (Y raises X)", dict(mode="reverse", b=0, cx=0.0, cy=0.6, r=0.4, g=np.sqrt(60 / N), d=0.0), 0.0),
      ("S3 Reverse causation, instruments act on Y", dict(mode="reverse", b=0, cx=0.0, cy=0.6, r=0.5, g=0.0, d=np.sqrt(60 / N) / 0.5), 0.0), ("S4 True effect, weak instruments", dict(mode="forward", b=0.3, cx=0.6, cy=0.6, r=0, g=np.sqrt(3 / N), d=0.0), 0.3),
      ("S5 True effect, strong instruments", dict(mode="forward", b=0.3, cx=0.6, cy=0.6, r=0, g=np.sqrt(60 / N), d=0.0), 0.3)]
rows = []
for name, p, truth in SC:
    obs, obs_se, ivw, ivw_se, share, F, sg = [], [], [], [], [], [], []
    for _ in range(REPS):
        maf = rng.uniform(0.1, 0.5, M); gam = np.full(M, p["g"]) * rng.choice([-1, 1], M); delta = np.full(M, p["d"]) * rng.choice([-1, 1], M) * (1 if p["d"] else 0)
        G1, X1, Y1 = sample(N, maf, p["mode"], p["b"], p["cx"], p["cy"], p["r"], gam, delta); G2, X2, Y2 = sample(N, maf, p["mode"], p["b"], p["cx"], p["cy"], p["r"], gam, delta)
        bx, sx = gwas(G1, X1); by, sy = gwas(G2, Y2); s = np.sign(bx); bx, by = bx * s, by * s            # orient every variant to the exposure-raising allele
        I = mr.ivw(bx, by, sy); ivw.append(I["beta"]); ivw_se.append(I["se_random"]); F.append(np.mean(mr.fstat(bx, sx)))
        ox = np.cov(X2, Y2)[0, 1] / X2.var(); obs.append(ox); obs_se.append(np.sqrt((Y2.var() - ox ** 2 * X2.var()) / (N * X2.var())))
        r2x = mr.r2_from_z(bx / sx, N); r2y = mr.r2_from_z(by / sy, N); share.append(np.mean(r2x > r2y)); sg.append(np.sum(r2x) > np.sum(r2y))
    obs, obs_se, ivw, ivw_se = map(np.array, (obs, obs_se, ivw, ivw_se))
    rows.append(dict(scenario=name, truth=truth, mean_F=np.mean(F), observational=obs.mean(), observational_excludes_0=float(np.mean(np.abs(obs / obs_se) > 1.96)), mr_ivw=ivw.mean(), mr_ivw_mean_se=ivw_se.mean(), mr_excludes_0=float(np.mean(np.abs(ivw / ivw_se) > 1.96)),
                     mr_covers_truth=float(np.mean(np.abs(ivw - truth) < 1.96 * ivw_se)), steiger_share_variants=float(np.mean(share)), steiger_direction_x_to_y=float(np.mean(sg)), reps=REPS))
T = pd.DataFrame(rows); T.to_csv("table_reconcile.csv", index=False); print(T.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3)); y = np.arange(len(T))[::-1]; w = 0.36
ax[0].barh(y + w / 2, T.observational, w, color="#a0aec0", label="Observational slope"); ax[0].barh(y - w / 2, T.mr_ivw, w, color="#E4572E", label="MR (IVW)"); ax[0].errorbar(T.mr_ivw, y - w / 2, xerr=1.96 * T.mr_ivw_mean_se, fmt="none", ecolor="#2E294E", elinewidth=1.2, capsize=3)
for t_, yy in zip(T.truth, y): ax[0].plot([t_, t_], [yy - 0.5, yy + 0.5], color="#1B998B", lw=2)
ax[0].set_yticks(y); ax[0].set_yticklabels([n.split(" ", 1)[0] for n in T.scenario]); ax[0].set(xlabel="Estimate of the effect of X on Y (green line: truth)", title="Same observational signal, different truths"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper right")
ax[1].bar(range(len(T)), T.steiger_direction_x_to_y * 100, color=["#1B998B", "#1B998B", "#E4572E", "#1B998B", "#1B998B"]); ax[1].set_xticks(range(len(T))); ax[1].set_xticklabels([n.split(" ", 1)[0] for n in T.scenario]); ax[1].set(ylabel="Replicates where Steiger says X → Y (%)", title="The direction check catches S3", ylim=(0, 105))
fig.tight_layout(); fig.savefig("figs11/fig6_reconcile.png", dpi=200); print("figure written")
