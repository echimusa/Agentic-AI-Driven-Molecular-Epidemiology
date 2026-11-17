"""
Chapter 10, Hands-on 8 - admixture mapping with a known truth.
In an admixed population (here African Americans, 21% +/- 14% European ancestry) a disease variant whose frequency differs between the ancestral populations makes cases carry an unusual share of one
ancestry at that locus. The lecture notes call this the ancestry odds ratio, Omega = (R pE + 1 - pE) / (R pA + 1 - pA), for a risk allele with odds ratio R and frequencies pE and pA in the European and African
source populations.
Part A  checks the formula on a simulated population of 4 million people for the six rows of the notes (R = 1.5): local ancestry of each chromosome is Bernoulli(theta) with theta the person's global ancestry,
        the risk allele follows the ancestry-specific frequency, and a person is a case with probability proportional to R per risk allele (about 2% of people are cases).
Part B  compares the power of (i) a case-control admixture-mapping test of local ancestry (Armitage trend test, adjusted for global ancestry, threshold P < 1e-5, roughly 0.05 over 5,000 ancestry blocks) with
        (ii) the test of the causal variant itself (adjusted for global ancestry, P < 5e-8), for 2,000 cases and 2,000 controls, R = 1.5, pA = 0.10 and pE = 0.10 + difference, 300 repeats per difference.
Part C  repeats B without the adjustment for global ancestry, with a disease whose risk also depends on global ancestry (an environmental effect of 0.5 log-odds per unit of European ancestry):
        the share of tests with P < 0.001 at a locus with no causal variant shows how unadjusted tests are confounded.
Output: table_admix_omega.csv, table_admix_power.csv, table_admix_confounding.csv, Figure 9.   Needs numpy, scipy, pandas, matplotlib only (no downloads).
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); import os; os.makedirs("figs10", exist_ok=True)
rng = np.random.default_rng(10)
def beta_params(mean, sd): k = mean * (1 - mean) / sd ** 2 - 1; return mean * k, (1 - mean) * k
a_, b_ = beta_params(0.21, 0.14)
# ---------------------------------------------------------------- Part A
R = 1.5; rows = []
for pE, pA in [(1.0, 0.0), (0.0, 1.0), (0.9, 0.1), (0.75, 0.25), (0.6, 0.4), (0.5, 0.5)]:
    n = 4_000_000; theta = np.full(n, 0.5); l1 = rng.random(n) < theta; l2 = rng.random(n) < theta
    k = (rng.random(n) < np.where(l1, pE, pA)).astype(np.int8) + (rng.random(n) < np.where(l2, pE, pA)).astype(np.int8)
    pc = 0.02 * R ** k.astype(float) / np.mean(0.02 * R ** k.astype(float)) * 0.02; case = rng.random(n) < np.minimum(pc, 1)
    e_case = (l1[case].sum() + l2[case].sum()); a_case = 2 * case.sum() - e_case; e_ctl = (l1[~case].sum() + l2[~case].sum()); a_ctl = 2 * (~case).sum() - e_ctl
    rows.append(dict(pE=pE, pA=pA, omega_formula=(R * pE + 1 - pE) / (R * pA + 1 - pA), omega_simulated=(e_case / a_case) / (e_ctl / a_ctl), cases=int(case.sum())))
A = pd.DataFrame(rows); A.to_csv("table_admix_omega.csv", index=False); print(A.round(3).to_string(index=False))
# ---------------------------------------------------------------- Part B and C
def partial_chi2(x, y, z):
    zc = np.column_stack([np.ones(len(z)), z]); bx = np.linalg.lstsq(zc, x, rcond=None)[0]; by = np.linalg.lstsq(zc, y, rcond=None)[0]; rx, ry = x - zc @ bx, y - zc @ by
    r = (rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)); return (len(y) - 2) * r * r
def chi2_plain(x, y): r = np.corrcoef(x, y)[0, 1]; return (len(y) - 1) * r * r
def sample(pE, pA, logR, n_pool=60_000, ncase=2000, nctl=2000, env=0.0, causal=True):
    theta = rng.beta(a_, b_, n_pool); l1 = rng.random(n_pool) < theta; l2 = rng.random(n_pool) < theta
    if causal: k = (rng.random(n_pool) < np.where(l1, pE, pA)).astype(float) + (rng.random(n_pool) < np.where(l2, pE, pA)).astype(float)
    else: k = np.zeros(n_pool)
    eta = -3.0 + logR * k + env * theta; case = rng.random(n_pool) < 1 / (1 + np.exp(-eta)); ic = np.flatnonzero(case); io = np.flatnonzero(~case)
    ic = rng.choice(ic, min(ncase, len(ic)), replace=False); io = rng.choice(io, nctl, replace=False); idx = np.concatenate([ic, io]); y = np.concatenate([np.ones(len(ic)), np.zeros(len(io))])
    return y, (l1[idx].astype(float) + l2[idx]), k[idx], theta[idx], l1, l2, theta
A_AM, A_SNP = 1e-5, 5e-8; crit_am, crit_snp = stats.chi2.isf(A_AM, 1), stats.chi2.isf(A_SNP, 1); rows = []
for delta in (0.05, 0.1, 0.2, 0.4, 0.6, 0.8):
    pA = 0.10; pE = min(pA + delta, 0.98); hit_am = hit_snp = 0; reps = 300
    for _ in range(reps):
        y, loc, k, th, *_ = sample(pE, pA, np.log(1.5)); hit_am += partial_chi2(loc, y, th) > crit_am; hit_snp += partial_chi2(k, y, th) > crit_snp
    rows.append(dict(difference=delta, pA=pA, pE=pE, power_admixture_mapping=hit_am / reps, power_causal_snp=hit_snp / reps, repeats=reps))
Pw = pd.DataFrame(rows); Pw.to_csv("table_admix_power.csv", index=False); print(Pw.round(3).to_string(index=False))
rows = []
for adj in (False, True):
    fp_loc = fp_snp = 0; reps = 400
    for _ in range(reps):
        y, loc, k, th, *_ = sample(0.5, 0.5, 0.0, env=0.5, causal=False)         # no causal variant; an ancestry-dependent disease risk
        g = (np.random.default_rng(int(rng.integers(1e9))).random(len(y)) < 0.3).astype(float) + (np.random.default_rng(int(rng.integers(1e9))).random(len(y)) < 0.3)     # an unrelated variant with no ancestry difference
        c_loc = partial_chi2(loc, y, th) if adj else chi2_plain(loc, y); c_g = partial_chi2(g, y, th) if adj else chi2_plain(g, y)
        fp_loc += c_loc > stats.chi2.isf(1e-3, 1); fp_snp += c_g > stats.chi2.isf(1e-3, 1)
    rows.append(dict(adjusted_for_global_ancestry=adj, share_local_ancestry_p_lt_0_001=fp_loc / reps, share_unrelated_snp_p_lt_0_001=fp_snp / reps, repeats=reps))
C = pd.DataFrame(rows); C.to_csv("table_admix_confounding.csv", index=False); print(C.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.1))
ax[0].plot(Pw.difference, Pw.power_admixture_mapping, "o-", color="#E4572E", label="admixture mapping (local ancestry, P < 1e-5)"); ax[0].plot(Pw.difference, Pw.power_causal_snp, "s-", color="#1B998B", label="the causal variant itself (P < 5e-8)")
ax[0].set(xlabel="Allele frequency difference between the source populations", ylabel="Power (2,000 cases, 2,000 controls, OR 1.5)", title="Admixture mapping needs a large frequency difference", ylim=(-0.02, 1.02)); ax[0].legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(0.0, 0.42))
ax[1].plot(A.omega_formula, A.omega_simulated, "o", color="#2E294E"); lim = [0.6, 1.6]; ax[1].plot(lim, lim, color="#a0aec0", lw=1); ax[1].set(xlabel="Ancestry odds ratio from the formula", ylabel="Ancestry odds ratio, simulated", title="The formula of the notes holds", xlim=lim, ylim=lim)
fig.tight_layout(); fig.savefig("figs10/fig9_admixture_mapping.png", dpi=200); print("figure written")
