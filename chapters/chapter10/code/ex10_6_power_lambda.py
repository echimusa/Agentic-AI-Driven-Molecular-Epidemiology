"""
Chapter 10, Hands-on 6 - how many people? Power at the genome-wide threshold, and how the genomic-inflation factor grows with the sample size.
Part A  Power of the allelic case-control test (1 degree of freedom) from the non-central chi-square distribution: for N_ca cases and N_co controls, control allele frequency p0 and odds ratio OR,
        p1 = OR p0 / (1 - p0 + OR p0), pbar = (N_ca p1 + N_co p0) / (N_ca + N_co), non-centrality = (p1 - p0)^2 / (pbar (1 - pbar) (1 / (2 N_ca) + 1 / (2 N_co))).
        Tables: cases needed (equal numbers of controls) for 80% power at P < 5e-8; the smallest odds ratio detectable in a 1,000 / 1,000 study and in the COVID-19 HGI file (9,986 cases, 1,877,672 controls);
        a Monte Carlo check of the formula (binomial allele counts, 4,000 repeats per scenario).
Part B  A null scan with ancestry imbalance, in two scenarios: cases 30% EAS and controls 10% EAS (the rest EUR), and cases 11% and controls 9% EAS; N/2 of each, and the allele frequencies are those of the 9,314 common chromosome 22 variants
        of Hands-on 5. The test statistics are simulated for N = 1,000, 10,000 and 100,000 people and the inflation factor is compared with 1 + N F, where F is the mean over variants of
        (p_case - p_control)^2 / (2 pbar (1 - pbar)) (the course notes give lambda = 1 + N FST for N/2 cases and N/2 controls).  Then the rescaling of the HGI lambda of Hands-on 2 is checked.
Input: chr22_group_geno_counts.tsv.gz and common.txt (Hands-on 5).  Output: table_power_cases_needed.csv, table_power_min_or.csv, table_power_mc.csv, table_lambda_scaling.csv, table_example_allelic.csv, Figure 7.
"""
import numpy as np, pandas as pd
from scipy import stats
ALPHA = 5e-8; CRIT = stats.chi2.isf(ALPHA, 1)
def ncp(nca, nco, p0, orr):
    p1 = orr * p0 / (1 - p0 + orr * p0); pb = (nca * p1 + nco * p0) / (nca + nco); return (p1 - p0) ** 2 / (pb * (1 - pb) * (1 / (2 * nca) + 1 / (2 * nco)))
def power(nca, nco, p0, orr, crit=CRIT): return stats.ncx2.sf(crit, 1, ncp(nca, nco, p0, orr))
def cases_needed(p0, orr, ratio=1.0, target=0.8):
    lo, hi = 10, 1e9
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if power(mid, mid * ratio, p0, orr) < target: lo = mid
        else: hi = mid
    return hi
def min_or(nca, nco, p0, target=0.8):
    lo, hi = 1.0001, 20.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if power(nca, nco, p0, mid) < target: lo = mid
        else: hi = mid
    return hi
# ---------------------------------------------------------------- the example of the lecture notes: 400 of 1,000 case alleles and 300 of 1,000 control alleles are A
a_c, t_c, a_o, t_o = 400, 600, 300, 700; n_ = a_c + t_c + a_o + t_o; chi_ex = n_ * (a_c * t_o - t_c * a_o) ** 2 / ((a_c + t_c) * (a_o + t_o) * (a_c + a_o) * (t_c + t_o))
pd.DataFrame([dict(case_alleles_a=a_c, case_alleles_t=t_c, control_alleles_a=a_o, control_alleles_t=t_o, chi_square=chi_ex, p_value=stats.chi2.sf(chi_ex, 1), odds_ratio=(a_c * t_o) / (t_c * a_o), cases=500, controls=500, power_at_5e8=power(500, 500, 0.3, (a_c * t_o) / (t_c * a_o)))]).to_csv("table_example_allelic.csv", index=False)
print(f"example of the notes: chi-square {chi_ex:.2f}, P {stats.chi2.sf(chi_ex, 1):.1e}, odds ratio {(a_c * t_o) / (t_c * a_o):.2f}, power at 5e-8 {power(500, 500, 0.3, (a_c * t_o) / (t_c * a_o)):.2f}")
MAFS = [0.01, 0.05, 0.20, 0.50]; ORS = [1.1, 1.2, 1.5, 2.0]
A = pd.DataFrame([dict(odds_ratio=o, **{f"maf_{m}": int(round(cases_needed(m, o))) for m in MAFS}) for o in ORS]); A.to_csv("table_power_cases_needed.csv", index=False); print("cases needed (equal controls), 80% power, P < 5e-8"); print(A.to_string(index=False))
B = pd.DataFrame([dict(maf=m, min_or_1000_1000=min_or(1000, 1000, m), min_or_hgi_r5=min_or(9986, 1877672, m), min_or_10000_10000=min_or(10000, 10000, m)) for m in MAFS]); B.to_csv("table_power_min_or.csv", index=False); print("smallest odds ratio with 80% power at P < 5e-8"); print(B.round(3).to_string(index=False))
rng = np.random.default_rng(10); rows = []
for nca, nco, p0, orr in [(1000, 1000, 0.3, 2.0), (5000, 5000, 0.2, 1.3), (10000, 10000, 0.1, 1.2), (20000, 20000, 0.3, 1.1), (3000, 30000, 0.15, 1.35)]:
    p1 = orr * p0 / (1 - p0 + orr * p0); reps = 4000; ca = rng.binomial(2 * nca, p1, reps); co = rng.binomial(2 * nco, p0, reps)
    a1, a2 = ca / (2 * nca), co / (2 * nco); pb = (ca + co) / (2 * nca + 2 * nco); chi = (a1 - a2) ** 2 / (pb * (1 - pb) * (1 / (2 * nca) + 1 / (2 * nco)))
    rows.append(dict(cases=nca, controls=nco, maf=p0, odds_ratio=orr, power_formula=power(nca, nco, p0, orr), power_monte_carlo=float((chi > CRIT).mean()), repeats=reps))
M = pd.DataFrame(rows); M.to_csv("table_power_mc.csv", index=False); print(M.round(3).to_string(index=False))
# ---------------------------------------------------------------- Part B
c = pd.read_csv("chr22_group_geno_counts.tsv.gz", sep="\t"); com = set(open("common.txt").read().split()); c = c[c.SNP.isin(com)]
def pf(g): n = c[f"{g}_n11"] + c[f"{g}_n12"] + c[f"{g}_n22"]; return ((2 * c[f"{g}_n11"] + c[f"{g}_n12"]) / (2 * n)).values
pe, pa = pf("EUR"), pf("EAS"); rows = []
for qca, qco in [(0.30, 0.10), (0.11, 0.09)]:
    pca = qca * pa + (1 - qca) * pe; pco = qco * pa + (1 - qco) * pe; pb = (pca + pco) / 2; keep = (pb > 0) & (pb < 1)
    F = float(np.mean((pca[keep] - pco[keep]) ** 2 / (2 * pb[keep] * (1 - pb[keep])))); print(f"cases {qca:.0%} EAS, controls {qco:.0%} EAS: {keep.sum():,} variants; F = {F:.6f}")
    for N in (1000, 10000, 100000):
        chis = []
        for rep in range(5):
            ca = rng.binomial(N, pca[keep]); co = rng.binomial(N, pco[keep]); q = (ca + co) / (2 * N); a1, a2 = ca / N, co / N; ok = (q > 0) & (q < 1)
            chis.append((a1[ok] - a2[ok]) ** 2 / (q[ok] * (1 - q[ok]) * (2 / N)))
        chi = np.concatenate(chis); rows.append(dict(cases_eas_share=qca, controls_eas_share=qco, people=N, F=F, lambda_median=float(np.median(chi) / stats.chi2.ppf(0.5, 1)), lambda_mean=float(chi.mean()), one_plus_N_F=1 + N * F, lambda_hgi_rescaled=np.nan))
L = pd.DataFrame(rows); n_eff = 4 / (1 / 9986 + 1 / 1877672); f_hgi = (1.045 - 1) / n_eff
L = pd.concat([L, pd.DataFrame([dict(cases_eas_share=np.nan, controls_eas_share=np.nan, people=int(round(n_eff)), F=f_hgi, lambda_median=1.045, lambda_mean=np.nan, one_plus_N_F=1 + f_hgi * n_eff, lambda_hgi_rescaled=1 + f_hgi * 2000)])], ignore_index=True)
L.to_csv("table_lambda_scaling.csv", index=False); print(L.round(4).to_string(index=False)); print(f"HGI release 5: effective N {n_eff:,.0f}, implied F {f_hgi:.2e}, lambda at 1,000 + 1,000 people {1 + f_hgi * 2000:.4f}")

# ---------------------------------------------------------------- Figure 7
import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.3)); cols = ["#2f855a", "#2b6cb0", "#dd6b20", "#c53030"]
for m, c_ in zip(MAFS, cols): ax[0].plot(A.odds_ratio, A[f"maf_{m}"], "o-", color=c_, label=f"risk allele frequency {m}")
ax[0].set_yscale("log"); ax[0].set(xlabel="Odds ratio", ylabel="Cases needed (and as many controls)", title="Cases for 80% power at P < 5e-8"); ax[0].legend(frameon=False, fontsize=8.5)
for (qa, qc), c_ in zip([(0.30, 0.10), (0.11, 0.09)], ["#E4572E", "#1B998B"]):
    d = L[(L.cases_eas_share == qa) & (L.controls_eas_share == qc)]; ax[1].plot(d.people, d.lambda_median, "o", color=c_, label=f"simulated, cases {qa:.0%} and controls {qc:.0%} EAS"); xs = np.logspace(3, 5, 20); ax[1].plot(xs, 1 + xs * d.F.iloc[0], "--", color=c_, lw=1)
ax[1].axhline(1.05, color="#718096", lw=0.8, ls=":"); ax[1].text(30000, 1.02, "λ = 1.05", color="#718096", fontsize=8.5); ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set(xlabel="People in the study", ylabel="Genomic-inflation factor λ", title="The same ancestry imbalance inflates more as N grows"); ax[1].legend(frameon=False, fontsize=8.2, loc="upper left")
fig.tight_layout(); fig.savefig("figs10/fig7_power_lambda.png", dpi=200); print("figure written")
