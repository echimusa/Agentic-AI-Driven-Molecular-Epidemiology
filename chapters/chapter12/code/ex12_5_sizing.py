"""
Chapter 12, Hands-on 5 - How many sequences per week does a surveillance aim need?
A. The detection formula of the Session 12 demo, n = ceil(log(1 - c) / log(1 - p)), its rule-of-three shortcut and a Monte Carlo check.
B. The smallest variant frequency that a weekly number of sequences can detect.
C. The same arithmetic applied to the weekly volume of the public Salmonella table (human isolates collected in 2022).
D. Four departures from the formula: an imperfect assay, a finite population, sampling by household, and the catchment of a wastewater sample.
Input: sal_surv.csv.gz (Hands-on 1). Output: table_sizing.csv, table_sizing_real.csv, table_sizing_departures.csv, table_catchment.csv, figs12/fig5_sizing.png.
"""
import math
import numpy as np, pandas as pd
from scipy import stats
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
detect_n = lambda p, c=0.95: math.ceil(math.log(1 - c) / math.log(1 - p))          # sequences needed to see a variant of frequency p with confidence c
p_min = lambda n, c=0.95: 1 - (1 - c) ** (1 / n)                                      # smallest frequency detectable with n sequences
# ---- A. the formula, the rule of three, a Monte Carlo check
rows = []
for p in (0.005, 0.01, 0.02, 0.05, 0.10):
    rows.append(dict(frequency_pct=p * 100, n_90=detect_n(p, 0.90), n_95=detect_n(p, 0.95), n_99=detect_n(p, 0.99), rule_of_three_95=round(3 / p), weeks_at_96_per_week_95=round(detect_n(p, 0.95) / 96, 2)))
A = pd.DataFrame(rows); A.to_csv("table_sizing.csv", index=False); print(A.to_string(index=False))
rng = np.random.default_rng(12); mc = {n: (rng.binomial(n, 0.02, 400_000) >= 1).mean() for n in (148, 149)}
print(f"n = 149 for 2% at 95%: exact {1 - 0.98 ** 149:.4f}, simulated {mc[149]:.4f}; n = 148: exact {1 - 0.98 ** 148:.4f}, simulated {mc[148]:.4f}")
print(f"capacity 96 per week detects frequencies down to {p_min(96) * 100:.2f}% at 95%")
# ---- C. real weekly volumes
m = load_surv(); h = m[(m.sector == "Human") & (m.year == 2022) & m.country.notna()]
cnt = h.country.value_counts().head(8); real = []
for ctry, n in cnt.items():
    w = n / 52.14; real.append(dict(country=ctry, isolates_2022=int(n), per_week=round(w, 1), smallest_frequency_pct_in_a_week=round(p_min(w) * 100, 2), weeks_to_detect_2pct=round(detect_n(0.02) / w, 2)))
R = pd.DataFrame(real); R.to_csv("table_sizing_real.csv", index=False); print(R.to_string(index=False))
# ---- D. departures from the formula, all for a 2% variant and 95% confidence
p, c = 0.02, 0.95; base = detect_n(p, c); dep = [dict(setting="Formula: random sample, perfect assay, infinite population", sequences_needed=base, note="Session 12 demo")]
s = 0.8; dep.append(dict(setting=f"Assay detects a variant carrier with probability {s}", sequences_needed=detect_n(p * s, c), note="replace p by p x sensitivity"))
N, K = 2000, 40
def hyper_n(N, K, c):
    for n in range(1, N + 1):
        if 1 - stats.hypergeom.cdf(0, N, K, n) >= c: return n
dep.append(dict(setting=f"Finite population of {N} infections, {K} carry the variant", sequences_needed=hyper_n(N, K, c), note="sampling without replacement"))
m_h, rho = 10, 0.3; a_ = p * (1 - rho) / rho; b_ = (1 - p) * (1 - rho) / rho                  # beta-binomial: household frequency varies with intraclass correlation rho
def det_prob_households(H): return 1 - stats.betabinom.pmf(0, m_h, a_, b_) ** H
Hneed = next(H for H in range(1, 2000) if det_prob_households(H) >= c)
dep.append(dict(setting=f"Households of {m_h}, intraclass correlation {rho}: {Hneed} households", sequences_needed=Hneed * m_h, note=f"design effect 1 + (m - 1) x rho = {1 + (m_h - 1) * rho:.1f} gives {round(base * (1 + (m_h - 1) * rho))}"))
dep.append(dict(setting=f"Same {base} sequences taken as {math.ceil(base / m_h)} households of {m_h}: detection probability", sequences_needed=base, note=f"{det_prob_households(math.ceil(base / m_h)) * 100:.1f}% instead of {(1 - (1 - p) ** base) * 100:.1f}%"))
D = pd.DataFrame(dep); D.to_csv("table_sizing_departures.csv", index=False); print(D.to_string(index=False))
cat = []
for q in (1e-4, 1e-5):
    for Npop in (10_000, 1_000_000): cat.append(dict(shedder_prevalence=q, catchment=Npop, expected_shedders=q * Npop, p_at_least_one=1 - (1 - q) ** Npop))
C_ = pd.DataFrame(cat); C_.to_csv("table_catchment.csv", index=False); print(C_.to_string(index=False))
# ---- figure
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.1, 1.0, 1.0]})
ps = np.logspace(-3, -0.7, 200)
for cc, col in ((0.90, "#a0aec0"), (0.95, "#2b6cb0"), (0.99, "#c53030")): ax[0].plot(ps * 100, np.log(1 - cc) / np.log(1 - ps), color=col, lw=2, label=f"{cc:.0%} confidence")
ax[0].scatter([2], [base], color="black", zorder=5); ax[0].annotate(f"2%: {base} sequences", (2, base), xytext=(3.5, 700), arrowprops=dict(arrowstyle="-", color="black"), fontsize=9)
ax[0].set(xscale="log", yscale="log", xlabel="Variant frequency in the sampled population (%)", ylabel="Sequences needed", title="A. Sequences needed to see a variant"); ax[0].legend(frameon=False, fontsize=8.8)
ax[1].barh(R.country[::-1], R.smallest_frequency_pct_in_a_week[::-1], color="#2b6cb0"); ax[1].set(xscale="log", xlim=(0.1, 200), xlabel="Smallest frequency detectable in a week\n(%, 95% confidence, log scale)", title="C. Weekly volume of public isolates, 2022")
for i_, v in enumerate(R.smallest_frequency_pct_in_a_week[::-1]): ax[1].text(v * 1.1, i_, f"{v:.2f}", va="center", fontsize=8.8)
short = ["Formula", "Assay sensitivity 0.8", "Finite population (2,000)", "35 households of 10\n(correlation 0.3)"]
ax[2].barh(short[::-1], D.sequences_needed[:4][::-1], color=["#2b6cb0", "#d69e2e", "#d69e2e", "#c53030"][::-1]); ax[2].set(xlabel="Sequences needed\n(2% variant, 95% confidence)", title="D. Departures from the formula")
for i_, v in enumerate(D.sequences_needed[:4][::-1]): ax[2].text(v + 6, i_, str(int(v)), va="center", fontsize=8.8)
ax[2].set_xlim(0, 420); ax[2].tick_params(axis="y", labelsize=8.8)
fig.tight_layout(); fig.savefig("figs12/fig5_sizing.png", dpi=200); plt.close(fig); print("figure written")
