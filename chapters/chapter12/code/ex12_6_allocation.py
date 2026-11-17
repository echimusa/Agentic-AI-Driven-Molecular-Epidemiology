"""
Chapter 12, Hands-on 6 - Where the sequences are taken decides when a new variant is seen.
A new variant starts at 0.1% of infections in one World Bank income group and doubles every week up to 50%. A fixed global capacity (5,000 sequences per week) is split across the four income groups in three ways:
  current   - in proportion to each group's share of public Salmonella isolates (Hands-on 2), as a stand-in for today's sequencing effort;
  population - in proportion to each group's share of the world population;
  floor     - at least 10% of capacity to every group, the remaining 60% in proportion to the current pattern.
In week t the chance that at least one of the n sequences from the group carries the variant is 1 - (1 - f_t)^n, so the distribution of the detection week is exact.
Input: table_income_groups.csv (Hands-on 2). Output: table_allocation.csv, figs12/fig6_allocation.png.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
g = pd.read_csv("table_income_groups.csv"); CAP = 5000; F0 = 0.001; DOUBLE_WEEKS = 1; CAPF = 0.5; WEEKS = 80
shares = {"current": (g.isolates / g.isolates.sum()).values, "population": (g.population / g.population.sum()).values}; shares["floor"] = 0.10 + 0.60 * shares["current"]
for k, v in shares.items(): assert abs(v.sum() - 1) < 1e-9, k
def detection(n, f0=F0):
    t = np.arange(WEEKS); f = np.minimum(CAPF, f0 * 2 ** (t / DOUBLE_WEEKS)); hit = 1 - (1 - f) ** n; S = np.cumprod(1 - hit); cdf = 1 - S
    med = int(np.argmax(cdf >= 0.5)) if (cdf >= 0.5).any() else None; p90 = int(np.argmax(cdf >= 0.9)) if (cdf >= 0.9).any() else None
    return med, p90, (f[med] if med is not None else np.nan), cdf
rows = []; curves = {}
for sc, sh in shares.items():
    for i, grp in enumerate(g.income):
        n = CAP * sh[i]; med, p90, f_at, cdf = detection(n); curves[(sc, grp)] = cdf
        rows.append(dict(allocation=sc, emerges_in=grp, sequences_per_week=round(n, 1), median_week_of_detection=med, p90_week=p90, local_frequency_at_median_pct=round(100 * f_at, 1)))
T = pd.DataFrame(rows); T.to_csv("table_allocation.csv", index=False); print(T.to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]})
cols = {"current": "#c53030", "population": "#2b6cb0", "floor": "#d69e2e"}; lab = [s.replace(" income", "").replace(" middle", "\nmiddle") for s in g.income]
x = np.arange(4); w = 0.26
for j, sc in enumerate(shares): ax[0].bar(x + (j - 1) * w, CAP * shares[sc], w, color=cols[sc], label=sc)
ax[0].set(yscale="log", xticks=x, xticklabels=lab, ylabel="Sequences per week (log scale)", title="A. Capacity by income group"); ax[0].legend(frameon=False, fontsize=8.8)
for j, sc in enumerate(shares): ax[1].bar(x + (j - 1) * w, [T[(T.allocation == sc) & (T.emerges_in == gr)].median_week_of_detection.iloc[0] for gr in g.income], w, color=cols[sc])
ax[1].set(xticks=x, xticklabels=lab, xlabel="Group where the variant emerges", ylabel="Median week of detection", title="B. Weeks until the variant is seen")
for gr, ls in ((g.income.iloc[3], "-"),): 
    for sc in shares: ax[2].plot(np.arange(WEEKS), curves[(sc, gr)] * 100, color=cols[sc], lw=2, ls=ls, label=sc)
ax[2].set(xlim=(0, 40), ylim=(0, 100), xlabel="Weeks since the variant emerged", ylabel="Cumulative probability of detection (%)", title="C. Emergence in the low-income group"); ax[2].axhline(90, color="grey", lw=0.8, ls=":")
fig.tight_layout(); fig.savefig("figs12/fig6_allocation.png", dpi=200); plt.close(fig); print("figure written")
