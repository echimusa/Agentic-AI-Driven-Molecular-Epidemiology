"""
Chapter 11, Hands-on 11 - a within-family design against dynastic effects (Yao et al. 2025, Section 7; Brumpton et al. 2020).
Parents pass one of their two alleles at each of 20 variants to each child, independently. A genetic score S sums the alleles. Two siblings per family, 20,000 families, 100 replicates.
  X = a S_child + c S_parents + noise;   Y = b X + d S_parents + U-free noise.   S_parents is the mean score of the two parents: it affects the child through the home environment (a dynastic effect), not through the child's own genes.
  F1  no dynastic effect:  b = 0.3, c = 0,   d = 0
  F2  dynastic effects only: b = 0,   c = 0.3, d = 0.3
  F3  a true effect and dynastic effects: b = 0.3, c = 0.3, d = 0.3
Two instrumental-variable estimates: population-based (first sibling only; ratio of cov(S, Y) to cov(S, X)) and within-family (sibling differences; ratio of cov(dS, dY) to cov(dS, dX)).
Output: Figure 10 and table_family.csv.   Needs numpy, pandas, matplotlib only (no downloads).
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
rng = np.random.default_rng(15); NF, M, REPS, A = 20_000, 20, 100, 0.2
def family(b, c, d):
    maf = rng.uniform(0.2, 0.5, M); hap = rng.random((NF, 2, 2, M)) < maf                      # families x parent x haplotype x variant
    def child(): pick = rng.integers(0, 2, (NF, 2)); return hap[np.arange(NF)[:, None], np.arange(2)[None, :], pick].sum(axis=1)       # one allele from each parent
    c1, c2 = child(), child(); par = hap.sum(axis=2).sum(axis=1)
    def sd(v): return (v - v.mean()) / v.std()
    S1, S2 = sd(c1.sum(axis=1).astype(float)), sd(c2.sum(axis=1).astype(float)); SP = sd((par.sum(axis=1)).astype(float)); out = []
    for S in (S1, S2): X = A * S + c * SP + rng.normal(size=NF); Y = b * X + d * SP + rng.normal(size=NF); out.append((S, X, Y))
    return out
rows = []
for name, (b, c, d) in (("F1 No dynastic effect", (0.3, 0.0, 0.0)), ("F2 Dynastic effects only (true effect 0)", (0.0, 0.3, 0.3)), ("F3 True effect and dynastic effects", (0.3, 0.3, 0.3))):
    pop, win = [], []
    for _ in range(REPS):
        (S1, X1, Y1), (S2, X2, Y2) = family(b, c, d); pop.append(np.cov(S1, Y1)[0, 1] / np.cov(S1, X1)[0, 1]); dS, dX, dY = S1 - S2, X1 - X2, Y1 - Y2; win.append(np.cov(dS, dY)[0, 1] / np.cov(dS, dX)[0, 1])
    rows.append(dict(scenario=name, true_effect=b, population=np.mean(pop), population_sd=np.std(pop), within_family=np.mean(win), within_family_sd=np.std(win), reps=REPS))
T = pd.DataFrame(rows); T.to_csv("table_family.csv", index=False); print(T.round(3).to_string(index=False))
fig, ax = plt.subplots(figsize=(7.6, 4.2)); x = np.arange(len(T)); ax.bar(x - 0.2, T.population, 0.4, yerr=T.population_sd, color="#E4572E", label="Population-based", error_kw=dict(ecolor="#2E294E", capsize=3)); ax.bar(x + 0.2, T.within_family, 0.4, yerr=T.within_family_sd, color="#1B998B", label="Within-family", error_kw=dict(ecolor="#2E294E", capsize=3))
for i, t in enumerate(T.true_effect): ax.plot([i - 0.45, i + 0.45], [t, t], color="#2E294E", lw=1.4, ls="--")
ax.set_xticks(x); ax.set_xticklabels(["F1", "F2", "F3"]); ax.set(ylabel="Estimate of the effect of X on Y (dashed: truth)", title="Dynastic effects bias the population estimate, not the within-family one"); ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(); fig.savefig("figs11/fig10_family.png", dpi=200); print("figure written")
