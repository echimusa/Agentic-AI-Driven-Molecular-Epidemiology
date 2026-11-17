"""
Chapter 13, Hands-on 7 - A partly sampled outbreak with an epidemiological attribute (the setting of Carson et al. 2025, as a toy).
Cases are of type A (mean 2.4 secondary cases) or type B (mean 1.2): a true ratio of 2. A new infection is of type A with probability 0.3, or, with assortative mixing, takes its infector's type with probability 0.8. Each case is sampled with a probability that may depend on its type. Two analyses of the same sampled cases:
  direct - count only transmission links whose two ends are both sampled, and compare the sampled offspring per sampled case of type A and type B (the estimate needs no tree reconstruction but uses only the observed links);
  collapsed - attach each sampled case to its nearest sampled ancestor, as an analysis that ignores the unsampled intermediate hosts would.
Each outbreak is simulated by generations until it has at least 1,500 cases; the last generation is left out of the estimates because its offspring are not yet generated. 200 outbreaks per scenario.
Output: table_partial.csv, figs13/fig7_partial.png.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
RB, RA, PA, TARGET, REPS = 1.2, 2.4, 0.3, 1500, 200
def outbreak(rng, assort=0.0):
    while True:
        inf = [-1]; typ = [0]; gen = [0]; cur = np.array([0]); g = 0
        while len(inf) < TARGET and len(cur):
            cnt = rng.poisson(np.where(np.array(typ)[cur] == 1, RA, RB)); tot = int(cnt.sum()); par = np.repeat(cur, cnt); ty = np.where(rng.random(tot) < assort, np.array(typ)[par], (rng.random(tot) < PA).astype(int)); g += 1
            new = np.arange(len(inf), len(inf) + tot); inf += par.tolist(); typ += ty.tolist(); gen += [g] * tot; cur = new
        if len(inf) >= 300: return np.array(inf), np.array(typ), np.array(gen)
def analyse(rng, piA, piB, assort=0.0):
    inf, typ, gen = outbreak(rng, assort); n = len(inf); smp = rng.random(n) < np.where(typ == 1, piA, piB); last = gen.max()
    nsa = np.full(n, -1); unsamp = np.zeros(n, int)
    for j in range(1, n):
        p = inf[j]; nsa[j], unsamp[j] = (p, 0) if smp[p] else (nsa[p], unsamp[p] + 1)
    has_inf = gen < last                                                             # infectors whose offspring are fully generated
    links = np.arange(1, n)[has_inf[inf[1:]]]; both = links[smp[links] & smp[inf[links]]]
    o = np.bincount(inf[both], minlength=n); c = np.bincount(nsa[(nsa >= 0) & smp], minlength=n)    # sampled children, direct and collapsed
    ok = smp & has_inf; res = {}
    for name, v in (("direct", o), ("collapsed", c)):
        mA = v[ok & (typ == 1)].mean() if (ok & (typ == 1)).any() else np.nan; mB = v[ok & (typ == 0)].mean() if (ok & (typ == 0)).any() else np.nan; res[name] = (mA, mB)
    pbar = PA * piA + (1 - PA) * piB; with_anc = smp & (nsa >= 0)
    return dict(share_links=len(both) / len(links), direct_links=len(both), ratio_direct=res["direct"][0] / res["direct"][1] if res["direct"][1] > 0 else np.nan, ratio_collapsed=res["collapsed"][0] / res["collapsed"][1] if res["collapsed"][1] > 0 else np.nan,
                RA_direct=res["direct"][0] / pbar, RB_direct=res["direct"][1] / pbar, RA_collapsed=res["collapsed"][0], RB_collapsed=res["collapsed"][1], share_unsampled=(unsamp[with_anc] > 0).mean(), mean_unsampled=unsamp[with_anc].mean(), n_cases=n, n_sampled=int(smp.sum()))
scen = [("pi = 1.0 (all cases sampled)", 1.0, 1.0, 0.0), ("pi = 0.5", 0.5, 0.5, 0.0), ("pi = 0.3", 0.3, 0.3, 0.0), ("pi = 0.1", 0.1, 0.1, 0.0), ("pi_A = 0.6, pi_B = 0.3", 0.6, 0.3, 0.0), ("Assortative mixing, pi = 0.3", 0.3, 0.3, 0.8), ("Assortative mixing, pi_A = 0.6, pi_B = 0.3", 0.6, 0.3, 0.8)]
rows = []; store = {}
for k, (name, a, b, asr) in enumerate(scen):
    rng = np.random.default_rng(1300 + k); D = pd.DataFrame([analyse(rng, a, b, asr) for _ in range(REPS)]); store[name] = D; q = lambda s: np.nanpercentile(s, [2.5, 50, 97.5])
    rd, rc = q(D.ratio_direct), q(D.ratio_collapsed)
    rows.append(dict(scenario=name, pi_A=a, pi_B=b, assortative=asr, share_of_links_observed=D.share_links.mean(), direct_links=D.direct_links.mean(), ratio_direct_median=rd[1], ratio_direct_lo=rd[0], ratio_direct_hi=rd[2], ratio_collapsed_median=rc[1], ratio_collapsed_lo=rc[0], ratio_collapsed_hi=rc[2],
                     share_with_unsampled_intermediate=D.share_unsampled.mean(), mean_unsampled_between=D.mean_unsampled.mean(), mean_cases=D.n_cases.mean()))
T = pd.DataFrame(rows); T.to_csv("table_partial.csv", index=False); print(T.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]})
lab = ["1.0", "0.5", "0.3", "0.1", "A 0.6\nB 0.3", "assort.\n0.3", "assort.\nA 0.6\nB 0.3"]; x = np.arange(7)
ax[0].bar(x, 100 * T.share_of_links_observed, color="#2b6cb0"); ax[0].plot(x[:4], [100 * p * p for p in (1.0, 0.5, 0.3, 0.1)], "k_", ms=22, mew=2); ax[0].set(xticks=x, xticklabels=lab, xlabel="Sampling probability", ylabel="True transmission links with both ends sampled (%)", title="A. How much of the tree is seen"); ax[0].text(3.6, 80, "bars: simulated\nmarks: pi squared", fontsize=8.8)
for k, (c1, c2) in enumerate(zip(("ratio_direct", "ratio_collapsed"), ("#2b6cb0", "#d69e2e"))):
    ax[1].errorbar(x + (k - 0.5) * 0.25, T[c1 + "_median"], yerr=[T[c1 + "_median"] - T[c1 + "_lo"], T[c1 + "_hi"] - T[c1 + "_median"]], fmt="o", color=c2, ms=6, capsize=3, label=["Direct links only", "Nearest sampled ancestor"][k])
ax[1].axhline(2.0, color="black", ls="--", lw=0.8); ax[1].set(xticks=x, xticklabels=lab, xlabel="Sampling probability", ylabel="Estimated ratio of reproduction numbers, A to B", title="B. The ratio (true value 2)", ylim=(0, 6.5)); ax[1].legend(frameon=False, fontsize=8.5, loc="upper left")
ax[2].bar(x, T.mean_unsampled_between, color="#c53030"); ax[2].set(xticks=x, xticklabels=lab, xlabel="Sampling probability", ylabel="Unsampled hosts between a sampled case\nand its nearest sampled ancestor", title="C. Hidden intermediates")
fig.tight_layout(); fig.savefig("figs13/fig7_partial.png", dpi=200); plt.close(fig); print("figure written")
