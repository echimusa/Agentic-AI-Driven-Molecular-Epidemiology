"""
Chapter 1, Hands-on 3 - Design matters: how the fraction sampled and the length of the study bias a clustering estimate of recent transmission.
An endemic, tuberculosis-like process is simulated month by month: new unrelated cases arrive from outside the chain (reactivation or importation), each case causes Poisson(R) new cases after a long
gamma delay (mean 36 months), and the genotype changes at a transmission with probability m. A study samples a fraction f of the cases diagnosed in a window of D months. Two truths are computed: the share of the window's
cases whose infector was also infected in the window (direct), and the share that belong to a transmission chain with an earlier case in the window (chain). Estimate = (cases in clusters - number of clusters) / cases sampled (the n-1 method), where a cluster is two or more cases with the same genotype.
Output: table_design_simulation.csv, table_design_mutation.csv, Figure 6 of the chapter (figs1/fig3_design_simulation.png).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
R, LAM, BURN, TOTAL = 0.75, 40, 240, 240 + 60; K, TH = 2.0, 18.0
w = np.diff(stats.gamma.cdf(np.arange(0, 121), a=K, scale=TH)); w = w / w.sum()                    # delay in months from infection to the secondary case, 1..120

def simulate(rng, m):
    month, parent, geno, chain = [], [], [], []; next_g = 0; next_c = 0; pending = {}; n_cases = 0
    for t in range(TOTAL):
        kids = np.concatenate(pending.pop(t)) if t in pending else np.array([], dtype=int); intro = rng.poisson(LAM); nk = len(kids)
        ids = np.arange(n_cases, n_cases + nk + intro); par = np.concatenate([kids, -np.ones(intro, dtype=int)]).astype(int)
        gen = -np.ones(len(ids), dtype=int); cha = -np.ones(len(ids), dtype=int)
        if nk:
            garr, carr = np.array(geno), np.array(chain); gen[:nk] = np.where(rng.random(nk) < m, -1, garr[kids]); cha[:nk] = carr[kids]            # a mutation at transmission starts a new genotype; the chain stays the same
        cha[nk:] = np.arange(next_c, next_c + intro); next_c += intro                                                                                    # every introduction starts a new chain
        newg = np.where(gen == -1)[0]; gen[newg] = np.arange(next_g, next_g + len(newg)); next_g += len(newg)
        month += [t] * len(ids); parent += list(par); geno += list(gen); chain += list(cha); n_cases += len(ids)
        n_off = rng.poisson(R, len(ids)); tot = n_off.sum()
        if tot:
            rep = np.repeat(ids, n_off); tm = t + rng.choice(len(w), size=tot, p=w) + 1
            for u in np.unique(tm):
                if u < TOTAL: pending.setdefault(int(u), []).append(rep[tm == u])
    return np.array(month), np.array(parent), np.array(geno), np.array(chain)

def evaluate(month, parent, geno, chain, D, f, rng):
    inwin = np.where((month >= BURN) & (month < BURN + D))[0]; has_inf = parent[inwin] >= 0
    direct = (has_inf & (month[np.where(has_inf, parent[inwin], 0)] >= BURN)).mean()                                  # infector also infected in the window
    chained = (len(inwin) - len(np.unique(chain[inwin]))) / len(inwin)                                                 # cases that share a transmission chain with an earlier case in the window
    n = max(2, int(round(f * len(inwin)))); s = rng.choice(inwin, n, replace=False); cnt = pd.Series(geno[s]).value_counts(); cl = cnt[cnt >= 2]; est = (cl.sum() - len(cl)) / n
    return direct, chained, est, len(inwin)

if __name__ == "__main__":
    rng = np.random.default_rng(1); REPS = 40; Ds = [12, 24, 36, 60]; Fs = [0.2, 0.4, 0.6, 0.8, 1.0]; rows = []
    for rep in range(REPS):
        mo, pa, ge, ch = simulate(rng, 0.03)
        for D in Ds:
            for f in Fs: dr, cn, e, nwin = evaluate(mo, pa, ge, ch, D, f, rng); rows.append(dict(rep=rep, D_months=D, f=f, truth_direct=dr, truth_chain=cn, estimate=e, cases_in_window=nwin))
    S = pd.DataFrame(rows); T = S.groupby(["D_months", "f"]).agg(truth_direct=("truth_direct", "mean"), truth_chain=("truth_chain", "mean"), estimate=("estimate", "mean"), est_sd=("estimate", "std"), cases=("cases_in_window", "mean")).reset_index()
    T["bias_vs_chain"] = T.estimate - T.truth_chain; T["relative_bias_vs_chain_pct"] = 100 * T.bias_vs_chain / T.truth_chain; T["relative_bias_vs_direct_pct"] = 100 * (T.estimate - T.truth_direct) / T.truth_direct
    T.round(4).to_csv("table_design_simulation.csv", index=False); print(f"cases per month about {S.cases_in_window.mean() / S.D_months.mean():.0f}; mutation probability 0.03 per transmission"); print(T.round(3).to_string(index=False))
    rows2 = []
    for mu in (0.0, 0.03, 0.1, 0.3):
        for rep in range(25):
            mo, pa, ge, ch = simulate(rng, mu); dr, cn, e, _ = evaluate(mo, pa, ge, ch, 36, 1.0, rng); rows2.append(dict(m=mu, truth_direct=dr, truth_chain=cn, estimate=e))
    M2 = pd.DataFrame(rows2).groupby("m").agg(truth_direct=("truth_direct", "mean"), truth_chain=("truth_chain", "mean"), estimate=("estimate", "mean"), sd=("estimate", "std")).reset_index(); M2["bias_vs_chain"] = M2.estimate - M2.truth_chain; M2.round(4).to_csv("table_design_mutation.csv", index=False); print(M2.round(3).to_string(index=False))
    fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5)); cols = {12: "#c53030", 24: "#dd6b20", 36: "#2b6cb0", 60: "#2f855a"}
    for D in Ds: x = T[T.D_months == D]; ax[0].plot(x.f, x.estimate, "o-", color=cols[D], label=f"{D} months"); ax[0].axhline(x.truth_chain.iloc[0], color=cols[D], ls=":", lw=1)
    ax[0].set(xlabel="Fraction of cases sampled", ylabel="Estimated recent transmission (n-1 method)", title="Estimate (solid) and chain truth (dotted)", ylim=(0, 0.6)); ax[0].legend(frameon=False, fontsize=8.5, title="Study length", title_fontsize=8.5)
    for D in Ds: x = T[T.D_months == D]; ax[1].plot(x.f, x.relative_bias_vs_chain_pct, "o-", color=cols[D])
    ax[1].axhline(0, color="black", lw=0.8); ax[1].set(xlabel="Fraction of cases sampled", ylabel="Relative bias against the chain truth (%)", title="Bias against sampling fraction")
    xs = np.arange(len(M2)); ax[2].bar(xs - 0.2, M2.truth_chain, 0.4, color="#a0aec0", label="Chain truth"); ax[2].bar(xs + 0.2, M2.estimate, 0.4, color="#2b6cb0", label="Estimate"); ax[2].plot(xs, M2.truth_direct, "k_", ms=18, mew=2, label="Direct truth"); ax[2].set_xticks(xs); ax[2].set_xticklabels([f"{m_:g}" for m_ in M2.m])
    ax[2].set(xlabel="Genotype change per transmission", ylabel="Recent transmission", title="Unstable markers (36 months, all sampled)"); ax[2].legend(frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig("figs1/fig3_design_simulation.png", dpi=200); plt.close(fig); print("figure written")
