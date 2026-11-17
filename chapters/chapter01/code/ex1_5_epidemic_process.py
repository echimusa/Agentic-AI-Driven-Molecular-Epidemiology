"""
Chapter 1, Hands-on 5 - Simulating the epidemic process behind a genomic data set (transmission, mutation, sampling).
Uses epi_sim.py. (a) How SNP distance relates to the number of transmission steps between two sampled cases; (b) how well a SNP threshold identifies direct transmission, as the sampling fraction rho changes;
(c) how long two importations must have diverged before their cases can be told apart; (d) the effect of superspreading (dispersion k) on transmission and on groups of identical genomes.
Output: table_sim_steps.csv, table_sim_threshold.csv, table_sim_introductions.csv, table_sim_superspreading.csv, figs/fig_sim_results.png.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt
from epi_sim import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(1)
def pairs_for(o, rho, rng):
    n = o["n"]; dep = depths(o["parent"]); s = np.where(rng.random(n) < rho)[0]
    if len(s) < 3: return []
    g = {i: sample_genome(o, i, rng.uniform(3, 10)) for i in s}; out = []
    for a in range(len(s)):
        for b in range(a + 1, len(s)):
            u, v = s[a], s[b]; out.append((steps_between(o["parent"], dep, u, v), len(g[u] ^ g[v])))
    return out
# (a) steps and SNPs with every case sequenced (R = 1.5, k = 1, 250 cases, 30 outbreaks)
def perf(M, t):
    d = M[:, 0] == 1; sel = M[:, 1] <= t; return 100 * (sel & d).sum() / max(d.sum(), 1), 100 * (sel & d).sum() / max(sel.sum(), 1), 100 * d.mean()
REP = 30; full = []
for _ in range(REP): o = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=250, min_cases=80); full += pairs_for(o, 1.0, rng)
A = np.array(full); rows = []
print(f"every case sequenced: {len(A):,} pairs, direct pairs {(A[:, 0] == 1).sum():,}")
for st in (1, 2, 3, 4):
    m = (A[:, 0] == st) if st < 4 else (A[:, 0] >= 4); x = A[m, 1]; rows.append(dict(steps="4 or more" if st == 4 else str(st), pairs=int(m.sum()), mean_snps=x.mean(), p_zero_snps=(x == 0).mean() * 100, p_le2=(x <= 2).mean() * 100, p_ge5=(x >= 5).mean() * 100))
S = pd.DataFrame(rows); S.round(2).to_csv("table_sim_steps.csv", index=False); print(S.round(2).to_string(index=False))
# (b1) sampling fraction at 250 cases
rows = []
for rho in (0.05, 0.1, 0.25, 0.5, 1.0):
    P = []; infector_sampled = []
    for _ in range(REP if rho < 1 else 0):
        o = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=250, min_cases=80); P += pairs_for(o, rho, rng)
    M = A if rho == 1.0 else np.array(P)
    for t in (1, 2): se, pp, dsh = perf(M, t); rows.append(dict(rho_pct=rho * 100, threshold=t, sensitivity=se, ppv=pp, direct_share_of_pairs=dsh, sampled_pairs=len(M)))
Tr = pd.DataFrame(rows); Tr.round(2).to_csv("table_sim_sampling.csv", index=False); print(Tr.round(1).to_string(index=False))
# (b2) outbreak size at 50% sampling
rows = []
for ncap in (50, 100, 250, 500):
    P = []
    for _ in range(REP): o = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=ncap, min_cases=int(0.8 * ncap)); P += pairs_for(o, 0.5, rng)
    M = np.array(P)
    for t in (0, 1, 2): se, pp, dsh = perf(M, t); rows.append(dict(cases=ncap, threshold=t, sensitivity=se, ppv=pp, direct_share_of_pairs=dsh, sampled_pairs=len(M)))
T = pd.DataFrame(rows); T.round(2).to_csv("table_sim_threshold.csv", index=False); print(T.round(1).to_string(index=False))
# (c) importations
rows = []
for tdiv in (0, 7, 14, 30, 60, 90, 180):
    wi, cr = [], []
    for _ in range(40):
        ctr = [0]; anc = frozenset(); m1 = rng.poisson(LAM_DEFAULT * tdiv); m2 = rng.poisson(LAM_DEFAULT * tdiv)
        r1 = frozenset(range(10 ** 6, 10 ** 6 + m1)); r2 = frozenset(range(2 * 10 ** 6, 2 * 10 ** 6 + m2))
        o1 = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=40, min_cases=20, root_genome=r1, mut_counter=ctr); o2 = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=40, min_cases=20, root_genome=r2, mut_counter=ctr)
        g1 = [sample_genome(o1, i, rng.uniform(3, 10)) for i in range(o1["n"])]; g2 = [sample_genome(o2, i, rng.uniform(3, 10)) for i in range(o2["n"])]
        for a in range(min(15, len(g1))):
            for b in range(a + 1, min(15, len(g1))): wi.append(len(g1[a] ^ g1[b]))
            for b in range(min(15, len(g2))): cr.append(len(g1[a] ^ g2[b]))
    wi, cr = np.array(wi), np.array(cr); rows.append(dict(divergence_days=tdiv, within_le2=100 * (wi <= 2).mean(), across_le2=100 * (cr <= 2).mean(), across_le5=100 * (cr <= 5).mean(), across_mean=cr.mean(), within_mean=wi.mean()))
I = pd.DataFrame(rows); I.round(2).to_csv("table_sim_introductions.csv", index=False); print(I.round(1).to_string(index=False))
# (d) superspreading
rows = []
for kk, label in ((0.3, "k = 0.3"), (1.0, "k = 1"), (None, "Poisson")):
    share, grp, nn = [], [], []
    for _ in range(40):
        o = simulate_outbreak(rng, R=1.5, k=kk, n_cap=250, min_cases=60); n = o["n"]; pr = o["processed"]; z = np.bincount(o["parent"][1:], minlength=n)[:pr]; zs = np.sort(z)[::-1]; share.append(zs[: max(1, int(0.2 * pr))].sum() / max(z.sum(), 1))
        g = [sample_genome(o, i, rng.uniform(3, 10)) for i in range(n)]; cnt = {}
        for x in g: cnt[x] = cnt.get(x, 0) + 1
        grp.append(max(cnt.values()) / n); nn.append(n)
    rows.append(dict(offspring=label, top20_share_of_transmission=100 * np.mean(share), largest_identical_group_pct=100 * np.mean(grp), mean_cases=np.mean(nn)))
Sp = pd.DataFrame(rows); Sp.round(2).to_csv("table_sim_superspreading.csv", index=False); print(Sp.round(1).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7)); data = [A[A[:, 0] == s, 1] for s in (1, 2, 3)] + [A[A[:, 0] >= 4, 1]]
bp = ax[0].boxplot(data, tick_labels=["1 (direct)", "2", "3", "4 or more"], showfliers=False, patch_artist=True)
for b in bp["boxes"]: b.set_facecolor("#bee3f8")
ax[0].set(xlabel="Transmission steps between two sampled cases", ylabel="SNP distance", title="SNP distance rises with steps, with overlap")
for t, c in [(1, "#2b6cb0"), (2, "#dd6b20")]:
    s = T[T.threshold == t]; ax[1].plot(s.cases, s.ppv, "o-", color=c, label=f"PPV: share of pairs at most {t} SNP{'s' if t > 1 else ''} apart that are direct"); ax[1].plot(s.cases, s.sensitivity, "s--", color=c, alpha=0.7, label=f"sensitivity: share of direct pairs at most {t} SNP{'s' if t > 1 else ''} apart")
ax[1].set_xscale("log"); ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_xticks([50, 100, 250, 500]); ax[1].set_xticklabels(["50", "100", "250", "500"]); ax[1].set(xlabel="Cases in the outbreak (half sequenced)", ylabel="Percent", title="Close pairs are mostly not direct", ylim=(0, 100)); ax[1].legend(frameon=False, fontsize=6.8, loc="center right")
ax[2].plot(I.divergence_days, I.across_le2, "o-", color="#c53030", label="Between importations: at most 2 SNPs"); ax[2].plot(I.divergence_days, I.across_le5, "o-", color="#dd6b20", label="Between importations: at most 5 SNPs"); ax[2].plot(I.divergence_days, I.within_le2, "s--", color="#2b6cb0", label="Within one outbreak: at most 2 SNPs")
ax[2].set(xlabel="Divergence of source lineages (days)", ylabel="Pairs (%)", title="Dispelling transmission needs divergence", ylim=(0, 100)); ax[2].legend(frameon=False, fontsize=7.4)
fig.tight_layout(); fig.savefig("figs/fig_sim_results.png", dpi=200); plt.close(fig); print("figure written")
