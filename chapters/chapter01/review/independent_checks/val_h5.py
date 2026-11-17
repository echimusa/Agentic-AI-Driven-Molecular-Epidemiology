import numpy as np, pandas as pd, sys
sys.path.insert(0, ".")
from epi_sim import *
def pairs_for(o, rho, rng):
    n = o["n"]; dep = depths(o["parent"]); s = np.where(rng.random(n) < rho)[0]
    if len(s) < 3: return []
    g = {i: sample_genome(o, i, rng.uniform(3, 10)) for i in s}; out = []
    for a in range(len(s)):
        for b in range(a + 1, len(s)):
            u, v = s[a], s[b]; out.append((steps_between(o["parent"], dep, u, v), len(g[u] ^ g[v])))
    return out
LAM_E, LAM_21 = 23.7 / 365, 21.05 / 365; print("LAM_DEFAULT in epi_sim:", round(LAM_DEFAULT, 4), "per genome per day =", round(LAM_DEFAULT * 365, 1), "per year")
def stats(lam, seed, nout=20, n=250):
    rng = np.random.default_rng(seed); P = []; Pd = []
    for _ in range(nout):
        o = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=n, min_cases=n, lam=lam); pr = np.array(pairs_for(o, 0.5, rng)); P.append(pr)
    P = np.vstack(P); st, d = P[:, 0], P[:, 1]; direct = st == 1
    row = dict(direct_mean_snp=d[direct].mean(), direct_identical=100 * (d[direct] == 0).mean(), direct_le2=100 * (d[direct] <= 2).mean(), step4_mean=d[st >= 4].mean())
    for t in (0, 1, 2): called = d <= t; row[f"sens_le{t}"] = 100 * (called & direct).sum() / direct.sum(); row[f"ppv_le{t}"] = 100 * (called & direct).sum() / called.sum()
    return row
out = []
for name, lam in (("23.7/yr (Chapter 5)", LAM_E), ("21.05/yr (Chapter 13)", LAM_21)):
    R = pd.DataFrame([stats(lam, 100 + s) for s in range(6)]); out.append((name, R.mean(), R.std()))
for name, m, s in out: print(f"\n{name}: mean (SD across 6 independent seeds x 20 outbreaks of 250 cases, half sequenced)"); print(pd.DataFrame({"mean": m.round(2), "sd": s.round(2)}).T.to_string())
# analytic expectation for a direct pair: mutations accumulated between the parent's sampling and the child's sampling along the two lineages
rng = np.random.default_rng(3); g = rng.gamma(4.0, 15.0 / 4.0, 400000); dp = rng.uniform(3, 10, 400000); dc = rng.uniform(3, 10, 400000)
for nm, lam in (("23.7/yr", LAM_E), ("21.05/yr", LAM_21)): print(f"analytic mean SNPs between a parent and a child, {nm}: {lam * (np.abs(g - dp) + dc).mean():.2f}")
print("share of generation intervals <= 60 days (horizon H):", round((g <= 60).mean(), 4))
# other numbers quoted in the chapter
d = pd.read_csv("typing_pfge_subset.csv.gz", dtype={"country": str, "pfge_xbai": str, "pfge_blni": str}, low_memory=False); q = d[d.pfge_xbai.notna() & d.snp_cluster.notna() & d.serotype.notna()]
print("PFGE isolates:", len(d), "| human & USA with PFGE:", int(((d.sector == "Human") & (d.country == "USA")).sum()), "| human isolates in the 4,811:", int((q.sector == "Human").sum()))
nest = q.groupby("snp_cluster").pfge_xbai.nunique(); print("SNP clusters with >1 PFGE pattern:", int((nest > 1).sum()), "of", q.snp_cluster.nunique()); top = q[q.pfge_xbai == q.pfge_xbai.value_counts().index[0]]; print("commonest PFGE pattern:", top.pfge_xbai.iloc[0], len(top), "isolates;", top.snp_cluster.nunique(), "SNP clusters;", top.serotype.nunique(), "serotypes")
