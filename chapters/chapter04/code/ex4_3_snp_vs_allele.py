"""
Chapter 4, Hands-on 3 - SNP distance against allele distance, and what recombination does to each.
A pair of isolates is simulated on a Listeria-sized core genome (2.86 Mb aligned, as in the real Pathogen Detection data of Hands-on 4) with a cgMLST scheme of 1,700 loci. Of the SNPs, a fraction (0.6) fall in scheme loci;
the SNP distance counts all of them, the allele distance counts the loci in which at least one SNP fell. Clonal divergence is d point mutations. Recombination adds events, Poisson(lambda) per pair, each placing 5-40 SNPs in
one locus (or in non-locus sequence). Pairs are called 'linked' when the distance is at most 10 under each rule. Output: table_snp_vs_allele.csv, table_snp_allele_calls.csv, Figure 3.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(4); NL, PL, N = 1700, 0.6, 3000

def pair(d, lam):
    loci = set(); snps = 0
    for _ in range(d):
        snps += 1
        if rng.random() < PL: loci.add(int(rng.integers(NL)))
    for _ in range(rng.poisson(lam)):
        m = int(rng.integers(5, 41)); snps += m
        if rng.random() < PL: loci.add(int(rng.integers(NL)))                                   # all m SNPs of the event fall in one locus: one allele change
    return snps, len(loci)

rows = []
for lam in (0.0, 0.25, 0.5, 1.0):
    for d in range(0, 31, 1):
        s, a = zip(*[pair(d, lam) for _ in range(N)]); s, a = np.array(s), np.array(a); rows.append(dict(lam=lam, true_clonal_distance=d, mean_snp=s.mean(), mean_allele=a.mean(), snp_le10=(s <= 10).mean() * 100, allele_le10=(a <= 10).mean() * 100, snp_le5=(s <= 5).mean() * 100, allele_le5=(a <= 5).mean() * 100))
R = pd.DataFrame(rows); R.round(3).to_csv("table_snp_vs_allele.csv", index=False)
for lam in (0.0, 0.5): x = R[R.lam == lam]; print(f"lambda {lam}: mean SNP / mean allele distance at true d = 5, 10, 20, 30:", [(round(x.mean_snp[x.true_clonal_distance == d].iloc[0], 1), round(x.mean_allele[x.true_clonal_distance == d].iloc[0], 1)) for d in (5, 10, 20, 30)])
sc = []
for lam in (0.0, 0.25, 0.5, 1.0):
    x = R[R.lam == lam].set_index("true_clonal_distance")
    sc.append(dict(lam=lam, snp_rule_linked_d5=x.snp_le10[5], allele_rule_linked_d5=x.allele_le10[5], snp_rule_linked_d8=x.snp_le10[8], allele_rule_linked_d8=x.allele_le10[8], snp_rule_linked_d15=x.snp_le10[15], allele_rule_linked_d15=x.allele_le10[15], snp_rule_linked_d20=x.snp_le10[20], allele_rule_linked_d20=x.allele_le10[20]))
C = pd.DataFrame(sc); C.round(1).assign(lam=C.lam).to_csv("table_snp_allele_calls.csv", index=False)       # keep λ exact: round(1) turned 0.25 into 0.2; print("percent of pairs called linked (distance <= 10) by true clonal distance d:"); print(C.round(1).to_string(index=False))
d_equal = R[(R.lam == 0.0)].assign(ratio=lambda t: t.mean_allele / t.mean_snp.replace(0, np.nan)); print(f"without recombination the allele distance is {d_equal.ratio.dropna().mean():.2f} of the SNP distance on average (the share of SNPs in loci is {PL})")
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6)); cols = {0.0: "#2b6cb0", 0.5: "#dd6b20", 1.0: "#c53030"}
for lam, c in cols.items(): x = R[R.lam == lam]; ax[0].plot(x.true_clonal_distance, x.mean_snp, "-", color=c, label=f"SNP distance, {lam:g} events"); ax[0].plot(x.true_clonal_distance, x.mean_allele, "--", color=c, label=f"allele distance, {lam:g} events")
ax[0].axhline(10, color="black", lw=0.8, ls=":"); ax[0].set(xlabel="True clonal distance (point mutations)", ylabel="Observed distance", title="Recombination inflates SNPs more than alleles"); ax[0].legend(frameon=False, fontsize=7.2, loc="upper left")
for lam, c in [(0.0, "#2b6cb0"), (0.5, "#dd6b20")]: x = R[R.lam == lam]; ax[1].plot(x.true_clonal_distance, x.snp_le10, "-", color=c, label=f"SNP <= 10, {lam:g} events"); ax[1].plot(x.true_clonal_distance, x.allele_le10, "--", color=c, label=f"allele <= 10, {lam:g} events")
ax[1].set(xlabel="True clonal distance (point mutations)", ylabel="Pairs called linked (%)", title="The same cut-off, different calls"); ax[1].legend(frameon=False, fontsize=7.6)
x = R[R.lam == 0.0]; ax[2].plot(x.mean_snp, x.mean_allele, color="#2b6cb0", lw=2); ax[2].plot([0, 30], [0, 30], color="black", ls=":", lw=1); ax[2].set(xlabel="Mean SNP distance", ylabel="Mean allele distance", title="Allele distance never exceeds SNP distance")
fig.tight_layout(); fig.savefig("figs4/fig3_snp_vs_allele.png", dpi=200); plt.close(fig); print("figure written")
