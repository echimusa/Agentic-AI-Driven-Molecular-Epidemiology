"""
Chapter 1, Hands-on 2 - SNP distance, time and the molecular clock on public Ebola virus genomes.
Data: 195 Ebola virus genomes of NCBI BioProject PRJNA257197 sampled in Sierra Leone between 27 May and 28 September 2014 (alignment from Chapter 13).
Question: how well does the number of SNPs separating two genomes tell us how close in time they were sampled, and how long ago they shared an ancestor?
Output: table_snp_distance.csv, table_snp_thresholds.csv, Figure 5 of the chapter (figs1/fig2_snp_distance.png).
"""
import json
import numpy as np, pandas as pd
from Bio import SeqIO
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
recs = list(SeqIO.parse("ebola_aln.fasta", "fasta")); ids = [r.id for r in recs]; L = len(recs[0].seq); A = np.array([list(str(r.seq).upper()) for r in recs]); code = {"A": 1, "C": 2, "G": 3, "T": 4}
M = np.zeros(A.shape, dtype=np.uint8)
for k, v in code.items(): M[A == k] = v
valid = M > 0; dates = pd.read_csv("ebola_dates.csv").set_index("acc").date; day = pd.to_datetime(dates.reindex(ids)).values.astype("datetime64[D]").astype(int); n = len(ids)
D = np.zeros((n, n), dtype=int)
for i in range(n): D[i] = ((M[i] != M) & valid[i] & valid).sum(axis=1)
iu = np.triu_indices(n, 1); snps = D[iu]; gap = np.abs(day[iu[0]] - day[iu[1]]); print(f"{n} genomes, alignment length {L:,}; {len(snps):,} pairs; median pairwise distance {np.median(snps):.0f} SNPs (range {snps.min()}-{snps.max()}); pairs within 2 SNPs: {(snps <= 2).mean():.1%}; identical pairs: {(snps == 0).mean():.1%}")
rho, p = stats.spearmanr(gap, snps); print(f"Spearman correlation of SNP distance with days between samples: {rho:.3f}")
bins = [(0, 7), (8, 30), (31, 60), (61, 90), (91, 150)]; rows = []
for a, b in bins:
    s = snps[(gap >= a) & (gap <= b)]; rows.append(dict(days_apart=f"{a}-{b}", pairs=len(s), median_snps=np.median(s), q25=np.percentile(s, 25), q75=np.percentile(s, 75), within_2_snps_pct=(s <= 2).mean() * 100, within_5_snps_pct=(s <= 5).mean() * 100))
T = pd.DataFrame(rows); T.round(2).to_csv("table_snp_distance.csv", index=False); print(T.round(2).to_string(index=False))
# SNP distance as a classifier of 'sampled within 30 days of each other'
y = (gap <= 30).astype(int); order = np.argsort(snps); 
def auc(score, y):   # probability that a random close pair has a smaller distance than a random distant pair
    r = stats.rankdata(-score); n1, n0 = y.sum(), (1 - y).sum(); return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
A_ = auc(snps, y); print(f"AUC of SNP distance for 'sampled within 30 days': {A_:.3f}")
thr = []
for t in (0, 1, 2, 3, 5, 10):
    pred = snps <= t; thr.append(dict(snp_threshold=t, pairs_called_close=int(pred.sum()), sensitivity_pct=(pred & (y == 1)).sum() / y.sum() * 100, specificity_pct=(~pred & (y == 0)).sum() / (1 - y).sum() * 100, ppv_pct=((pred & (y == 1)).sum() / max(pred.sum(), 1)) * 100))
TT = pd.DataFrame(thr); TT.round(1).to_csv("table_snp_thresholds.csv", index=False); print(TT.round(1).to_string(index=False))
rate_site_year = 1.11e-3; per_genome_year = rate_site_year * L; per_day = per_genome_year / 365.25; print(f"clock rate {rate_site_year:.2e} per site per year (Chapter 13) = {per_genome_year:.1f} SNPs per genome per year = {per_day:.4f} per day")
print("expected time to the common ancestor of two genomes d SNPs apart (days) = d / (2 x rate):", {d_: round(d_ / (2 * per_day), 1) for d_ in (1, 2, 5, 10)})
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5)); mx = int(np.percentile(snps, 99))
ax[0].hist(snps[gap <= 14], bins=np.arange(0, mx + 2), density=True, alpha=0.65, color="#2b6cb0", label=f"Sampled within 14 days (n={int((gap <= 14).sum()):,})"); ax[0].hist(snps[gap > 60], bins=np.arange(0, mx + 2), density=True, alpha=0.55, color="#dd6b20", label=f"More than 60 days apart (n={int((gap > 60).sum()):,})")
ax[0].set(xlabel="SNPs between two genomes", ylabel="Share of pairs", title="Genomes close in time are close in sequence"); ax[0].legend(frameon=False, fontsize=8.5)
data = [snps[(gap >= a) & (gap <= b)] for a, b in bins]; bp = ax[1].boxplot(data, tick_labels=[f"{a}-{b}" for a, b in bins], patch_artist=True, showfliers=False)
for b_ in bp["boxes"]: b_.set(facecolor="#bee3f8")
for m_ in bp["medians"]: m_.set(color="black")
ax[1].set(xlabel="Days between the two samples", ylabel="SNPs between the two genomes", title=f"SNP distance grows with time (rho = {rho:.2f})")
fpr, tpr = [], []
for t in range(0, mx + 1): pred = snps <= t; tpr.append((pred & (y == 1)).sum() / y.sum()); fpr.append((pred & (y == 0)).sum() / (1 - y).sum())
ax[2].plot(fpr, tpr, color="#2b6cb0", lw=2); ax[2].plot([0, 1], [0, 1], color="black", lw=0.8, ls="--")
for t in (0, 2, 5): ax[2].scatter([fpr[t]], [tpr[t]], color="#c53030", zorder=3); ax[2].annotate(f"{t} SNPs", (fpr[t], tpr[t]), xytext=(6, -10), textcoords="offset points", fontsize=8.5)
ax[2].set(xlabel="1 - specificity", ylabel="Sensitivity", title=f"A SNP threshold as a test (AUC {A_:.2f})")
fig.tight_layout(); fig.savefig("figs1/fig2_snp_distance.png", dpi=200); plt.close(fig); print("figure written")
