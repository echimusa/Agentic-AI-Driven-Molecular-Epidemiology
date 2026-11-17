"""
Chapter 2, Hands-on 3 - Population stratification: false positives, adjustment and design alternatives.
Starts from the uploaded Chapter_2_3.py (600 people of ancestry A and 400 of ancestry B; a SNP with frequency 0.10 in A and 0.45 in B; a higher disease rate in B from an unmeasured environmental factor;
no causal SNP effect) and repeats it 1,000 times to measure the type I error of four analyses:
  naive logistic regression; adjustment for two principal components of 100 background SNPs; restriction to one ancestry; and a stratified (Mantel-Haenszel-like) analysis within ancestry.
It then varies how strongly ancestry drives disease risk, and shows that adjustment with too few ancestry-informative markers fails.
Output: table_stratification.csv, table_stratification_strength.csv, Figure 3.
"""
import numpy as np, pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
NA, NB, REPS = 600, 400, 1000

def simulate(rng, riskA=-2.5, riskB=-0.5, n_informative=30, n_background=100):
    anc = np.r_[np.zeros(NA, int), np.ones(NB, int)]; n = len(anc); p = np.where(anc == 0, 0.10, 0.45); snp = rng.binomial(2, p)
    prob = 1 / (1 + np.exp(-np.where(anc == 0, riskA, riskB))); disease = (rng.uniform(size=n) < prob).astype(int)          # true SNP effect is zero
    G = rng.binomial(2, 0.2, (n, n_background)); G[anc == 1, :n_informative] = rng.binomial(2, 0.6, (int(anc.sum()), n_informative))
    return anc, snp, disease, G

def pval(y, X):
    try: return sm.Logit(y, sm.add_constant(X)).fit(disp=0).pvalues[1]
    except Exception: return np.nan

def analyses(rng, **kw):
    anc, snp, y, G = simulate(rng, **kw); pcs = PCA(n_components=2).fit_transform(G)
    p_naive = pval(y, snp); p_pc = pval(y, np.column_stack([snp, pcs])); m = anc == 0; p_restrict = pval(y[m], snp[m])
    p_strat = pval(y, np.column_stack([snp, anc]))                                           # ancestry label as a covariate: the ideal adjustment
    return dict(naive=p_naive, pcs=p_pc, restrict=p_restrict, labelled=p_strat), snp[anc == 0].mean() / 2, snp[anc == 1].mean() / 2

rng = np.random.default_rng(99); res = []
for r in range(REPS): d, fa, fb = analyses(rng); res.append(d)
P = pd.DataFrame(res); names = {"naive": "Naive", "pcs": "Adjusted for 2 PCs", "labelled": "Adjusted for ancestry label", "restrict": "Restricted to ancestry A"}
T = pd.DataFrame([dict(analysis=names[k], false_positive_rate_pct=(P[k] < 0.05).mean() * 100, median_p=P[k].median(), lambda_gc=np.median(np.nan_to_num(__import__("scipy.stats", fromlist=["chi2"]).chi2.isf(P[k].dropna(), 1))) / 0.4549) for k in ["naive", "pcs", "labelled", "restrict"]])
T.round(3).to_csv("table_stratification.csv", index=False); print(f"{REPS} simulated studies of {NA + NB} people; true SNP effect = 0; frequency {fa:.2f} in A and {fb:.2f} in B (last replicate)"); print(T.round(3).to_string(index=False))
# one replicate with the odds ratios, as in the uploaded script
x = [analyses(np.random.default_rng(99))]
anc, snp, y, G = simulate(np.random.default_rng(99)); pcs = PCA(n_components=2).fit_transform(G); naive = sm.Logit(y, sm.add_constant(snp)).fit(disp=0); adj = sm.Logit(y, sm.add_constant(np.column_stack([snp, pcs]))).fit(disp=0)
print(f"one study (seed 99): naive OR {np.exp(naive.params[1]):.2f}, P = {naive.pvalues[1]:.2g}; PC-adjusted OR {np.exp(adj.params[1]):.2f}, P = {adj.pvalues[1]:.2g}; PC1 vs ancestry correlation {abs(np.corrcoef(pcs[:, 0], anc)[0, 1]):.3f}")
# how strongly ancestry must drive risk, and how many informative markers PCs need
rows = []
for gap in (0.0, 0.5, 1.0, 1.5, 2.0):
    pp = [analyses(rng, riskA=-2.5, riskB=-2.5 + gap)[0] for _ in range(300)]; Q = pd.DataFrame(pp); rows.append(dict(setting="risk gap (log odds)", value=gap, naive=(Q.naive < 0.05).mean() * 100, pcs=(Q.pcs < 0.05).mean() * 100))
for k in (0, 2, 5, 10, 30):
    pp = [analyses(rng, n_informative=k)[0] for _ in range(300)]; Q = pd.DataFrame(pp); rows.append(dict(setting="informative background SNPs", value=k, naive=(Q.naive < 0.05).mean() * 100, pcs=(Q.pcs < 0.05).mean() * 100))
S = pd.DataFrame(rows); S.round(2).to_csv("table_stratification_strength.csv", index=False); print(S.round(1).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5), gridspec_kw={"width_ratios": [0.85, 1, 1]}); cols = ["#c53030", "#2b6cb0", "#2f855a", "#805ad5"]
ax[0].bar(range(4), T.false_positive_rate_pct, color=cols); ax[0].axhline(5, color="black", ls="--", lw=1); ax[0].set_xticks(range(4)); ax[0].set_xticklabels(["Naive", "2 PCs", "Ancestry\nlabel", "Restricted"], fontsize=9); ax[0].set(ylabel="False positive rate at P < 0.05 (%)", title="Type I error over 1,000 studies")
for i, v in enumerate(T.false_positive_rate_pct): ax[0].text(i, v + 2, f"{v:.1f}", ha="center", fontsize=9)
g = S[S.setting == "risk gap (log odds)"]; ax[1].plot(g.value, g.naive, "o-", color="#c53030", label="Naive"); ax[1].plot(g.value, g.pcs, "o-", color="#2b6cb0", label="2 PCs"); ax[1].axhline(5, color="black", ls="--", lw=1); ax[1].set(xlabel="Extra disease risk in ancestry B (log odds)", ylabel="False positive rate (%)", title="The bias grows with the confounder"); ax[1].legend(frameon=False, fontsize=8.5)
k = S[S.setting == "informative background SNPs"]; ax[2].plot(k.value, k.pcs, "o-", color="#2b6cb0", label="Adjusted for 2 PCs"); ax[2].plot(k.value, k.naive, "o-", color="#c53030", label="Naive"); ax[2].axhline(5, color="black", ls="--", lw=1); ax[2].set(xlabel="Background SNPs that differ by ancestry (of 100)", ylabel="False positive rate (%)", title="PCs need informative markers"); ax[2].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("figs2/fig3_stratification.png", dpi=200); plt.close(fig); print("figure written")
