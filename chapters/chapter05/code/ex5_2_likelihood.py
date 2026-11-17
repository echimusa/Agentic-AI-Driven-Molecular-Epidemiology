"""
Chapter 5, Hands-on 2 - Likelihood from first principles, checked against IQ-TREE.
(1) For two sequences the maximum-likelihood distance under JC69 equals the closed-form Jukes-Cantor distance, which links the formulas of the distance section to the likelihood.
(2) Felsenstein's pruning algorithm is written from scratch for JC69, K80, HKY85 and GTR (time-reversible rate matrix, P(t) = exp(Qt) by eigen-decomposition) and used to fit 8 influenza H3N2 HA sequences on a fixed topology;
the log-likelihoods are compared with those of IQ-TREE for the same topology and models, and AIC and BIC are computed.
(3) The model-selection tables of IQ-TREE for the two full datasets are read back.
Output: table_likelihood_models.csv, table_modelfinder.csv, Figure 2.
"""
import re, subprocess, io
import numpy as np, pandas as pd
from scipy.optimize import minimize, minimize_scalar
from Bio import Phylo, AlignIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo_utils import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
BASES = "ACGT"
def rate_matrix(model, par, pi):
    r = np.ones((4, 4))                                                # exchangeabilities, order A C G T
    if model in ("K80", "HKY"): r[0, 2] = r[2, 0] = r[1, 3] = r[3, 1] = par[0]               # transitions A<->G and C<->T at rate kappa
    if model == "GTR":
        ac, ag, at, cg, ct = par; gt = 1.0
        for (i, j), v in {(0, 1): ac, (0, 2): ag, (0, 3): at, (1, 2): cg, (1, 3): ct, (2, 3): gt}.items(): r[i, j] = r[j, i] = v
    Q = r * pi[None, :]; np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); return Q / -(pi * np.diag(Q)).sum()
def P_of_t(Q, pi, ts):
    s = np.sqrt(pi); S = (s[:, None] * Q) / s[None, :]; S = (S + S.T) / 2; w, U = np.linalg.eigh(S); return [(U * np.exp(w * t)) @ U.T * (1 / s)[:, None] * s[None, :] for t in ts]
def loglik(tree, order, tips, counts, model, par, pi, bl):
    Q = rate_matrix(model, par, pi); Ps = dict(zip([id(c) for c in order], P_of_t(Q, pi, bl))); L = {}
    for cl in tree.get_terminals(): L[id(cl)] = tips[cl.name]
    for cl in [c for c in tree.find_clades(order="postorder") if not c.is_terminal()]:
        v = np.ones_like(next(iter(L.values())))
        for c in cl.clades: v = v * (Ps[id(c)] @ L[id(c)])
        L[id(cl)] = v
    site = pi @ L[id(tree.root)]; return float((counts * np.log(site)).sum())
def fit(tree, tips, counts, pi, model):
    order = [c for c in tree.find_clades() if c is not tree.root]; nb = len(order); npar = {"JC": 0, "K80": 1, "HKY": 1, "GTR": 5}[model]; pi_use = np.full(4, 0.25) if model in ("JC", "K80") else pi
    def f(th):
        bl = np.exp(th[:nb]); par = np.exp(th[nb:]); return -loglik(tree, order, tips, counts, model, par, pi_use, bl)
    th0 = np.concatenate([np.log(np.full(nb, 0.01)), np.zeros(npar)]); best = None
    for start in (0.0, 0.7):
        th = th0.copy(); th[nb:] = start; r = minimize(f, th, method="L-BFGS-B", bounds=[(-12, 2)] * nb + [(-4, 4)] * npar)
        if best is None or r.fun < best.fun: best = r
    return -best.fun, nb, npar, np.exp(best.x[nb:]), np.exp(best.x[:nb])
# ---------------------------------------------------------------- (1) two sequences
p_obs = 0.10; n = 1000; k = int(p_obs * n)
def two_seq_ll(d): pdiff = 0.75 * (1 - np.exp(-4 * d / 3)); return k * np.log(0.25 * pdiff / 0.75 * 1) + (n - k) * np.log(0.25 * (1 - pdiff) / 0.25 * 1) if False else k * np.log(pdiff / 3 * 1) + (n - k) * np.log(1 - pdiff)
ds = np.linspace(0.001, 0.4, 400); ll2 = np.array([two_seq_ll(d) for d in ds]); dml = minimize_scalar(lambda d: -two_seq_ll(d), bounds=(1e-6, 1), method="bounded").x
print(f"two sequences, p = {p_obs}: closed-form JC69 distance {jc(p_obs):.6f}; maximum-likelihood distance {dml:.6f}")
Q = rate_matrix("JC", [], np.full(4, 0.25)); Pt = P_of_t(Q, np.full(4, 0.25), [0.2])[0]; closed = 0.25 + 0.75 * np.exp(-4 * 0.2 / 3); print(f"JC69 P_ii(t = 0.2): eigen-decomposition {Pt[0, 0]:.6f}, closed form {closed:.6f}")
# ---------------------------------------------------------------- (2) 8 taxa vs IQ-TREE
names, X = read_aln("h3n2_ha_aln.fasta"); meta = pd.read_csv("h3n2_meta.csv").set_index("accession"); order_by_year = sorted(names, key=lambda a: (int(meta.year[a]), a)); sel = [order_by_year[i] for i in (0, 6, 12, 18, 24, 30, 36, 44)] if False else [order_by_year[i] for i in (1, 9, 14, 22, 27, 35, 43, 52)]
ix = [names.index(s) for s in sel]; sub = X[ix]; keep = np.all(np.isin(sub, list("ACGT")), axis=0); sub = sub[:, keep]; print(f"8 taxa ({sorted(set(int(meta.year[s]) for s in sel))}), {keep.sum()} gap-free sites of {X.shape[1]}")
with open("sub8.fasta", "w") as f: [f.write(f">{s}\n{''.join(sub[i])}\n") for i, s in enumerate(sel)]
ml = Phylo.read("h3n2.treefile", "newick"); [ml.prune(t) for t in [t.name for t in ml.get_terminals()] if t not in sel]; Phylo.write(ml, "sub8.nwk", "newick")
cols = ["".join(sub[:, j]) for j in range(sub.shape[1])]; pats = pd.Series(cols).value_counts(); counts = pats.values.astype(float); tips = {}
for i, s in enumerate(sel):
    M = np.zeros((4, len(pats)))
    for j, pat in enumerate(pats.index): M[BASES.index(pat[i]), j] = 1
    tips[s] = M
freq = np.array([np.mean(sub == b) for b in BASES]); print("empirical base frequencies A C G T:", freq.round(3).tolist(), "| patterns:", len(pats))
iq_model = {"JC": "JC", "K80": "K2P", "HKY": "HKY+F", "GTR": "GTR+F"}; rows = []; Nsites = int(keep.sum())
for m in ("JC", "K80", "HKY", "GTR"):
    lnl, nb, npar, par, bl = fit(ml, tips, counts, freq, m); kfree = None
    out = subprocess.run(["iqtree2", "-s", "sub8.fasta", "-te", "sub8.nwk", "-m", iq_model[m], "--prefix", f"sub8_{m}", "-redo", "-nt", "1", "-quiet"], capture_output=True, text=True); txt = open(f"sub8_{m}.iqtree").read(); iq = float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", txt).group(1)); kfree = int(re.search(r"Number of free parameters \(#branches \+ #model parameters\): (\d+)", txt).group(1))
    rows.append(dict(model=m, free_parameters=kfree, lnL_pruning=lnl, lnL_iqtree=iq, difference=lnl - iq, AIC=2 * kfree - 2 * lnl, BIC=kfree * np.log(Nsites) - 2 * lnl, kappa=(par[0] if m in ("K80", "HKY") else np.nan)))
T = pd.DataFrame(rows); T["dBIC"] = T.BIC - T.BIC.min(); T["dAIC"] = T.AIC - T.AIC.min(); T.round(4).to_csv("table_likelihood_models.csv", index=False); print(T.round(3).to_string(index=False))
# ---------------------------------------------------------------- (3) ModelFinder tables of the full datasets
mf = []
for d, pre in [("H3N2 HA", "h3n2"), ("Ebola", "ebola")]:
    txt = open(f"{pre}.iqtree").read(); a = txt.index("List of models sorted by BIC scores:"); lines = txt[a:].split("\n")[3:11]
    for ln in lines:
        t = ln.split()
        nums = [float(x) for x in t[1:] if re.fullmatch(r"-?\d+(\.\d+)?(e[+-]?\d+)?", x)]
        if len(nums) >= 7: mf.append(dict(dataset=d, model=t[0], lnL=nums[0], AIC=nums[1], BIC=nums[5]))
MF = pd.DataFrame(mf); MF["dBIC"] = MF.groupby("dataset").BIC.transform(lambda s: s - s.min()); MF.round(2).to_csv("table_modelfinder.csv", index=False); print(MF.groupby("dataset").head(4).round(1).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.6)); ax[0].plot(ds, ll2 - ll2.max(), color="#2b6cb0", lw=2); ax[0].axvline(jc(p_obs), color="#c53030", ls="--"); ax[0].text(jc(p_obs) + 0.008, -40, "JC69 formula\nand ML agree", color="#c53030", fontsize=8.5); ax[0].set(xlabel="Distance d (substitutions per site)", ylabel="Log-likelihood minus maximum", title="Two sequences, 100 of 1,000 sites differ", ylim=(-60, 2))
x = np.arange(len(T)); ax[1].bar(x - 0.2, T.dBIC, 0.4, color="#2b6cb0", label="Delta BIC"); ax[1].bar(x + 0.2, T.dAIC, 0.4, color="#dd6b20", label="Delta AIC"); ax[1].set_xticks(x); ax[1].set_xticklabels(T.model); ax[1].set(ylabel="Difference from the best model", title="8 sequences: AIC and BIC"); ax[1].legend(frameon=False)
for d, c in [("H3N2 HA", "#2b6cb0"), ("Ebola", "#dd6b20")]: s = MF[MF.dataset == d].head(6); ax[2].plot(range(len(s)), s.dBIC, "o-", color=c, label=d)
ax[2].set_xticks(range(6)); ax[2].set_xticklabels([f"{i + 1}" for i in range(6)]); ax[2].set(xlabel="Rank of the model by BIC", ylabel="Delta BIC", title="ModelFinder on the full datasets"); ax[2].legend(frameon=False)
fig.tight_layout(); fig.savefig("figs5/fig2_likelihood.png", dpi=200); plt.close(fig); print("figure written")
