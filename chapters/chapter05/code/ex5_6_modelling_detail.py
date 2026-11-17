"""
Chapter 5, Hands-on 5 - Modelling in detail: rate variation, codon-position partitions and model adequacy (60 human H3N2 HA coding sequences, 1,701 sites, ML topology fixed).
Part 1: a ladder of 16 substitution models (JC, K2P, HKY+F, GTR+F, each with no rate variation, +I, +G4 and +I+G4), fitted in IQ-TREE on the fixed ML topology; AIC, BIC, gamma shape and invariant-site fraction.
Part 2: the same gene fitted with one HKY+F+G4 model, and with a separate HKY+F+G4 for each codon position under proportional branch lengths (-spp); relative evolutionary rate and constant sites of each codon position.
Part 3: model adequacy by parametric simulation. Alignments are simulated 200 times along the fitted tree under (a) JC with equal rates and (b) HKY+F+G4; four statistics of the real alignment are compared with their simulated
distributions: share of constant sites, number of parsimony-informative sites, mean transition share among pairwise differences, and the variance of GC content among sequences.
Needs: IQ-TREE 2 on the path, h3n2_cds.fasta (written by ex5_6a_extract_cds.py), h3n2_cp1-3.fasta, parts.nex, h3n2.treefile. Output: table_model_ladder.csv, table_codon_positions.csv, table_adequacy.csv, Figure 5.
"""
import re, subprocess
import numpy as np, pandas as pd
from scipy.special import gammainc, gammaincinv
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo_utils import read_aln, encode, pairwise
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
def iq(prefix, args):
    subprocess.run(["iqtree2", "-te", "h3n2.treefile", "--prefix", prefix, "-redo", "-nt", "1", "-quiet"] + args, capture_output=True, text=True); return open(f"{prefix}.iqtree").read()
def num(pat, txt, cast=float, default=np.nan):
    m = re.search(pat, txt); return cast(m.group(1)) if m else default
def summary(txt):
    return dict(lnL=num(r"Log-likelihood of the tree: (-?[\d.]+)", txt), k=num(r"Number of free parameters \(#branches \+ #model parameters\): (\d+)", txt, int), AIC=num(r"Akaike information criterion \(AIC\) score: ([\d.]+)", txt), BIC=num(r"Bayesian information criterion \(BIC\) score: ([\d.]+)", txt),
                alpha=num(r"Gamma shape alpha: ([\d.]+)", txt), pinv=num(r"Proportion of invariable sites: ([\d.]+)", txt), treelen=num(r"Total tree length \(sum of branch lengths\): ([\d.]+)", txt), const_pct=num(r"Number of constant sites: \d+ \(= ([\d.]+)%", txt), kappa=num(r"A-G: ([\d.]+)", txt))
# ---------------------------------------------------------------- part 1: model ladder
rows = []
for base in ("JC", "K2P", "HKY+F", "GTR+F"):
    for suf in ("", "+I", "+G4", "+I+G4"):
        m = base + suf; s = summary(iq("lad_" + re.sub(r"\W", "_", m), ["-s", "h3n2_cds.fasta", "-m", m])); s["model"] = m; rows.append(s)
LD = pd.DataFrame(rows)[["model", "k", "lnL", "AIC", "BIC", "alpha", "pinv", "kappa"]]; LD["dAIC"] = LD.AIC - LD.AIC.min(); LD["dBIC"] = LD.BIC - LD.BIC.min(); LD.round(3).to_csv("table_model_ladder.csv", index=False); print(LD.round(2).to_string(index=False))
# ---------------------------------------------------------------- part 2: codon positions
whole = summary(iq("cds_whole", ["-s", "h3n2_cds.fasta", "-m", "HKY+F+G4"])); part = summary(iq("cds_part", ["-s", "h3n2_cds.fasta", "-spp", "parts.nex"])); print(f"\nwhole gene HKY+F+G4: lnL {whole['lnL']:.2f}, k {whole['k']}, AIC {whole['AIC']:.1f}, BIC {whole['BIC']:.1f}, alpha {whole['alpha']:.3f}, kappa {whole['kappa']:.2f}")
print(f"codon positions (-spp):  lnL {part['lnL']:.2f}, k {part['k']}, AIC {part['AIC']:.1f}, BIC {part['BIC']:.1f}")
cp = []
for p in (1, 2, 3): s = summary(iq(f"cp{p}", ["-s", f"h3n2_cp{p}.fasta", "-m", "HKY+F+G4"])); s["position"] = p; cp.append(s)
CP = pd.DataFrame(cp)[["position", "treelen", "alpha", "kappa", "const_pct", "lnL"]]; CP["relative_rate"] = CP.treelen / CP.treelen.mean(); CP.round(4).to_csv("table_codon_positions.csv", index=False); print(CP.round(3).to_string(index=False))
# ---------------------------------------------------------------- part 3: adequacy
names, X = read_aln("h3n2_cds.fasta"); E = encode(X); n, L = E.shape
def stats(E):
    const = np.mean([len(set(E[:, j][E[:, j] >= 0])) <= 1 for j in range(E.shape[1])]); inf = 0
    for j in range(E.shape[1]):
        c = np.bincount(E[:, j][E[:, j] >= 0], minlength=4); inf += (np.sort(c)[-2] >= 2)
    p, P, Q = pairwise(E); iu = np.triu_indices(E.shape[0], 1); d = P[iu] + Q[iu]; ts = float((P[iu].sum()) / max(d.sum(), 1e-12)); gc = np.array([np.mean((E[i][E[i] >= 0] == 1) | (E[i][E[i] >= 0] == 2)) for i in range(E.shape[0])])
    return dict(constant_sites=const, informative_sites=inf, transition_share=ts, gc_variance=float(gc.var()))
obs = stats(E)
BASES = "AGCT"                                   # the code of phylo_utils: A 0, G 1, C 2, T 3
def rate_matrix(kappa, pi):
    r = np.ones((4, 4)); r[0, 1] = r[1, 0] = r[2, 3] = r[3, 2] = kappa; Q = r * pi[None, :]; np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); return Q / -(pi * np.diag(Q)).sum()
def P_of_t(Q, pi, t):
    s = np.sqrt(pi); S = (s[:, None] * Q) / s[None, :]; S = (S + S.T) / 2; w, U = np.linalg.eigh(S); return (U * np.exp(w * t)) @ U.T * (1 / s)[:, None] * s[None, :]
def gamma_rates(a, K=4):
    q = gammaincinv(a, np.arange(0, K + 1) / K); q[-1] = np.inf; cdf = gammainc(a + 1, q * 1.0); return K * np.diff(cdf)       # category means of a Gamma(a, a) distribution
def fitted(prefix):
    txt = open(prefix + ".iqtree").read(); pi_t = np.array([num(rf"pi\({b}\) = ([\d.]+)", txt) for b in "ACGT"]); pi_t = np.where(np.isnan(pi_t), 0.25, pi_t); tr = Phylo.read(prefix + ".treefile", "newick")
    return tr, pi_t[[0, 2, 1, 3]], num(r"A-G: ([\d.]+)", txt, default=1.0), num(r"Gamma shape alpha: ([\d.]+)", txt, default=np.nan)      # reorder A C G T -> A G C T
def simulate(tr, pi, kappa, alpha, rng):
    Q = rate_matrix(kappa, pi); rates = np.ones(L) if np.isnan(alpha) else None
    if rates is None: gr = gamma_rates(alpha); cat = rng.integers(4, size=L); rates = gr[cat]
    cats = np.unique(rates); state = {tr.root: rng.choice(4, L, p=pi)}; seqs = {}
    for cl in tr.find_clades(order="preorder"):
        for c in cl.clades:
            t = max(c.branch_length or 0.0, 0.0); child = np.empty(L, dtype=int)
            for r in cats:
                idx = np.where(rates == r)[0]; cum = P_of_t(Q, pi, t * r).cumsum(1); u = rng.random(len(idx)); child[idx] = (u[:, None] > cum[state[cl][idx]]).sum(1).clip(0, 3)
            state[c] = child
    tips = {c.name: state[c] for c in tr.get_terminals()}; return np.array([tips[nm] for nm in names])
rng = np.random.default_rng(8); res = {}; S = 200
for label, pre in (("JC (equal rates)", "lad_JC"), ("HKY+F+G4", "lad_HKY_F_G4")):
    tr, pi, kappa, alpha = fitted(pre); sims = [stats(simulate(tr, pi, kappa, alpha, rng)) for _ in range(S)]; res[label] = pd.DataFrame(sims)
rowsA = []
for label, D in res.items():
    for k in obs: lo, hi = D[k].quantile([0.025, 0.975]); p2 = 2 * min((D[k] <= obs[k]).mean(), (D[k] >= obs[k]).mean()); rowsA.append(dict(model=label, statistic=k, observed=obs[k], sim_mean=D[k].mean(), lo=lo, hi=hi, p_two_sided=min(p2, 1.0), adequate=bool(lo <= obs[k] <= hi)))
AD = pd.DataFrame(rowsA); AD.round(4).to_csv("table_adequacy.csv", index=False); print(AD.round(4).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7)); sel = ["JC", "JC+G4", "K2P+G4", "HKY+F", "HKY+F+G4", "HKY+F+I+G4", "GTR+F", "GTR+F+G4", "GTR+F+I+G4"]; Ls = LD.set_index("model").loc[sel]; x = np.arange(len(sel))
ax[0].bar(x - 0.2, Ls.dBIC, 0.4, color="#2b6cb0", label="Delta BIC"); ax[0].bar(x + 0.2, Ls.dAIC, 0.4, color="#dd6b20", label="Delta AIC"); ax[0].set_xticks(x); ax[0].set_xticklabels(sel, rotation=60, ha="right", fontsize=8); ax[0].set_yscale("symlog", linthresh=10); ax[0].set(ylabel="Difference from the best model", title="Model ladder (H3N2 HA coding region)"); ax[0].legend(frameon=False)
ax[1].bar(CP.position, CP.relative_rate, color=["#2b6cb0", "#2f855a", "#dd6b20"]); ax[1].set_xticks([1, 2, 3]); ax[1].set_xticklabels(["1st", "2nd", "3rd"]); ax[1].axhline(1, color="black", lw=0.8, ls=":"); ax[1].set(xlabel="Codon position", ylabel="Relative evolutionary rate (tree length / mean)", title="Codon positions evolve at different rates")
for label, c in (("JC (equal rates)", "#c53030"), ("HKY+F+G4", "#2b6cb0")): ax[2].hist(res[label].constant_sites, bins=20, alpha=0.55, color=c, label=label)
ax[2].axvline(obs["constant_sites"], color="black", lw=2); ax[2].set(xlabel="Share of constant sites", ylabel="Simulated alignments", title="Adequacy: the black line is the real alignment"); ax[2].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("figs5/fig5_modelling_detail.png", dpi=200); plt.close(fig); print("figure written")
