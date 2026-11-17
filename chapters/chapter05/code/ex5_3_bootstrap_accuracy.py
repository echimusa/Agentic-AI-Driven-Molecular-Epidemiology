"""
Chapter 5, Hands-on 3 - Support values and accuracy.
Part A (real data): the nonparametric bootstrap written from scratch (resample alignment columns, rebuild a neighbour-joining tree, count how often each split of the ML tree recurs, using site patterns for speed) is compared with the
support values of IQ-TREE (ultrafast bootstrap for the H3N2 tree, SH-aLRT for the Ebola tree). Branches too short to carry a substitution are counted.
Part B (simulation): 30-taxon coalescent trees are simulated with HKY sequence evolution (kappa 6) at two divergence levels and four alignment lengths, trees are rebuilt with NJ and FastTree and compared with the truth by the normalised
Robinson-Foulds distance, and the bootstrap support of NJ splits is checked for calibration (is a split with 95% support correct 95% of the time?).
Run: python ex5_3_bootstrap_accuracy.py A | B. Output: table_support_real.csv, table_accuracy_sim.csv, table_calibration_sim.csv, Figure 3.
"""
import sys, subprocess, io, tempfile, os
import numpy as np, pandas as pd
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo_utils import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})

def nj_splits(D):
    n = D.shape[0]; F = np.zeros((2 * n, 2 * n)); F[:n, :n] = D; active = list(range(n)); mem = {i: frozenset([i]) for i in range(n)}; out = set(); nxt = n
    while len(active) > 3:
        m = len(active); sub = F[np.ix_(active, active)]; r = sub.sum(1); Q = (m - 2) * sub - r[:, None] - r[None, :]; np.fill_diagonal(Q, np.inf); i, j = np.unravel_index(np.argmin(Q), Q.shape); a, b = active[i], active[j]
        for k in active:
            if k not in (a, b): F[nxt, k] = F[k, nxt] = (F[a, k] + F[b, k] - F[a, b]) / 2
        s = mem[a] | mem[b]; mem[nxt] = s; c = s if 0 not in s else frozenset(range(n)) - s
        if 1 < len(c) < n - 1: out.add(c)
        active = [k for k in active if k not in (a, b)] + [nxt]; nxt += 1
    return out

def jc_from_patterns(E, w):
    """JC69 distance matrix from encoded sites E (taxa x patterns) weighted by w; sites with a gap or ambiguity in either sequence are skipped for that pair."""
    n = E.shape[0]; same = np.zeros((n, n)); valid = np.zeros((n, n)); V = (E >= 0).astype(float)
    for b in range(4): Z = (E == b).astype(float); same += (Z * w) @ Z.T
    valid = (V * w) @ V.T; p = 1 - same / np.maximum(valid, 1); return np.clip(-0.75 * np.log(np.maximum(1 - 4.0 / 3.0 * np.minimum(p, 0.74), 1e-12)), 0, None)

def part_A():
    rows = []; rng = np.random.default_rng(11)
    for name, aln, tf, lab in [("H3N2 HA", "h3n2_ha_aln.fasta", "h3n2.treefile", "UFBoot"), ("Ebola", "ebola_aln.fasta", "ebola.treefile", "SH-aLRT")]:
        names, X = read_aln(aln); E = encode(X); L = X.shape[1]; pat, inv, cnt = np.unique(E.T, axis=0, return_inverse=True, return_counts=True); Ep = pat.T; ml = Phylo.read(tf, "newick"); idx = {t: i for i, t in enumerate(names)}; full = frozenset(range(len(names)))
        info = {}; n_internal = 0; short = 0
        for cl in ml.get_nonterminals():
            if cl is ml.root: continue
            n_internal += 1
            if cl.branch_length is not None and cl.branch_length < 0.5 / L: short += 1; continue
            s = frozenset(idx[t.name] for t in cl.get_terminals()); s = full - s if 0 in s else s
            if 1 < len(s) < len(names) - 1:
                try: info[s] = float(cl.confidence) if cl.confidence is not None else (float(cl.name) if cl.name else np.nan)
                except Exception: info[s] = np.nan
        B = 100; hits = {s: 0 for s in info}
        for _ in range(B):
            w = rng.multinomial(L, cnt / cnt.sum()).astype(float); sp = nj_splits(jc_from_patterns(Ep, w))
            for s in info:
                if s in sp: hits[s] += 1
        nb = np.array([hits[s] for s in info]) / B * 100; iq = np.array([info[s] for s in info])
        rows.append(dict(dataset=name, internal_branches=n_internal, branches_shorter_than_one_substitution=short, resolved_splits=len(info), iqtree_support_type=lab, iqtree_ge95_pct=100 * np.nanmean(iq >= 95), iqtree_ge80_pct=100 * np.nanmean(iq >= 80), iqtree_lt50_pct=100 * np.nanmean(iq < 50), nj_boot_ge95_pct=100 * (nb >= 95).mean(), nj_boot_ge70_pct=100 * (nb >= 70).mean(), nj_boot_lt50_pct=100 * (nb < 50).mean(), corr_nj_vs_iqtree=np.corrcoef(nb[~np.isnan(iq)], iq[~np.isnan(iq)])[0, 1]))
        np.save(f"support_{name.split()[0]}.npy", np.vstack([nb, iq]))
    R = pd.DataFrame(rows); R.round(2).to_csv("table_support_real.csv", index=False); print(R.round(2).T.to_string())

HKY_K = 6.0; PI = np.array([0.3, 0.2, 0.2, 0.3])
def hky_P(t):
    Rm = np.ones((4, 4)); Rm[0, 2] = Rm[2, 0] = Rm[1, 3] = Rm[3, 1] = HKY_K; Q = Rm * PI[None, :]; np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); Q = Q / -(PI * np.diag(Q)).sum()
    s = np.sqrt(PI); S = (s[:, None] * Q) / s[None, :]; S = (S + S.T) / 2; w, U = np.linalg.eigh(S); return (U * np.exp(w * t)) @ U.T * (1 / s)[:, None] * s[None, :]
def sim_tree(n, height, rng):
    """Kingman coalescent tree scaled to a given expected height; returns (children dict, branch length dict, root id)."""
    nodes = list(range(n)); t = 0.0; kids = {}; age = {i: 0.0 for i in nodes}; nxt = n; k = n
    while len(nodes) > 1:
        t += rng.exponential(1.0 / (k * (k - 1) / 2)); a, b = rng.choice(len(nodes), 2, replace=False); A, Bn = nodes[a], nodes[b]; kids[nxt] = (A, Bn); age[nxt] = t; nodes = [x for x in nodes if x not in (A, Bn)] + [nxt]; nxt += 1; k -= 1
    scale = height / (2 * (1 - 1 / n)); bl = {c: (age[p] - age[c]) * scale for p, (a, b) in kids.items() for c in (a, b)}; return kids, bl, nxt - 1
def sim_aln(n, height, L, rng):
    kids, bl, root = sim_tree(n, height, rng); seqs = {root: rng.choice(4, L, p=PI)}; true = set(); mem = {}
    def desc(x): return frozenset([x]) if x < n else desc(kids[x][0]) | desc(kids[x][1])
    for p in sorted(kids, reverse=True):
        for c in kids[p]:
            P = hky_P(bl[c]); cum = P.cumsum(1); u = rng.random(L); seqs[c] = (u[:, None] > cum[seqs[p]]).sum(1).clip(0, 3)
    for x in kids:
        if x != root: s = desc(x); c = s if 0 not in s else frozenset(range(n)) - s; (1 < len(c) < n - 1) and true.add(c)
    return np.array([seqs[i] for i in range(n)]), true
def fasttree_splits(E, n):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        for i in range(n): f.write(f">t{i}\n" + "".join("ACGT"[b] for b in E[i]) + "\n")
        fn = f.name
    out = subprocess.run(["FastTree", "-nt", "-gtr", "-quiet", "-nosupport", fn], capture_output=True, text=True).stdout; os.unlink(fn); tr = Phylo.read(io.StringIO(out), "newick"); idx = {f"t{i}": i for i in range(n)}; res = set()
    for cl in tr.get_nonterminals():
        if cl is tr.root: continue
        s = frozenset(idx[t.name] for t in cl.get_terminals()); c = s if 0 not in s else frozenset(range(n)) - s
        if 1 < len(c) < n - 1: res.add(c)
    return res
def part_B():
    rng = np.random.default_rng(21); n = 30; rows = []; cal = []
    for H in (0.002, 0.02):
        for L in (250, 1000, 4000, 16000):
            rn, rf_ft, inf = [], [], []
            for rep in range(20):
                E, true = sim_aln(n, H, L, rng); D = jc_from_patterns(E.T, np.ones(E.shape[1])) if False else jc_from_patterns(E, np.ones(L)); sj = nj_splits(D); sf = fasttree_splits(E, n)
                rn.append(rf(sj, true)); rf_ft.append(rf(sf, true)); inf.append(int(sum(len(set(E[:, j])) > 1 and sorted(np.bincount(E[:, j], minlength=4))[-2] >= 2 for j in range(L))))
                if L in (250, 1000) and rep < 20:
                    pat, cnt = np.unique(E, axis=1, return_counts=True); hit = {s: 0 for s in sj}
                    for _ in range(100):
                        w = rng.multinomial(L, cnt / cnt.sum()).astype(float); sp = nj_splits(jc_from_patterns(pat, w))
                        for s in sj:
                            if s in sp: hit[s] += 1
                    for s in sj: cal.append(dict(height=H, sites=L, support=hit[s], correct=s in true))
            rows.append(dict(tree_height=H, sites=L, informative_sites=np.mean(inf), nj_rf=np.mean(rn), nj_rf_sd=np.std(rn), fasttree_rf=np.mean(rf_ft), fasttree_rf_sd=np.std(rf_ft)))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).round(4).to_csv("table_accuracy_sim.csv", index=False); C = pd.DataFrame(cal); C["bin"] = pd.cut(C.support, [-1, 49, 69, 94, 100], labels=["<50", "50-69", "70-94", "95-100"]); G = C.groupby(["height", "bin"], observed=True).agg(splits=("correct", "size"), correct_pct=("correct", lambda s: 100 * s.mean())).reset_index(); G.round(2).to_csv("table_calibration_sim.csv", index=False); print(G.round(1).to_string(index=False))
def part_C():
    """Check the calibration result (are the 95% splits at low divergence arbitrary resolutions of ties?) and draw Figure 3."""
    rng = np.random.default_rng(33); n = 30; rows = []
    for rep in range(20):
        E, true = sim_aln(n, 0.002, 1000, rng); sj = nj_splits(jc_from_patterns(E, np.ones(1000))); pat, cnt = np.unique(E, axis=1, return_counts=True); hit = {s: 0 for s in sj}
        for _ in range(100):
            w = rng.multinomial(1000, cnt / cnt.sum()).astype(float); sp = nj_splits(jc_from_patterns(pat, w))
            for s in sj:
                if s in sp: hit[s] += 1
        cols = [frozenset(np.where(E[:, j] == E[0, j])[0]) for j in range(E.shape[1]) if len(set(E[:, j])) == 2]; colsplits = set()
        for c in cols:
            c = c if 0 not in c else frozenset(range(n)) - c
            colsplits.add(c)
        for s in sj: rows.append(dict(support=hit[s], supported=s in colsplits, correct=s in true))
    V = pd.DataFrame(rows); hi = V[V.support >= 95]; print("H=0.002, 1,000 sites: splits with bootstrap >= 95:", len(hi))
    for flag, g in hi.groupby("supported"): print(f"  split supported by at least one alignment column = {flag}: n = {len(g)}, correct {100 * g.correct.mean():.1f}%")
    print(f"  share of all NJ splits with no supporting column: {100 * (~V.supported).mean():.1f}%"); V.to_csv("calibration_check_raw.csv", index=False)
    A = np.load("support_H3N2.npy"); Bq = np.load("support_Ebola.npy"); acc = pd.read_csv("table_accuracy_sim.csv"); cal = pd.read_csv("table_calibration_sim.csv"); cal["height"] = [0.002] * 4 + [0.02] * 4
    fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7))
    for (nm, M, c) in [("H3N2 (UFBoot)", A, "#2b6cb0"), ("Ebola (SH-aLRT)", Bq, "#dd6b20")]:
        x = np.sort(M[1][~np.isnan(M[1])]); ax[0].plot(x, np.arange(1, len(x) + 1) / len(x) * 100, color=c, lw=2, label=nm); x2 = np.sort(M[0]); ax[0].plot(x2, np.arange(1, len(x2) + 1) / len(x2) * 100, color=c, lw=1.4, ls="--", label=nm.split()[0] + " (own NJ bootstrap)")
    ax[0].set(xlabel="Support of a split (%)", ylabel="Splits at or below this support (%)", title="Support in the real trees"); ax[0].legend(frameon=False, fontsize=7.6)
    for H, c in [(0.002, "#dd6b20"), (0.02, "#2b6cb0")]:
        a = acc[acc.tree_height == H]; ax[1].plot(a.informative_sites, a.nj_rf, "o-", color=c, label=f"NJ, height {H}"); ax[1].plot(a.informative_sites, a.fasttree_rf, "s--", color=c, label=f"FastTree, height {H}")
    ax[1].set_xscale("log"); ax[1].set(xlabel="Informative sites in the alignment", ylabel="Robinson-Foulds distance to the true tree", title="Accuracy rises with information", ylim=(0, 1)); ax[1].legend(frameon=False, fontsize=7.6)
    w = 0.38
    for k, (H, c) in enumerate([(0.002, "#dd6b20"), (0.02, "#2b6cb0")]): g = cal[cal.height == H]; ax[2].bar(np.arange(len(g)) + (k - 0.5) * w, g.correct_pct, w, color=c, label=f"height {H}")
    ax[2].plot([-0.5, 3.5], [25, 25], alpha=0); ax[2].plot([0, 1, 2, 3], [25, 60, 82, 97], "k:", lw=1, label="nominal support (bin centre)"); ax[2].set_xticks(range(4)); ax[2].set_xticklabels(["<50", "50-69", "70-94", "95-100"]); ax[2].set(xlabel="NJ bootstrap support (%)", ylabel="Splits that are correct (%)", title="Calibration of bootstrap support", ylim=(0, 100)); ax[2].legend(frameon=False, fontsize=7.6)
    fig.tight_layout(); fig.savefig("figs5/fig3_support_accuracy.png", dpi=200); plt.close(fig); print("figure written")
if __name__ == "__main__":
    {"A": part_A, "B": part_B, "C": part_C}[sys.argv[1]]()
