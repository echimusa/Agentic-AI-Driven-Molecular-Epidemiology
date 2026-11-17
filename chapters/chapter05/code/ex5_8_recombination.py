"""
Chapter 5, Hands-on 7 - Recombination, branch lengths and the tree (a simulation).
Alignments of 30 taxa and 10,000 sites are simulated along a coalescent tree (HKY, kappa 4, expected root height 0.004 substitutions per site) with R recombination events. Each event copies a block of 400 sites into every taxon below a randomly chosen branch from a donor
sequence that has diverged by a further 0.05 substitutions per site, which brings about 20 SNPs into one gene of one clade. Three trees are built with FastTree (GTR): from the clonal alignment (no events), from the recombinant alignment, and from the recombinant alignment after
masking blocks found by a screen written here (a simplified, Gubbins-like idea; it is not Gubbins). The screen reconstructs the substitutions on every branch with Fitch parsimony and flags windows of 500 sites in which a branch has an excess of substitutions (Poisson tail below 1e-4).
Compared: tree length against the true clonal tree, Robinson-Foulds distance to the true tree, the share of imported bases masked (sensitivity) and the share of clonal bases masked (false positives).
Output: table_recombination.csv, Figure 6. Needs FastTree on the path (about 6 minutes).
"""
import numpy as np, pandas as pd, subprocess, io, os, sys, tempfile
from scipy import stats
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo_utils import *
from phylo_utils2 import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
KAPPA, PIV = 4.0, np.array([0.3, 0.2, 0.2, 0.3]); NT, LEN, BLOCK, DONOR, WIN = 30, 10000, 400, 0.05, 500
def hky_P(t):
    Rm = np.ones((4, 4)); Rm[0, 2] = Rm[2, 0] = Rm[1, 3] = Rm[3, 1] = KAPPA; Q = Rm * PIV[None, :]; np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); Q = Q / -(PIV * np.diag(Q)).sum()
    s = np.sqrt(PIV); S = (s[:, None] * Q) / s[None, :]; S = (S + S.T) / 2; w, U = np.linalg.eigh(S); return (U * np.exp(w * t)) @ U.T * (1 / s)[:, None] * s[None, :]
def coalescent(n, height, rng):
    nodes = list(range(n)); t = 0.0; kids = {}; age = {i: 0.0 for i in nodes}; nxt = n; k = n
    while len(nodes) > 1:
        t += rng.exponential(1.0 / (k * (k - 1) / 2)); a, b = rng.choice(len(nodes), 2, replace=False); A, B = nodes[a], nodes[b]; kids[nxt] = (A, B); age[nxt] = t; nodes = [x for x in nodes if x not in (A, B)] + [nxt]; nxt += 1; k -= 1
    sc = height / age[nxt - 1]; return kids, {v: a * sc for v, a in age.items()}, nxt - 1
def evolve(kids, age, root, L, rng, seq0=None):
    seq = {root: rng.choice(4, L, p=PIV) if seq0 is None else seq0}; stack = [root]
    while stack:
        v = stack.pop()
        for c in kids.get(v, ()):
            P = hky_P((age[v] - age[c])); cum = P.cumsum(1); u = rng.random(L); seq[c] = (u[:, None] > cum[seq[v]]).sum(1).clip(0, 3).astype(np.int8); stack.append(c)
    return seq
def descendants(kids, v): return [v] if v not in kids else [x for c in kids[v] for x in descendants(kids, c)]
def fasttree(E, names):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        for i, nm in enumerate(names): f.write(f">{nm}\n{''.join('AGCTN'[x] if x >= 0 else 'N' for x in E[i])}\n")
        fn = f.name
    out = subprocess.run(["FastTree", "-nt", "-gtr", "-quiet", "-nosupport", fn], capture_output=True, text=True).stdout; os.unlink(fn); return Phylo.read(io.StringIO(out), "newick")
def tree_length(tr): return sum(c.branch_length or 0 for c in tr.find_clades() if c is not tr.root)
def screen(tr, names, E, alpha=1e-4):
    """Flag, for every branch of the inferred tree, windows holding more substitutions than a uniform spread would give; returns a boolean mask (taxa x sites) over the taxa below the branch."""
    tb = binarize(Phylo.read(io.StringIO(tr.format("newick")), "newick")); _, ch = fitch(tb, names, E); idx = {n: i for i, n in enumerate(names)}; mask = np.zeros(E.shape, bool); starts = np.arange(0, E.shape[1] - WIN + 1, WIN // 5)
    for cl, pos in ch.items():
        if len(pos) < 4: continue
        counts = np.array([((pos >= s) & (pos < s + WIN)).sum() for s in starts]); lam = len(pos) * WIN / E.shape[1]; p = stats.poisson.sf(counts - 1, lam)
        hit = (p < alpha) & (counts >= 4)
        if hit.any():
            cols = np.zeros(E.shape[1], bool)
            for s in starts[hit]: cols[s:s + WIN] = True
            rows = [idx[t.name] for t in cl.get_terminals()]; mask[np.ix_(rows, np.where(cols)[0])] = True
    return mask
def one(R, seed):
    rng = np.random.default_rng(seed); kids, age, root = coalescent(NT, 0.004, rng); seq = evolve(kids, age, root, LEN, rng); E = np.array([seq[i] for i in range(NT)]).astype(np.int8); names = [f"t{i}" for i in range(NT)]
    true = Phylo.read(io.StringIO(truetree(kids, age, root)), "newick"); sp_true = splits(true, names); tl_true = sum(age[v] - age[c] for v, cs in kids.items() for c in cs); clonal = E.copy(); truth = np.zeros(E.shape, bool)
    inner = [v for v in kids if v != root] + list(range(NT))
    for _ in range(R):
        v = inner[rng.integers(len(inner))]; taxa = [x for x in descendants(kids, v) if x < NT]; a = int(rng.integers(0, LEN - BLOCK)); donor_kids, donor_age, droot = {v: kids.get(v, ())}, age, v
        base = seq[v][a:a + BLOCK]; P = hky_P(DONOR); cum = P.cumsum(1); donor = (rng.random(BLOCK)[:, None] > cum[base]).sum(1).clip(0, 3).astype(np.int8)
        for tx in taxa: E[tx, a:a + BLOCK] = donor; truth[tx, a:a + BLOCK] = True
    out = {}
    for label, M in [("clonal", clonal), ("recombinant", E)]:
        tr = fasttree(M, names); out[label] = (tree_length(tr) / tl_true, rf(splits(tr, names, min_len=0.5 / LEN), sp_true)); trees = tr
    rec_tree = fasttree(E, names); mask = screen(rec_tree, names, E); Em = E.copy(); Em[mask] = -1; mt = fasttree(Em, names); out["masked"] = (tree_length(mt) / tl_true, rf(splits(mt, names, min_len=0.5 / LEN), sp_true))
    sens = (mask & truth).sum() / max(truth.sum(), 1) if R else np.nan; fpr = (mask & ~truth).sum() / (~truth).sum()
    return out, sens, fpr, int(truth.sum() / BLOCK)
def truetree(kids, age, root):
    def f(v): return f"t{v}" if v not in kids else "(" + ",".join(f"{f(c)}:{age[v] - age[c]:.8f}" for c in kids[v]) + ")"
    return f(root) + ";"
def make_figure(D, M):
    fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.3)); cols = {"clonal": "#2f855a", "recombinant": "#c53030", "masked": "#2b6cb0"}
    for lab, c in cols.items():
        a = D[D.tree == lab].groupby("events"); ax[0].errorbar(a.tree_length_ratio.mean().index, a.tree_length_ratio.mean(), a.tree_length_ratio.std(), color=c, marker="o", lw=2, capsize=3, label=lab); ax[1].errorbar(a.rf_to_truth.mean().index, a.rf_to_truth.mean(), a.rf_to_truth.std(), color=c, marker="o", lw=2, capsize=3, label=lab)
    ax[0].axhline(1, color="black", ls=":", lw=1); ax[0].set(xlabel="Recombination events per alignment", ylabel="Inferred tree length / true clonal tree length", title="Recombination inflates the tree"); ax[0].legend(frameon=False)
    ax[1].set(xlabel="Recombination events per alignment", ylabel="Robinson-Foulds distance to the true tree", title="but hardly moves the topology"); w = 0.36; ev = [3, 10]; x = np.arange(2); mm = M.set_index("events")
    b1 = ax[2].bar(x - w / 2, [mm.sensitivity[e] * 100 for e in ev], w, color="#2b6cb0", label="imported bases masked (sensitivity)"); b2 = ax[2].bar(x + w / 2, [mm.false_positive[e] * 100 for e in ev], w, color="#dd6b20", label="other bases masked")
    for bars in (b1, b2):
        for r in bars: ax[2].text(r.get_x() + r.get_width() / 2, r.get_height() + 1.0, f"{r.get_height():.0f}%", ha="center", fontsize=10)
    ax[2].set(xticks=x, xticklabels=[f"{e} events" for e in ev], ylabel="% of bases", ylim=(0, 112), yticks=[0, 20, 40, 60, 80, 100], title="What the screen masks"); ax[2].legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    plt.tight_layout(); os.makedirs("figs5", exist_ok=True); plt.savefig("figs5/fig6_recombination.png", dpi=170)

if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "figure":          # redraw Figure 6 from the saved replicates
        D = pd.read_csv("table_recombination_raw.csv"); M = pd.read_csv("table_recombination.csv"); make_figure(D, M[M.tree == "masked"]); sys.exit()
    rows = []; reps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    for R in (0, 3, 10):
        for r in range(reps):
            out, sens, fpr, cells = one(R, 1000 * R + r)
            for lab, (tl, d) in out.items(): rows.append(dict(events=R, rep=r, tree=lab, tree_length_ratio=tl, rf_to_truth=d, sensitivity=sens, false_positive=fpr))
        print("events", R, "done", flush=True)
    D = pd.DataFrame(rows); D.round(5).to_csv("table_recombination_raw.csv", index=False)
    T = D.groupby(["events", "tree"]).agg(tree_length_ratio=("tree_length_ratio", "mean"), rf_to_truth=("rf_to_truth", "mean")).reset_index(); M = D[D.tree == "masked"].groupby("events").agg(sensitivity=("sensitivity", "mean"), false_positive=("false_positive", "mean")).reset_index()
    T = T.merge(M, on="events", how="left"); T.loc[T.tree != "masked", ["sensitivity", "false_positive"]] = np.nan; T.round(4).to_csv("table_recombination.csv", index=False); print(T.round(4).to_string(index=False))
    make_figure(D, M)
