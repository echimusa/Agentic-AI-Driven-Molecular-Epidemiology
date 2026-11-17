"""Shared helpers for Chapter 5: reading alignments, distances (p, JC69, K80), tree splits and Robinson-Foulds distance, and a small tree drawing function."""
import io, re
import numpy as np
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix

def read_aln(path):
    a = AlignIO.read(path, "fasta"); names = [r.id for r in a]; X = np.array([list(str(r.seq).upper()) for r in a]); return names, X

PUR = set("AG"); PYR = set("CT")
def encode(X):
    m = {"A": 0, "G": 1, "C": 2, "T": 3}; E = np.full(X.shape, -1, dtype=np.int8)
    for k, v in m.items(): E[X == k] = v
    return E

def pairwise(E):
    """Proportion of differing sites (p), of transitions (P) and of transversions (Q) for every pair, ignoring sites with a gap or ambiguity in either sequence (pairwise deletion)."""
    n = E.shape[0]; p = np.zeros((n, n)); P = np.zeros((n, n)); Q = np.zeros((n, n)); ok = E >= 0
    for i in range(n):
        for j in range(i + 1, n):
            v = ok[i] & ok[j]; a, b = E[i][v], E[j][v]; L = max(v.sum(), 1); diff = a != b; ts = diff & ((a // 2) == (b // 2))      # A,G are 0,1 and C,T are 2,3, so a shared half of the code is a transition
            p[i, j] = p[j, i] = diff.sum() / L; P[i, j] = P[j, i] = ts.sum() / L; Q[i, j] = Q[j, i] = (diff.sum() - ts.sum()) / L
    return p, P, Q

def jc(p): return -0.75 * np.log(1 - 4.0 / 3.0 * p)
def k80(P, Q): return -0.5 * np.log(1 - 2 * P - Q) - 0.25 * np.log(1 - 2 * Q)

def to_dm(names, D): return DistanceMatrix(names, [list(D[i, :i + 1]) for i in range(len(names))])

def splits(tree, taxa, min_len=0.0):
    """Non-trivial bipartitions of an unrooted tree as frozensets of the side NOT containing taxon 0; branches shorter than min_len are treated as absent (collapsed)."""
    idx = {t: i for i, t in enumerate(taxa)}; out = set(); full = frozenset(range(len(taxa)))
    for cl in tree.get_nonterminals():
        if cl is tree.root: continue
        if cl.branch_length is not None and cl.branch_length < min_len: continue
        s = frozenset(idx[t.name] for t in cl.get_terminals())
        if 1 < len(s) < len(taxa) - 1:
            if 0 in s: s = full - s
            if 1 < len(s) < len(taxa) - 1: out.add(s)
    return out

def rf(a, b):
    """Normalised Robinson-Foulds distance between two split sets: the share of splits found in only one tree."""
    return len(a ^ b) / max(len(a) + len(b), 1)

def draw_tree(ax, tree, colors=None, lw=0.6, label_tips=False):
    """Rectangular phylogram drawn with matplotlib; colors maps tip name to colour."""
    pos = {}; y = [0]
    def ypos(cl):
        if cl.is_terminal(): pos[cl] = [None, y[0]]; y[0] += 1; return pos[cl][1]
        ys = [ypos(c) for c in cl.clades]; pos[cl] = [None, (min(ys) + max(ys)) / 2]; return pos[cl][1]
    ypos(tree.root)
    def xpos(cl, x0):
        x = x0 + (cl.branch_length or 0.0); pos[cl][0] = x
        for c in cl.clades: xpos(c, x)
    xpos(tree.root, 0.0)
    for cl in pos:
        x, yv = pos[cl]
        for c in cl.clades: ax.plot([x, x], [yv, pos[c][1]], color="#4a5568", lw=lw); ax.plot([x, pos[c][0]], [pos[c][1]] * 2, color=(colors or {}).get(c.name, "#4a5568") if c.is_terminal() else "#4a5568", lw=lw)
    for cl in tree.get_terminals():
        if colors and cl.name in colors: ax.plot(pos[cl][0], pos[cl][1], "o", ms=3, color=colors[cl.name])
    return pos
