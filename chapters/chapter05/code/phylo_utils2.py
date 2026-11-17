"""Helpers added for the extended Chapter 5: Fitch parsimony on a rooted binary tree (score and the branches on which each site changes), binarisation of a tree and IQ-TREE evaluation of a fixed tree."""
import io, os, re, subprocess, tempfile
import numpy as np
from Bio import Phylo

BIT = np.array([1, 2, 4, 8], dtype=np.uint8)          # A G C T as bits; the encoding of phylo_utils is A=0 G=1 C=2 T=3, gap or ambiguity -1

def to_sets(E):
    """Encoded sites (taxa x sites, -1 for gap or ambiguity) as 4-bit state sets; a gap or ambiguity is the full set (compatible with everything)."""
    S = np.full(E.shape, 15, dtype=np.uint8); ok = E >= 0; S[ok] = BIT[E[ok]]; return S

def binarize(tree):
    """Midpoint-root the tree and resolve every polytomy into zero-length bifurcations, so that the Fitch algorithm applies."""
    tree.root_at_midpoint()
    for cl in list(tree.get_nonterminals()):
        while len(cl.clades) > 2:
            a, b = cl.clades[-2], cl.clades[-1]; new = type(cl)(branch_length=0.0, clades=[a, b]); cl.clades = cl.clades[:-2] + [new]
    return tree

def fitch(tree, names, E):
    """Fitch (1971) parsimony on a rooted binary tree. Returns the total score and, for every non-root clade, the array of sites on which the base changes on the branch above it
    (traceback with an arbitrary choice at the root, which gives a minimum-length assignment)."""
    S = to_sets(E); idx = {n: i for i, n in enumerate(names)}; n_sites = E.shape[1]; sets = {}; cost = 0
    def up(cl):
        nonlocal cost
        if cl.is_terminal(): sets[cl] = S[idx[cl.name]].copy(); return
        for c in cl.clades: up(c)
        a, b = (sets[c] for c in cl.clades); inter = a & b; empty = inter == 0; sets[cl] = np.where(empty, a | b, inter); cost += int(empty.sum())
    up(tree.root)
    low = lambda m: np.where(m & 1, 1, np.where(m & 2, 2, np.where(m & 4, 4, 8))).astype(np.uint8)       # the lowest set bit of each mask
    assign = {tree.root: low(sets[tree.root])}; changes = {}
    def down(cl):
        for c in cl.clades:
            par = assign[cl]; keep = (sets[c] & par) != 0; assign[c] = np.where(keep, par, low(sets[c])); changes[c] = np.where(assign[c] != par)[0]; down(c)
    down(tree.root); return cost, changes

def iqtree_loglik(aln, tree, model, prefix):
    """Log-likelihood of a fixed topology (branch lengths optimised) under a model, from IQ-TREE; the tree is written to a temporary file."""
    with open(prefix + ".in.nwk", "w") as f: Phylo.write(tree, f, "newick")
    subprocess.run(["iqtree2", "-s", aln, "-te", prefix + ".in.nwk", "-m", model, "-nt", "1", "--prefix", prefix, "-redo", "-quiet"], capture_output=True)
    txt = open(prefix + ".iqtree").read(); return float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", txt).group(1))
