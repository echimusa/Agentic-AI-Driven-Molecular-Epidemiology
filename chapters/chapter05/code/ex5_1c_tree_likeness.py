"""
Chapter 5, Hands-on 1, parts C and D - How tree-like are the distances, and how do the trees score?
Part C: the four-point condition (additive distances: neighbour-joining is exact) and the three-point condition (ultrametric distances: UPGMA is exact) are tested on random quartets and triples of the real alignments (195 Ebola virus genomes, 60 H3N2 HA sequences), and the observed JC69
distances are compared with the path lengths of the maximum-likelihood tree.
Part D: the maximum-likelihood, neighbour-joining and UPGMA trees are scored by Fitch parsimony (written here) and by log-likelihood under the model of the ML tree (IQ-TREE, topology fixed, branch lengths optimised).
Output: table_tree_likeness.csv, table_tree_scores.csv. Needs h3n2.treefile, ebola.treefile and the alignments; IQ-TREE 2 on the path for part D.
"""
import numpy as np, pandas as pd, os, re
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from phylo_utils import *
from phylo_utils2 import *
sets = {"H3N2 HA (60 sequences)": ("h3n2_ha_aln.fasta", "h3n2.treefile", "TPM2+F+G4"), "Ebola virus (195 genomes)": ("ebola_aln.fasta", "ebola.treefile", "GTR+F+G4")}
rng = np.random.default_rng(5); likeness, scores = [], []
for name, (aln, tf, model) in sets.items():
    names, X = read_aln(aln); E = encode(X); L = X.shape[1]; n = len(names); p, P, Q = pairwise(E); D = p * L          # distances in substitutions (p-distance times the alignment length)
    # ---- part C: four-point and three-point conditions on random quartets and triples
    m = 200000; q = np.array([rng.choice(n, 4, replace=False) for _ in range(m)]); i, j, k, l = q.T
    S = np.stack([D[i, j] + D[k, l], D[i, k] + D[j, l], D[i, l] + D[j, k]], 1); S.sort(1); ex4 = S[:, 2] - S[:, 1]               # largest minus middle sum: zero for tree-like distances
    t3 = np.array([rng.choice(n, 3, replace=False) for _ in range(m)]); a, b, c = t3.T; T = np.stack([D[a, b], D[a, c], D[b, c]], 1); T.sort(1); ex3 = T[:, 2] - T[:, 1]    # largest minus middle distance: zero for ultrametric distances
    ml = Phylo.read(tf, "newick"); names_set = set(names); dm = {}; leaves = ml.get_terminals(); tips = [t.name for t in leaves]; idx = [names.index(t) for t in tips]
    pat = np.zeros((n, n))
    for ii, ti in enumerate(leaves):
        for jj in range(ii + 1, len(leaves)): pat[idx[ii], idx[jj]] = pat[idx[jj], idx[ii]] = ml.distance(ti, leaves[jj])
    iu = np.triu_indices(n, 1); obs = jc(p)[iu]; tre = pat[iu]; r = np.corrcoef(obs, tre)[0, 1]; sl = np.polyfit(tre, obs, 1)[0]; mae = np.mean(np.abs(obs - tre)) * L
    likeness.append(dict(dataset=name, taxa=n, quartets=m, four_point_zero_pct=100 * (ex4 < 0.5).mean(), four_point_median=np.median(ex4), four_point_p95=np.percentile(ex4, 95), triples=m, three_point_zero_pct=100 * (ex3 < 0.5).mean(), three_point_median=np.median(ex3), three_point_p95=np.percentile(ex3, 95), corr_obs_vs_tree=r, slope=sl, mean_abs_residual_subs=mae))
    # ---- part D: scores of ML, NJ and UPGMA trees
    dJ = jc(p); trees = {"ML (IQ-TREE)": ml, "NJ (JC69)": DistanceTreeConstructor().nj(to_dm(names, np.nan_to_num(dJ))), "UPGMA (JC69)": DistanceTreeConstructor().upgma(to_dm(names, np.nan_to_num(dJ)))}; sp_ml = splits(ml, names, min_len=0.5 / L)
    for label, tr in trees.items():
        sp = splits(tr, names, min_len=0.5 / L); tb = binarize(Phylo.read(__import__("io").StringIO(tr.format("newick")), "newick")); cost, _ = fitch(tb, names, E)
        pre = f"score_{name.split()[0]}_{label.split()[0]}"; ll = iqtree_loglik(aln, tr, model, pre); scores.append(dict(dataset=name, tree=label, rf_to_ml=rf(sp, sp_ml), parsimony_score=cost, loglik=ll))
S_ = pd.DataFrame(scores); S_["parsimony_excess"] = S_.parsimony_score - S_.groupby("dataset").parsimony_score.transform("min"); S_["loglik_deficit"] = S_.groupby("dataset").loglik.transform("max") - S_.loglik
pd.DataFrame(likeness).round(3).to_csv("table_tree_likeness.csv", index=False); S_.round(2).to_csv("table_tree_scores.csv", index=False)
print(pd.DataFrame(likeness).round(3).T.to_string()); print(S_.round(2).to_string(index=False))
