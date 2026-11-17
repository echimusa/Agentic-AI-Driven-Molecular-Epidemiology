"""
Chapter 5, Hands-on 6 - Ascertainment bias: a tree from variable sites only.
A core-genome SNP alignment (the output of snippy-core in Chapter 4) holds only the variable sites. A model that assumes every site could have been constant is then wrong: the likelihood must be conditioned on the fact that each site is variable, which IQ-TREE does with the +ASC option.
For each dataset (195 Ebola virus genomes, 60 H3N2 HA sequences) a tree is built from (i) the full alignment, (ii) the variable sites only without a correction, and (iii) the variable sites with +ASC, and (iv) the variable sites with the counts of the constant sites supplied (-fconst), under the same substitution model with gamma rates.
Compared: tree length per site (the scale a molecular clock would read), the substitutions that tree length claims for the whole genome, and the share of the resolved splits of the full-alignment tree that each tree recovers.
Output: table_ascertainment.csv. Needs IQ-TREE 2 on the path.
"""
import numpy as np, pandas as pd, subprocess, re, os
from Bio import Phylo
from phylo_utils import *
sets = {"H3N2 HA (60 sequences)": ("h3n2_ha_aln.fasta", "HKY+F+G4"), "Ebola virus (195 genomes)": ("ebola_aln.fasta", "HKY+F+G4")}
def run(aln, model, prefix, extra=()):
    subprocess.run(["iqtree2", "-s", aln, "-m", model, "-nt", "1", "--prefix", prefix, "-seed", "1", "-redo", "-quiet", *extra], capture_output=True)
    tr = Phylo.read(prefix + ".treefile", "newick"); txt = open(prefix + ".iqtree").read(); tl = float(re.search(r"Total tree length \(sum of branch lengths\): ([\d.]+)", txt).group(1)); ll = float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", txt).group(1)); return tr, tl, ll
rows = []
for name, (aln, model) in sets.items():
    names, X = read_aln(aln); L = X.shape[1]; ok = np.isin(X, list("ACGT")); var = np.array([len(set(X[:, s][ok[:, s]])) >= 2 for s in range(L)]); nv = int(var.sum())
    const = [int(sum(1 for s in np.where(~var)[0] if set(X[:, s][ok[:, s]]) == {b})) for b in "ACGT"]          # constant sites of the full alignment by base (A, C, G, T), the counts IQ-TREE's -fconst needs
    with open("snp_" + aln, "w") as f:
        for i, nme in enumerate(names): f.write(f">{nme}\n{''.join(X[i, var])}\n")
    tag = name.split()[0]; full, tl_f, ll_f = run(aln, model, f"asc_{tag}_full"); sp_full = splits(full, names, min_len=0.5 / L)
    variants = [("full alignment", model, aln, (), L), ("variable sites, no correction", model, "snp_" + aln, (), nv), ("variable sites, +ASC", model.replace("+F", "+F+ASC"), "snp_" + aln, (), nv), ("variable sites, constant-site counts (-fconst)", model, "snp_" + aln, ("-fconst", ",".join(map(str, const))), nv)]
    for label, mod, aln_used, extra, n_sites in variants:
        tr, tl, ll = (full, tl_f, ll_f) if label == "full alignment" else run(aln_used, mod, f"asc_{tag}_{label.split()[-1].strip('(-)').replace('+', '')}", extra); sp = splits(tr, names, min_len=1e-6)
        rows.append(dict(dataset=name, alignment=label, sites=n_sites, tree_length_per_site=tl, tree_length_ratio_to_full=tl / tl_f, substitutions_claimed_for_the_genome=tl * L, resolved_splits_of_full_tree_found=len(sp & sp_full) / len(sp_full)))
R = pd.DataFrame(rows); R.round(4).to_csv("table_ascertainment.csv", index=False); print(R.round(4).to_string(index=False))
