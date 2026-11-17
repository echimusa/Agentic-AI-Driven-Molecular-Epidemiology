"""
Chapter 5, Hands-on 3, part C - The Ebola tree with the full model set and the ultrafast bootstrap.
The Ebola tree of the other hands-on exercises came from a restricted model set (JC, K80, HKY and GTR, with equal or gamma rates) and the SH-aLRT. Here the whole of ModelFinder (all models, including invariant sites) and 1,000 ultrafast bootstrap replicates are used:
    iqtree2 -s ebola_aln.fasta -m MFP -B 1000 -alrt 1000 -nt 1 -seed 1 --prefix ebola_full
This script reads the outputs (ebola_full.iqtree, ebola_full.treefile; shipped, or regenerated with the command above, about 5 minutes) and compares them with the restricted run (ebola.treefile).
Output: table_ebola_full_support.csv. Compared: the model chosen and the BIC gaps, the Robinson-Foulds distance between the two trees, and the support of the resolved branches (at least one substitution) and of the unresolved ones.
"""
import re, numpy as np, pandas as pd
from Bio import Phylo
from phylo_utils import *
names, X = read_aln("ebola_aln.fasta"); L = X.shape[1]; iq = open("ebola_full.iqtree").read()
tab = re.findall(r"^(\S+)\s+(-?[\d.]+)\s+([\d.]+)\s+[+-]\s+[\d.e+-]+\s+([\d.]+)\s+[+-]\s+[\d.e+-]+\s+([\d.]+)\s+[+-]\s+[\d.e+-]+\s*$", iq.split("List of models sorted by BIC scores:")[1].split("AIC, w-AIC")[0], flags=re.M)
bic = {m: float(b) for m, ll, a, c, b in tab}; best = min(bic, key=bic.get); gaps = {m: bic[m] - bic[best] for m in ("GTR+F+G4", "HKY+F+G4", "HKY+F+I", "GTR+F+I", "GTR+F+I+G4", "HKY+F+I+G4") if m in bic}
ref = Phylo.read("ebola.treefile", "newick"); full = Phylo.read("ebola_full.treefile", "newick"); s_ref = splits(ref, names, min_len=0.5 / L); s_full = splits(full, names, min_len=0.5 / L); idx = {n: i for i, n in enumerate(names)}; allset = frozenset(range(len(names)))
rows = []
for cl in full.get_nonterminals():
    if cl is full.root: continue
    lab = cl.confidence if cl.confidence is not None else cl.name; m = re.match(r"([\d.]+)/([\d.]+)", str(lab)) if lab is not None else None
    sh, ufb = (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan); rows.append(dict(length=cl.branch_length or 0, sh=sh, ufb=ufb))
B = pd.DataFrame(rows); res = B[B.length >= 0.5 / L]; unr = B[B.length < 0.5 / L]
out = dict(best_model_bic=best, loglik=float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", iq).group(1)), bic_gap_GTR_F_G4=gaps.get("GTR+F+G4"), bic_gap_HKY_F_G4=gaps.get("HKY+F+G4"), bic_gap_GTR_F_I_G4=gaps.get("GTR+F+I+G4"), rf_restricted_vs_full=rf(s_ref, s_full), splits_restricted=len(s_ref), splits_full=len(s_full), shared=len(s_ref & s_full),
           internal_branches=len(B), resolved_branches=len(res), ufboot_ge95_resolved_pct=100 * (res.ufb >= 95).mean(), ufboot_ge80_resolved_pct=100 * (res.ufb >= 80).mean(), shalrt_ge80_resolved_pct=100 * (res.sh >= 80).mean(), both_ufboot95_shalrt80_resolved_pct=100 * ((res.ufb >= 95) & (res.sh >= 80)).mean(),
           unresolved_branches=len(unr), ufboot_ge95_unresolved_pct=100 * (unr.ufb >= 95).mean(), ufboot_lt50_unresolved_pct=100 * (unr.ufb < 50).mean(), corr_ufboot_shalrt_resolved=float(np.corrcoef(res.ufb, res.sh)[0, 1]))
pd.DataFrame([out]).round(3).to_csv("table_ebola_full_support.csv", index=False); print(pd.DataFrame([out]).round(3).T.to_string()); print(sorted(bic.items(), key=lambda t: t[1])[:6])
