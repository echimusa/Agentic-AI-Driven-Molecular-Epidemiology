"""
Chapter 9, Hands-on 3 - How well do resistance determinants predict the susceptibility test result?
Uses every Salmonella isolate in NCBI Pathogen Detection that has both an AMRFinderPlus genotype and a broth-microdilution result.
Rule: a determinant for the drug's class predicts resistance (R); no determinant predicts susceptible (S). Intermediate results are left out.
Output: table_concordance.csv, Figure 6.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from amr9 import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})

m = add_sector(load_metadata()); a = m[m.AST_phenotypes.notna()].copy(); cat = load_catalog()
a["genes"] = a.AMR_genotypes_core.map(parse_genotype); a["ast"] = a.AST_phenotypes.map(parse_ast)
print(f"{len(a):,} isolates with a susceptibility result; sector: {a.sector.value_counts().to_dict()}; years {int(a.year.min())}-{int(a.year.max())}; countries: {a.country.value_counts().head(3).to_dict()}")
rows, tests = [], []
for drug, rule in RULES.items():
    obs = a.ast.map(lambda d: d.get(drug)); ok = obs.isin(["S", "R"])
    if ok.sum() < 30: continue
    pred = a.genes[ok].map(lambda g: rule(g, cat)); y = (obs[ok] == "R")
    tp, fn, fp, tn = int((pred & y).sum()), int((~pred & y).sum()), int((pred & ~y).sum()), int((~pred & ~y).sum())
    sens, spec, ppv, npv = tp / max(tp + fn, 1), tn / max(tn + fp, 1), tp / max(tp + fp, 1), tn / max(tn + fn, 1)
    rows.append(dict(drug=drug, n=int(ok.sum()), n_R=tp + fn, tp=tp, fn=fn, fp=fp, tn=tn, agreement=(tp + tn) / ok.sum(), sens=sens, sens_lo=wilson(tp, tp + fn)[0], sens_hi=wilson(tp, tp + fn)[1],
                     spec=spec, spec_lo=wilson(tn, tn + fp)[0], spec_hi=wilson(tn, tn + fp)[1], ppv=ppv, npv=npv))
    tests.append(((pred == y).sum(), len(y)))
t = pd.DataFrame(rows); t.to_csv("table_concordance.csv", index=False)
agree = sum(x for x, _ in tests); total = sum(n for _, n in tests); print(f"overall: {agree:,} of {total:,} drug-isolate tests agree ({agree / total:.1%})")
print(t[["drug", "n", "n_R", "agreement", "sens", "spec", "ppv", "fn", "fp"]].round(3).sort_values("agreement").to_string(index=False))

# ciprofloxacin counted as non-susceptible when I or R (the reduced-susceptibility result that qnr genes produce)
obs = a.ast.map(lambda d: d.get("ciprofloxacin")); ok = obs.isin(["S", "I", "R"]); pred = a.genes[ok].map(lambda g: RULES["ciprofloxacin"](g, cat)); y = obs[ok].isin(["I", "R"])
print(f"ciprofloxacin I or R: n={ok.sum()}, sensitivity {(pred & y).sum() / y.sum():.3f}, specificity {(~pred & ~y).sum() / (~y).sum():.3f}; I results: {(obs[ok] == 'I').sum()}")

# which determinants sit behind the errors?
for drug, kind in [("streptomycin", "fn"), ("gentamicin", "fn"), ("ampicillin", "fp"), ("sulfamethoxazole", "fp"), ("azithromycin", "fn")]:
    if drug not in RULES: continue
    o = a.ast.map(lambda d: d.get(drug)); ok = o.isin(["S", "R"]); sub = a[ok]; pred = sub.genes.map(lambda g: RULES[drug](g, cat)); y = o[ok] == "R"
    bad = sub[(~pred & y)] if kind == "fn" else sub[(pred & ~y)]
    tops = bad.genes.explode().value_counts().head(4).to_dict() if len(bad) else {}
    print(f"{drug}: {len(bad)} {'resistant isolates without a determinant' if kind == 'fn' else 'susceptible isolates with a determinant'}; most common elements among them: {tops}")

# ---------------------------------------------------------------- Figure 5
d = t.sort_values("sens").reset_index(drop=True); y = np.arange(len(d))
fig, ax = plt.subplots(1, 2, figsize=(11.6, 5.4), sharey=True)
for a_, col, lo, hi, ttl in [(ax[0], "sens", "sens_lo", "sens_hi", "Sensitivity: resistant isolates with a determinant"), (ax[1], "spec", "spec_lo", "spec_hi", "Specificity: susceptible isolates without one")]:
    few = d.n_R < 20
    a_.errorbar(d[col], y, xerr=[np.clip(d[col] - d[lo], 0, None), np.clip(d[hi] - d[col], 0, None)], fmt="none", ecolor="#a0aec0", capsize=2.5, zorder=1)
    a_.scatter(d[col][~few], y[~few], color="#2b6cb0", zorder=3, s=36, label="20 or more resistant isolates"); a_.scatter(d[col][few], y[few], facecolor="white", edgecolor="#2b6cb0", zorder=3, s=36, label="fewer than 20 resistant isolates")
    a_.set(title=ttl, xlim=(0.0, 1.02) if col == "sens" else (0.9, 1.005), xlabel="Proportion"); a_.axvline(0.95, color="#c53030", ls=":", lw=1)
ax[0].set_yticks(y); ax[0].set_yticklabels([f"{r.drug}  (R = {r.n_R:,} of {r.n:,})" for r in d.itertuples()], fontsize=8.8)
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=8.5, loc="lower center", ncol=2); fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig("figs9/fig6_genotype_phenotype.png", dpi=200); plt.close(fig)
