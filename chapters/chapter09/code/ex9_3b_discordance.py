"""
Chapter 9, Hands-on 3 (extension) - why do the genotype and the susceptibility result disagree?
For each drug of ex9_3_genotype_phenotype.py, every resistant isolate that the rule calls susceptible (a false negative, 'very major error') and every susceptible isolate that the rule calls
resistant (a false positive, 'major error') is examined for what it carries: no core determinant at all, determinants of other classes only, and which serovar and period it comes from.
Azithromycin is followed in detail because its sensitivity is 33%.  Input: sal_ast_subset.csv.gz and the AMRFinderPlus catalogue (amrf_db_catalog/).  Output: table_discordance.csv, table_azithromycin_period.csv.
"""
import numpy as np, pandas as pd
from amr9 import *
a = pd.read_csv("sal_ast_subset.csv.gz", low_memory=False); cat = load_catalog(); a["genes"] = a.AMR_genotypes_core.map(parse_genotype); a["ast"] = a.AST_phenotypes.map(parse_ast)
a["year"] = pd.to_datetime(a.collection_date, errors="coerce", format="mixed").dt.year; rows = []; az = None
for drug, rule in RULES.items():
    obs = a.ast.map(lambda d: d.get(drug)); ok = obs.isin(["S", "R"])
    if ok.sum() < 30: continue
    s = a[ok].copy(); s["pred"] = s.genes.map(lambda g: rule(g, cat)); s["y"] = obs[ok] == "R"; fn = s[~s.pred & s.y]; fp = s[s.pred & ~s.y]
    rows.append(dict(drug=drug, n=len(s), n_R=int(s.y.sum()), fn=len(fn), fn_no_core_determinant=int((fn.genes.map(len) == 0).sum()), fn_other_determinants_only=int((fn.genes.map(len) > 0).sum()), fp=len(fp),
                     fn_top_serovar=fn.serovar.value_counts().index[0] if len(fn) else "", fn_top_serovar_share=round(float(fn.serovar.value_counts().iloc[0] / len(fn)), 3) if len(fn) else np.nan))
    if drug == "azithromycin": az = (s, fn, fp)
D = pd.DataFrame(rows); D.to_csv("table_discordance.csv", index=False); print(D.to_string(index=False))
tot_fn = int(D.fn.sum()); tot_none = int(D.fn_no_core_determinant.sum()); print(f"\nall drugs: {tot_fn} false negatives, {tot_none} ({100 * tot_none / tot_fn:.1f}%) with no core determinant at all")
s, fn, fp = az; s["period"] = pd.cut(s.year, [1999, 2015, 2019, 2024], labels=["2000-2015", "2016-2019", "2020-2024"])
P = s.groupby("period", observed=True).apply(lambda d: pd.Series(dict(n=len(d), n_R=int(d.y.sum()), sensitivity=(d.pred & d.y).sum() / max(d.y.sum(), 1), tp=int((d.pred & d.y).sum()), fn=int((~d.pred & d.y).sum()))), include_groups=False).reset_index(); P.round(4).to_csv("table_azithromycin_period.csv", index=False); print(P.round(3).to_string(index=False))
mut = pd.read_csv("amrf_db_catalog/AMRProt-mutation.tsv", sep="\t"); acr = mut[(mut["#taxgroup"] == "Salmonella") & mut.standard_mutation_symbol.str.startswith("acrB_R717")]
carry = int(fn.genes.map(lambda g: any(x.startswith("acrB_R717") for x in g)).sum())
print(f"\nazithromycin: {len(fn)} resistant isolates without a macrolide determinant; the catalogue lists {len(acr)} Salmonella acrB R717 variants ({', '.join(acr.standard_mutation_symbol)}) and {carry} of the {len(fn)} carry one;"
      f" the {len(fn)} carry on median {int(fn.genes.map(len).median())} other determinants; {int((fn.serovar == 'Infantis').sum())} are Infantis; years {int(fn.year.min())}-{int(fn.year.max())};"
      f" {len(fp)} susceptible isolates carry a macrolide determinant")
pd.DataFrame([dict(drug="azithromycin", resistant_without_macrolide_determinant=len(fn), acrB_R717_variants_in_catalogue=len(acr), false_negatives_carrying_one=carry, fn_database_versions=";".join(f"{k}:{v}" for k, v in fn.refgene_db_version.value_counts().items()), median_other_determinants=int(fn.genes.map(len).median()), infantis=int((fn.serovar == "Infantis").sum()), first_year=int(fn.year.min()), last_year=int(fn.year.max()), susceptible_with_macrolide_determinant=len(fp))]).to_csv("table_azithromycin_detail.csv", index=False)
