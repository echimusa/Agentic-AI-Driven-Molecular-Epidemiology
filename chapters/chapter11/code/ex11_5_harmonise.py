"""
Chapter 11, Hands-on 5 - why alleles must be harmonised (the third question of the Session 11 quiz).
The 927 body mass index variants of Hands-on 1 are found in both files. Four ways of preparing them are compared on the same data:
  A  no harmonisation: the outcome effect is taken as it stands, as if both studies reported the same allele
  B  align alleles that are listed in opposite order, but keep every palindromic (A/T or C/G) variant as it stands
  C  align by allele order and align palindromic variants by allele frequency, keeping those with a frequency near 0.5
  D  the rule of Hands-on 1: as C, and drop palindromic variants whose frequency is within 0.42-0.58 in either study
Output: table_harmonise.csv (the estimate under each way), table_harmonise_counts.csv (what the alleles of the 927 variants look like).  Needs numpy, pandas, scipy and mr11.py; the cached inputs of Hands-on 1.
"""
import numpy as np, pandas as pd
import mr11 as mr
X = pd.read_csv("bmi_exposure_snps.csv").rename(columns={"Allele1": "ea", "Allele2": "oa", "Freq1": "eaf_x", "Effect": "bx", "StdErr": "sx", "P-value": "px"}); X[["ea", "oa"]] = X[["ea", "oa"]].apply(lambda c: c.str.upper())
Y = pd.read_csv("hgi_for_bmi.csv"); D = X.drop_duplicates("rsid").merge(Y.drop_duplicates("rsid"), on="rsid"); print(f"{len(D)} variants in both files")
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}; pal = (D.ea.map(comp) == D.oa).values; same = ((D.ea == D.oa1) & (D.oa == D.oa2)).values; flip = ((D.ea == D.oa2) & (D.oa == D.oa1)).values
amb = pal & (((D.eaf_x > 0.42) & (D.eaf_x < 0.58)) | ((D.eaf_y > 0.42) & (D.eaf_y < 0.58))).values
counts = pd.DataFrame([dict(category="same allele order", variants=int((same & ~pal).sum())), dict(category="opposite allele order, not palindromic", variants=int((flip & ~pal).sum())), dict(category="palindromic (A/T or C/G)", variants=int(pal.sum())),
                       dict(category="  of which frequency within 0.42 to 0.58 in either study", variants=int(amb.sum())), dict(category="not matched by alleles in either order or strand", variants=int((~(same | flip | pal)).sum())), dict(category="all variants in both files", variants=len(D))])
counts.to_csv("table_harmonise_counts.csv", index=False); print(counts.to_string(index=False))
bx, sx, by0, sy = D.bx.values, D.sx.values, D.by.values, D.sy.values
def est(mask, sign, label):
    b, s = bx[mask], sx[mask]; y = by0[mask] * sign[mask]; e = mr.ivw(b, y, sy[mask]); se = e["se_random"]
    return dict(way=label, variants=int(mask.sum()), log_or=e["beta"], se=se, odds_ratio=np.exp(e["beta"]), ci_low=np.exp(e["beta"] - 1.96 * se), ci_high=np.exp(e["beta"] + 1.96 * se), mean_F=float(np.mean(mr.fstat(b, s))),
                share_same_sign_as_estimate=float(np.mean(np.sign(y / b) == np.sign(e["beta"]))))
one = np.ones(len(D)); sgn_B = np.where(flip & ~pal, -1, 1); keep, sgn_D = mr.harmonise(D.ea.values, D.oa.values, D.eaf_x.values, D.oa1.values, D.oa2.values, D.eaf_y.values)
pal_sign = np.where((D.eaf_x > 0.5).values == (D.eaf_y > 0.5).values, 1, -1); sgn_C = np.where(pal, pal_sign, np.where(flip, -1, 1)); ok = same | flip | pal
R = pd.DataFrame([est(np.ones(len(D), bool), one, "A  no harmonisation"), est(ok, sgn_B, "B  align order, keep palindromes as they stand"), est(ok, sgn_C, "C  align by order and by frequency, keep all palindromes"), est(keep, sgn_D, "D  rule of Hands-on 1 (drop ambiguous palindromes)")])
R.to_csv("table_harmonise.csv", index=False); print(R.round(3).to_string(index=False))
