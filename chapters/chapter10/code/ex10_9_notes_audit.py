"""
Chapter 10, Hands-on 9 - recompute the worked examples printed in the course sources.
The lecture notes on population structure, GWAS and admixture, the genetic-epidemiology slides of K. Van Steen and the student slides on GWAS print small worked examples. This script recomputes each printed
number from the counts around it and records whether it agrees to the printed precision (half a unit of the last printed digit, or the stated tolerance).  Nothing is downloaded; the inputs are typed from the sources.
Output: table_notes_audit.csv, table_notes_audit_summary.csv.
"""
import math
import numpy as np, pandas as pd
from scipy import stats
rows = []
def add(src, item, printed, recomputed, tol, note=""): rows.append(dict(source=src, item=item, printed=printed, recomputed=round(float(recomputed), 6), agrees=bool(abs(printed - recomputed) <= tol), tolerance=tol, note=note))
N2, N3, N5, VS, GB = "Notes: population structure", "Notes: GWAS and stratification", "Notes: admixture mapping", "Van Steen slides", "Student GWAS slides"
# ---------------------------------------------------------------- Hardy-Weinberg example: 298 MM, 489 MN, 213 NN
mm, mn, nn = 298, 489, 213; n = mm + mn + nn; p = (2 * mm + mn) / (2 * n); q = 1 - p
add(N2, "frequency of allele M", 0.54, p, 0.005, "0.5425 rounds to 0.54"); se = math.sqrt(p * q / (2 * n))
add(N2, "95% interval, lower limit", 0.52, p - 1.96 * se, 0.005); add(N2, "95% interval, upper limit", 0.56, p + 1.96 * se, 0.005)
for g, o, e in (("MM", mm, p * p * n), ("MN", mn, 2 * p * q * n), ("NN", nn, q * q * n)): add(N2, f"expected count {g}", {"MM": 294.3, "MN": 496.4, "NN": 209.3}[g], e, 0.05)
D = mm / n - p * p; add(N2, "D = P(MM) - p^2", 0.004, D, 0.0005)
chi_notes = n * (0.004 / (0.54 * 0.46)) ** 2; chi = sum((o - e) ** 2 / e for o, e in ((mm, p * p * n), (mn, 2 * p * q * n), (nn, q * q * n)))
add(N2, "chi-square printed (from the rounded D and p)", 0.25, chi_notes, 0.01, "0.259 from the rounded inputs: the printed 0.25 is truncated"); add(N2, "chi-square from the observed and expected counts", 0.25, chi, 0.05, "0.22 from unrounded values; the notes round D and p first")
add(N2, "P value of the printed chi-square", 0.61, stats.chi2.sf(0.25, 1), 0.01, "0.617")
# ---------------------------------------------------------------- linkage disequilibrium example
pA, pB, pAB = 0.25, 0.40, 0.25; Dl = pAB - pA * pB; r2 = Dl ** 2 / (pA * (1 - pA) * pB * (1 - pB)); dmax = min(pA * (1 - pB), (1 - pA) * pB)
add(N2, "r2 for pA 0.25, pB 0.4, pAB 0.25", 0.5, r2, 0.005); add(N2, "D' for the same haplotype frequencies", 1.0, Dl / dmax, 0.005)
add(N2, "between-population share of variation, 5-7%, against 93-95% within (Rosenberg et al.)", 100.0, 7 + 93, 0.0, "the two ranges are complements")
# ---------------------------------------------------------------- inflation factor
med = stats.chi2.ppf(0.5, 1); add(N3, "median of the chi-square distribution with 1 df", 0.455, med, 0.0005, "0.4549")
add(N3, "T2D example: lambda from median chi-square 0.48", 1.05, 0.48 / med, 0.011, "1.055; the median is itself rounded to two digits, which allows +/- 0.011"); add(N3, "T2D example: P value at the median chi-square 0.48", 0.49, stats.chi2.sf(0.48, 1), 0.005); add(N3, "T2D example: -log10 P at the median", 0.31, -math.log10(stats.chi2.sf(0.48, 1)), 0.005)
add(N3, "lactase SNP: P value of chi-square 25.9 (printed 4 x 10^-7)", 4e-7, stats.chi2.sf(25.9, 1), 0.5e-7, "3.6 x 10^-7")
add(N3, "lactase SNP: corrected chi-square 25.9 / 1.11", 23.3, 25.9 / 1.11, 0.05); add(N3, "lactase SNP: corrected P (printed 1 x 10^-6)", 1e-6, stats.chi2.sf(25.9 / 1.11, 1), 0.5e-6, "1.4 x 10^-6")
add(N3, "lactase SNP: 1,123 tall + 1,056 short individuals", 2179, 1123 + 1056, 0)
# ---------------------------------------------------------------- admixture toy example and ancestry odds ratios
R = 1.5; ctl = {("A", "POP1"): 0.05, ("A", "POP2"): 0.45, ("C", "POP1"): 0.45, ("C", "POP2"): 0.05}
cas = {k: v * (R if k[0] == "A" else 1.0) for k, v in ctl.items()}; z = sum(cas.values()); cas = {k: v / z for k, v in cas.items()}
printed = {("A", "POP1"): 0.06, ("A", "POP2"): 0.54, ("C", "POP1"): 0.36, ("C", "POP2"): 0.04}
for k in ctl: add(N5, f"case probability of SNP1={k[0]} with {k[1]} ancestry", printed[k], cas[k], 0.005)
add(N5, "P(SNP1 = A) in cases", 0.60, cas[("A", "POP1")] + cas[("A", "POP2")], 0.005); add(N5, "P(POP2) in cases", 0.58, cas[("A", "POP2")] + cas[("C", "POP2")], 0.005)
for pE, pAf, om in ((1.0, 0.0, 1.50), (0.0, 1.0, 0.67), (0.9, 0.1, 1.38), (0.75, 0.25, 1.22), (0.6, 0.4, 1.08), (0.5, 0.5, 1.00)): add(N5, f"ancestry odds ratio, pE {pE:.0%}, pA {pAf:.0%}", om, (R * pE + 1 - pE) / (R * pAf + 1 - pAf), 0.005)
# ---------------------------------------------------------------- familial risk ratios and the first GWAS
add(VS, "cystic fibrosis: sibling risk 0.25 over population risk 0.0004", 500, 0.25 / 0.0004, 0.5, "625; 500 would need a population risk of 0.0005"); add(VS, "Huntington disease: 0.5 over 0.0001", 5000, 0.5 / 0.0001, 0.5)
add(GB, "WTCCC: associations in the seven diseases, as the slide counts them by disease (0, 1, 9, 7, 3, 4) against its stated total of 23", 23, 0 + 1 + 9 + 7 + 3 + 4, 0, "the slide's own parts add up to 24; the paper reports 24 signals at P < 5 x 10^-7 (1, 1, 9, 3, 7, 3)")
add(GB, "WTCCC paper: signals by disease (bipolar 1, coronary artery 1, Crohn 9, rheumatoid 3, type 1 diabetes 7, type 2 diabetes 3)", 24, 1 + 1 + 9 + 3 + 7 + 3, 0)
A = pd.DataFrame(rows); A.to_csv("table_notes_audit.csv", index=False); S = A.groupby("source", sort=False).agg(figures=("agrees", "size"), agree=("agrees", "sum")).reset_index(); S["disagree"] = S.figures - S.agree; S.to_csv("table_notes_audit_summary.csv", index=False)
print(A[~A.agrees][["source", "item", "printed", "recomputed", "note"]].to_string(index=False)); print(S.to_string(index=False)); print(f"{int(A.agrees.sum())} of {len(A)} figures agree")
