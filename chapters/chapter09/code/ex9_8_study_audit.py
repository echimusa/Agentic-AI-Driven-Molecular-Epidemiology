"""
Chapter 9, reading cross-sector studies - recompute the figures printed in three of the sources of this chapter from the counts around them.
  Klaas et al. (2026) Appl Environ Microbiol 92(2):e01905-25, doi 10.1128/aem.01905-25   dairy farms, wastewater and E. coli, western Canada
  Kitamura et al. (2026) Nat Commun, doi 10.1038/s41467-026-78530-9 (article in press)  Salmonella, Campylobacter and Enterococcus in Japan
  McDermott (course slides)                                                              genomes with phenotype data in the US NARMS programme
  Muloi et al. (2023) Lancet Microbe 4:e1056-e1062, doi 10.1016/S2666-5247(23)00284-7   colistin and mcr-1 screening
Nothing is downloaded: the counts are typed from the papers and the slides. Writes table_study_audit.csv with the reported value, the recomputed value and whether they agree to the printed
precision (tolerance 0.05 for a percentage, exact for a count, 0.5 point where the paper says 'approximately').
"""
import math
import pandas as pd
rows = []
def add(study, item, reported, recomputed, tol, note=""):
    rows.append(dict(study=study, item=item, reported=reported, recomputed=round(float(recomputed), 4), agrees=bool(abs(reported - recomputed) <= tol), note=note))
K = "Klaas 2026"; N = "Kitamura 2026"; M = "McDermott slides"; L = "Muloi 2023"
# ---------------------------------------------------------------- Klaas et al. 2026
tot_src = dict(calf=109, milking=219, dry=207, pit=152, soil=76, water=179, influent=30, effluent=44); sub_src = dict(calf=60, milking=60, dry=60, pit=60, soil=61, water=60, influent=30, effluent=30)
add(K, "isolates recovered: the eight sources add up to 1,016", 1016, sum(tot_src.values()), 0)
add(K, "isolates sequenced: the eight sources add up to 421", 421, sum(sub_src.values()), 0)
add(K, "susceptible to all fourteen drugs (%)", 82.19, 100 * 346 / 421, 0.05, "346 of 421")
add(K, "resistant isolates = 421 - 346", 75, 421 - 346, 0)
add(K, "tetracycline resistance (%)", 13.06, 100 * 55 / 421, 0.05, "55 of 421")
add(K, "ampicillin resistance (%)", 10.45, 100 * 44 / 421, 0.05, "44 of 421")
add(K, "resistant to more than one drug, of the 75 (%)", 72.00, 100 * 54 / 75, 0.05, "54 of 75")
add(K, "phenotypic MDR, of the 75 (%)", 45.33, 100 * 34 / 75, 0.05, "34 of 75")
add(K, "isolates with at least one ARG (%)", 17.8, 100 * 75 / 421, 0.05, "75 of 421")
add(K, "genotypic MDR, of the 75 with an ARG (%)", 57.33, 100 * 43 / 75, 0.05, "43 of 75")
add(K, "genotypic MDR, of all 421 (%), as in the abstract", 10.2, 100 * 43 / 421, 0.05, "43 of 421")
add(K, "isolates with a beta-lactam gene (%), as in the abstract", 9.5, 100 * 40 / 421, 0.05, "40 of 421")
add(K, "isolates with a plasmid replicon (%)", 72.4, 100 * 305 / 421, 0.05, "305 of 421; the paper gives the count only")
pairs = math.comb(421, 2)
add(K, "pairwise comparisons exceed 88,000", 88000, min(pairs, 88000), 0, f"{pairs:,} pairs")
add(K, "pairs within 100 SNPs, of all pairs (%): less than 0.25%", 0.25, max(0.25, 100 * 207 / pairs), 0.0005, f"207 of {pairs:,} = {100 * 207 / pairs:.3f}%")
add(K, "close pairs from the same sample source (%)", 55.07, 100 * 114 / 207, 0.05, "114 of 207")
add(K, "close pairs from different One Health domains (%): the paper says approximately 25%", 25.0, 100 * 46 / 207, 0.5, "46 of 207")
add(K, "close pairs: same source + different domains + same domain, different source = 207", 207, 114 + 46 + 47, 0)
add(K, "isolates with a close relative (%)", 50.59, 100 * 213 / 421, 0.05, "213 of 421")
add(K, "isolates with a close relative in another source (%)", 27.08, 100 * 114 / 421, 0.05, "114 of 421")
# ---------------------------------------------------------------- Kitamura et al. 2026
add(N, "Salmonella isolates = human 623 + food 563 + animal 427", 1613, 623 + 563 + 427, 0)
add(N, "Campylobacter isolates = 261 + 57 + 364", 682, 261 + 57 + 364, 0, "the split between sectors is read from Table 1")
add(N, "Enterococcus isolates = hospital 38 + chicken meat 34", 72, 38 + 34, 0)
# ---------------------------------------------------------------- McDermott slides: genomes with phenotype data in NARMS
H = dict(Salmonella=8629, Campylobacter=1041, Ecoli=617); FP = dict(Salmonella=7206, Campylobacter=9571, Ecoli=3079, Enterococcus=1772); RM = dict(Salmonella=9041, Campylobacter=4770, Ecoli=3160, Enterococcus=325); DV = dict(Salmonella=827, Ecoli=842, Spseudintermedius=855)
tot = dict(Salmonella=25703, Campylobacter=15382, Ecoli=7698, Enterococcus=2097, Spseudintermedius=855)
for b in tot: add(M, f"{b}: the column adds up", tot[b], H.get(b, 0) + FP.get(b, 0) + RM.get(b, 0) + DV.get(b, 0), 0)
add(M, "human genomes = 8,629 + 1,041 + 617", 10287, sum(H.values()), 0); add(M, "food-producing animal genomes", 21628, sum(FP.values()), 0); add(M, "retail meat genomes", 17296, sum(RM.values()), 0); add(M, "diagnostic laboratory genomes", 2524, sum(DV.values()), 0)
add(M, "all genomes with phenotype data", 51735, sum(tot.values()), 0); add(M, "all genomes, summed over the four sources", 51735, 10287 + 21628 + 17296 + 2524, 0)
# ---------------------------------------------------------------- Muloi et al. 2023
add(L, "isolates screened for mcr-1 and sequenced (%)", 4.8, 100 * 135 / 2824, 0.05, "135 of 2,824; the paper gives the counts only")
A = pd.DataFrame(rows); A.to_csv("table_study_audit.csv", index=False)
S = A.groupby("study", sort=False).agg(figures=("agrees", "size"), agree=("agrees", "sum")).reset_index(); S["disagree"] = S.figures - S.agree; S.to_csv("table_study_audit_summary.csv", index=False)
print(A[~A.agrees].to_string(index=False)); print(S.to_string(index=False)); print(f"{int(A.agrees.sum())} of {len(A)} figures agree")
