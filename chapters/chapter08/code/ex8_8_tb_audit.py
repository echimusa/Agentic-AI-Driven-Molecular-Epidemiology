"""
Chapter 8, reading a published study - recompute the numbers of Li et al. (2024), Int J Infect Dis 140:124-131 (doi 10.1016/j.ijid.2023.10.015), from the figures printed in its text.
Nothing is downloaded: the counts below are typed from the paper (Results, Table 1 and Methods). The script recomputes each reported percentage and mean from the counts around it
and tests the main comparison, and writes table_tb_audit.csv with the reported value, the recomputed value and whether they agree to the printed precision.
"""
import numpy as np, pandas as pd
from scipy import stats
TOTAL = 3204; SHANDONG_NEW = 1449; PUBLIC = 1755
rows = []
def add(item, reported, recomputed, note, tol):
    rows.append(dict(item=item, reported=reported, recomputed=round(recomputed, 4), agrees=bool(abs(reported - recomputed) <= tol), note=note))
add("New genomes + public genomes = all genomes", TOTAL, SHANDONG_NEW + PUBLIC, "1,449 + 1,755", 0)
add("Share of all genomes from Shandong (%)", np.nan, 100 * SHANDONG_NEW / TOTAL, "not stated in the paper", np.inf)
add("Lineage 2 share of all genomes (%)", 85.7, 100 * 2745 / TOTAL, "2,745 of 3,204", 0.05)
add("Sublineage 2.2.1 share of all genomes (%)", 80.9, 100 * 2591 / TOTAL, "2,591 of 3,204", 0.05)
add("Sublineage 2.2.1 share of lineage 2 (%)", 94.4, 100 * 2591 / 2745, "2,591 of 2,745", 0.05)
add("Lineage 4 share of all genomes (%)", 13.8, 100 * 443 / TOTAL, "443 of 3,204", 0.05)
add("Pansusceptible share (%)", 45.4, 100 * 1455 / TOTAL, "1,455 of 3,204", 0.05)
add("MDR share of all genomes (%)", 27.6, 100 * 1315 / TOTAL, "1,315 of 3,204", 0.05)
add("Pre-XDR share of all genomes (%)", 13.4, 100 * 430 / TOTAL, "430 of 3,204", 0.05)
add("Lineage 2 share of MDR (%)", 86.8, 100 * 1142 / 1315, "1,142 of 1,315", 0.05)
add("Lineage 4 share of MDR (%)", 12.8, 100 * 168 / 1315, "168 of 1,315", 0.05)
add("Pyrazinamide resistance (%)", 12.8, 100 * 339 / TOTAL, "339 of 3,204", 0.05)
add("Ethambutol resistance (%)", 12.8, 100 * 855 / TOTAL, "855 of 3,204", 0.05)
add("Streptomycin resistance (%)", 12.8, 100 * 1101 / TOTAL, "1,101 of 3,204", 0.05)
add("Mean cluster size at 5 SNPs", 3.05, 1384 / 471, "1,384 strains in 471 clusters", 0.005)
add("Mean cluster size at 10 SNPs", 3.17, 1596 / 499, "1,596 strains in 499 clusters", 0.005)
add("Mean cluster size at 12 SNPs", 3.18, 1631 / 504, "1,631 strains in 504 clusters", 0.005)
add("Clustered share, lineage 2 (%)", 53.3, 100 * 1462 / 2745, "1,462 of 2,745", 0.05)
add("Clustered share, sublineage 2.2.1 (%)", 54.0, 100 * 1400 / 2591, "1,400 of 2,591", 0.05)
add("Clustered share, lineage 4 (%)", 29.8, 100 * 132 / 443, "132 of 443", 0.05)
add("Clustered strains in all (n)", 1596, 1462 + 132 + 2, "L2 1,462 + L4 132 + L3 2 (L1 0)", 0)
add("Regional + cross-regional clusters = all clusters", 499, 413 + 86, "413 + 86", 0)
# a count cannot exceed the number it is a part of: MDR isolates said to belong to sublineage 4.8, which has 3 isolates in all
rows.append(dict(item="MDR isolates assigned to L4.8 (cannot exceed the 3 L4.8 isolates)", reported=80, recomputed=3, agrees=bool(80 <= 3), note="80 reported; L4.8 n = 3 (L4.4 has 200)"))
A = pd.DataFrame(rows); A.to_csv("table_tb_audit.csv", index=False)
chi = stats.chi2_contingency([[1462, 2745 - 1462], [132, 443 - 132]], correction=False); T = pd.DataFrame([dict(test="clustered, lineage 2 vs lineage 4 (Pearson chi-square, no correction)", chi2=chi[0], p=chi[1], reported="P < 0.001")]); T.to_csv("table_tb_chi2.csv", index=False)
print(A.to_string(index=False)); print(T.round(6).to_string(index=False))
