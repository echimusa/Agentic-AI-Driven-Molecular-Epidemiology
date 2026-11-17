"""
Chapter 13, Hands-on 8 - Two databases, one outbreak: an audit of the published numbers of the Rhode Island HIV study (Khanna et al. 2025).
The study joins a contact-tracing database (CTDB) and a database of HIV-1 sequences (GDB). The script recomputes the percentages and the Jaccard coefficients from the counts printed in the paper, flags the places where
the abstract and the results use different counts, derives the counts implied by the abstract's cascade percentages, and works out what overlap two independent samples of one population would give.
The counts are those printed in the article (abstract, results and Table 1); no data are redistributed. Output: table_overlap_audit.csv, table_cascade_derived.csv, table_independence.csv, figs13/fig8_integration.png.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
CTDB, CTDB_abs, GDB, BOTH = 2417, 2418, 2527, 894              # persons in the contact-tracing database (results text, abstract), in the sequence database, in both
LCT, LGDB, LBOTH = 1373, 3526, 59                              # partner-naming links, molecular links, links in both
jc = lambda a, b, ab: ab / (a + b - ab)
checks = [("Persons in both databases as a share of CTDB persons", BOTH, CTDB, "37%", "results"), ("Persons in both databases as a share of GDB persons", BOTH, GDB, "35%", "results"), ("Contact-tracing links also found as molecular links", LBOTH, LCT, "4%", "results"), ("Molecular links also found as contact-tracing links", LBOTH, LGDB, "2%", "results"),
          ("GDB persons in molecular clusters", 1204, GDB, "48%", "abstract, results"), ("Persons in both databases who cluster", 583, BOTH, "65%", "results"), ("Index cases in both databases who named a partner", 494, BOTH, "55%", "results"),
          ("Partners named by these 494: contacts attempted", 687, 1056, "65%", "results"), ("Of attempted: reached", 631, 687, "92%", "results"), ("Of reached: HIV-tested", 446, 631, "71%", "results"), ("Reached, not tested and already diagnosed", 24, 185, "13%", "results"), ("Of those tested: HIV-diagnosed", 121, 446, "27%", "results")]
rows = []
for lab, a, b, pub, src in checks: v = 100 * a / b; rows.append(dict(quantity=lab, numerator=a, denominator=b, published=pub, recomputed_pct=round(v, 2), source=src))
rows.append(dict(quantity="Jaccard coefficient of persons", numerator=BOTH, denominator=CTDB + GDB - BOTH, published="0.22", recomputed_pct=round(jc(CTDB, GDB, BOTH), 4), source="abstract, results"))
rows.append(dict(quantity="Jaccard coefficient of persons, with the abstract's 2,418", numerator=BOTH, denominator=CTDB_abs + GDB - BOTH, published="0.22", recomputed_pct=round(jc(CTDB_abs, GDB, BOTH), 4), source="abstract"))
rows.append(dict(quantity="Jaccard coefficient of links", numerator=LBOTH, denominator=LCT + LGDB - LBOTH, published="0.012", recomputed_pct=round(jc(LCT, LGDB, LBOTH), 4), source="abstract, results"))
A = pd.DataFrame(rows); A.to_csv("table_overlap_audit.csv", index=False); print(A.to_string(index=False))
# the abstract states a cascade for the 1,342 named partners with different percentages; the implied counts are derived here and are not printed in the paper
named = 1342; steps = [("Named partners", 1.0), ("Contacts attempted", 0.66), ("Reached", 0.93), ("Newly HIV-tested", 0.71), ("Newly diagnosed", 0.27)]; n = named; rows = []
for lab, p in steps: n = n * p if lab != "Named partners" else named; rows.append(dict(step=lab, percent_of_previous_step=100 * p if lab != "Named partners" else np.nan, implied_count=n))
C = pd.DataFrame(rows); C.to_csv("table_cascade_derived.csv", index=False); print(C.round(1).to_string(index=False))
# overlap expected if two databases were independent samples of one population
N_ind = CTDB * GDB / BOTH; c_eq = 2 * jc(CTDB, GDB, BOTH) / (1 + jc(CTDB, GDB, BOTH)); rows = [dict(quantity="Population size for which independent samples of 2,417 and 2,527 persons overlap by 894", value=N_ind), dict(quantity="Coverage of CTDB persons in that population", value=CTDB / N_ind), dict(quantity="Coverage of GDB persons in that population", value=GDB / N_ind),
                                                                                              dict(quantity="Equal coverage c that gives a Jaccard coefficient of 0.22 (c / (2 - c) = JC)", value=c_eq)]
I = pd.DataFrame(rows); I.to_csv("table_independence.csv", index=False); print(I.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]})
st = ["Named by 494\nindex cases", "Attempted", "Reached", "HIV-tested", "Newly\ndiagnosed"]; vals = [1056, 687, 631, 446, 121]
ax[0].bar(range(5), vals, color="#2b6cb0"); [ax[0].text(i, v + 20, f"{v:,}", ha="center", fontsize=9) for i, v in enumerate(vals)]; ax[0].set(xticks=range(5), xticklabels=st, ylabel="Partners", title="A. Partner naming cascade"); ax[0].tick_params(axis="x", labelsize=8.5)
ax[1].bar([0, 1], [CTDB - BOTH, GDB - BOTH], bottom=0, color=["#a0aec0", "#d69e2e"], label="only one database"); ax[1].bar([0, 1], [BOTH, BOTH], bottom=[CTDB - BOTH, GDB - BOTH], color="#2f855a", label="in both")
ax[1].bar([2, 3], [LCT - LBOTH, LGDB - LBOTH], color=["#a0aec0", "#d69e2e"]); ax[1].bar([2, 3], [LBOTH, LBOTH], bottom=[LCT - LBOTH, LGDB - LBOTH], color="#2f855a"); ax[1].set(xticks=range(4), xticklabels=["CTDB\npersons", "GDB\npersons", "CTDB\nlinks", "GDB\nlinks"], ylabel="Persons or links", title="B. Overlap: persons 0.22, links 0.012"); ax[1].legend(frameon=False, fontsize=8.8)
cl = [48, 65, 71, 88]; ax[2].bar(range(4), cl, color=["#a0aec0", "#2b6cb0", "#2b6cb0", "#c53030"]); [ax[2].text(i, v + 1.5, f"{v}%", ha="center", fontsize=9) for i, v in enumerate(cl)]; ax[2].set(xticks=range(4), xticklabels=["All GDB\npersons", "In both\ndatabases", "Gave partner\ndata", "Named partners\nin both"], ylim=(0, 100), ylabel="Persons in molecular clusters (%)", title="C. Clustering rises with overlap"); ax[2].tick_params(axis="x", labelsize=8.5)
fig.tight_layout(); fig.savefig("figs13/fig8_integration.png", dpi=200); plt.close(fig); print("figure written")
