"""
Chapter 10 - the two new drawn diagrams (no data): the ladder of questions from families to variants, and the four words of population genetics (structure, stratification, admixture, admixture mapping).
Output: figs10/diagram_ladder_of_questions.png, figs10/diagram_structure_stratification_admixture.png
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
os.makedirs("figs10", exist_ok=True); INK = "#1a202c"; GREY = "#4a5568"; EDGE = "#c8c5bb"
def box(ax, x, y, w, h, fc="white", ec=EDGE, lw=1.4): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=lw))
# ---------------------------------------------------------------- the ladder of questions
import textwrap
fig, ax = plt.subplots(figsize=(12.4, 4.7)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.7)
ax.text(0.15, 4.45, "Host genetics climbs a ladder of questions, and each rung has its own design", fontsize=13.5, fontweight="bold", color=INK, va="center")
rungs = [("1  Does it run in families?", "Familial aggregation", "Risk in relatives against the population: the familial risk ratio", "#1B998B"), ("2  Genes or shared environment?", "Twins and adoptees", "Concordance and heritability; deaths from infection in adoptees", "#3D7EA6"),
         ("3  How is it inherited?", "Segregation analysis", "Penetrance; dominant, recessive and additive models", "#F3A712"), ("4  Where is the gene?", "Linkage in families, association in populations", "GWAS tests millions of variants, one at a time", "#E4572E"), ("5  Which variant?", "Fine-mapping and function", "Credible set; functional follow-up", "#2E294E")]
for i, (q, d, m, col) in enumerate(rungs):
    y = 3.55 - i * 0.78; box(ax, 0.15, y, 12.0, 0.66); ax.add_patch(FancyBboxPatch((0.15, y), 0.12, 0.66, boxstyle="round,pad=0.0,rounding_size=0.05", fc=col, ec=col))
    ax.text(0.5, y + 0.33, q, fontsize=10.4, fontweight="bold", color=INK, va="center"); ax.text(4.35, y + 0.33, "\n".join(textwrap.wrap(d, 30)), fontsize=9.6, color=INK, va="center", linespacing=1.3); ax.text(7.7, y + 0.33, "\n".join(textwrap.wrap(m, 46)), fontsize=9.0, color=GREY, va="center", linespacing=1.3)
fig.savefig("figs10/diagram_ladder_of_questions.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# ---------------------------------------------------------------- four words
fig, ax = plt.subplots(figsize=(12.4, 4.9)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.9)
ax.text(0.15, 4.6, "Structure, stratification and admixture are different things", fontsize=13.5, fontweight="bold", color=INK, va="center")
cards = [("Population structure", "Genetic differences between populations because of geographic ancestry.", "Measured by: FST, clusters, principal components", "#3D7EA6"), ("Population stratification", "Ancestry differs between cases and controls, so frequencies differ without any disease effect.", "Detected by: Q-Q plot, λ; handled by components or a mixed model", "#E4572E"),
         ("Admixture", "Recent ancestry from two or more continents; genomes are mosaics of ancestral segments.", "Measured by: global and local ancestry", "#1B998B"), ("Admixture mapping", "A locus where cases carry an unusual share of one ancestry, because a risk allele is more common in that ancestry.", "Tested by: local ancestry in cases against controls, or in cases alone", "#2E294E")]
for i, (t, d, m, col) in enumerate(cards):
    x = 0.15 + (i % 2) * 6.1; y = 2.4 - (i // 2) * 2.05; box(ax, x, y, 5.9, 1.85); ax.add_patch(FancyBboxPatch((x, y + 1.85 - 0.1), 5.9, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col))
    ax.text(x + 0.2, y + 1.45, t, fontsize=11.5, fontweight="bold", color=INK, va="center")
    import textwrap; ax.text(x + 0.2, y + 0.95, "\n".join(textwrap.wrap(d, 62)), fontsize=9.6, color=INK, va="center", linespacing=1.4); ax.text(x + 0.2, y + 0.3, "\n".join(textwrap.wrap(m, 66)), fontsize=8.9, color=GREY, va="center")
fig.savefig("figs10/diagram_structure_stratification_admixture.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("diagrams written")
