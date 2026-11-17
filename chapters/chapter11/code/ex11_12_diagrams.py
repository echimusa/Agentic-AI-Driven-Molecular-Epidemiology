"""
Chapter 11 - the two new drawn diagrams (no data): five ways an instrument can fail with the check for each, and the four causal estimands that an instrumental-variable estimate can identify (Table 2 of Yao et al. 2025, in my words).
Output: figs11/diagram_failures_and_checks.png, figs11/diagram_estimands.png
"""
import os, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
os.makedirs("figs11", exist_ok=True); INK = "#1a202c"; GREY = "#4a5568"; EDGE = "#c8c5bb"
def box(ax, x, y, w, h, fc="white", ec=EDGE, lw=1.4): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=lw))
W = lambda t, n: "\n".join(textwrap.wrap(t, n))
fig, ax = plt.subplots(figsize=(12.4, 5.9)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 5.9)
ax.text(0.15, 5.68, "An instrument can fail in five ways, and each has its own check", fontsize=13.5, fontweight="bold", color=INK, va="center")
rows = [("1  Weak", "The variants barely move the exposure (relevance)", "F-statistic; debiased estimators; larger exposure GWAS", "#1B998B"), ("2  Confounded", "Ancestry, assortative mating or family effects link the variant to the outcome (independence)", "Principal components, mixed models, within-family design", "#E4572E"),
        ("3  Pleiotropic", "The variant reaches the outcome by another route, or is in LD with a variant that does (exclusion)", "Q, Egger intercept, median, mode, outlier tests, multivariable MR", "#F3A712"), ("4  Backwards", "The outcome or a feedback loop moves the exposure, or the variants really act on the outcome", "Steiger test, bidirectional MR, choose instruments from the right GWAS", "#3D7EA6"),
        ("5  Selected or overlapping", "Variants chosen and estimated in one sample; exposure and outcome samples share people", "Select in an independent sample; overlap-robust methods", "#2E294E")]
for i, (t, d, c, col) in enumerate(rows):
    y = 4.35 - i * 0.86; box(ax, 0.15, y, 12.0, 0.74); ax.add_patch(FancyBboxPatch((0.15, y), 0.12, 0.74, boxstyle="round,pad=0.0,rounding_size=0.05", fc=col, ec=col))
    ax.text(0.5, y + 0.37, t, fontsize=10.6, fontweight="bold", color=INK, va="center"); ax.text(3.45, y + 0.37, W(d, 50), fontsize=9.4, color=INK, va="center", linespacing=1.25); ax.text(8.5, y + 0.37, W(c, 43), fontsize=9.4, color=GREY, va="center", linespacing=1.25)
ax.text(3.45, 5.3, "What goes wrong", fontsize=9.2, color=GREY, va="center", style="italic"); ax.text(8.5, 5.3, "How to check", fontsize=9.2, color=GREY, va="center", style="italic")
fig.savefig("figs11/diagram_failures_and_checks.png", dpi=200, bbox_inches="tight"); plt.close(fig)
fig, ax = plt.subplots(figsize=(12.4, 4.9)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.9)
ax.text(0.15, 4.6, "The same instrumental-variable ratio can answer four different causal questions", fontsize=13.5, fontweight="bold", color=INK, va="center")
cards = [("Constant effect", "The effect of one unit of the exposure, the same for everyone.", "Needs: an effect that is the same in every person and linear (the ALICE model)", "#3D7EA6"), ("Effect on the treated", "The effect among people who actually have the higher exposure.", "Needs: the effect does not depend on the instrument among those people (additive homogeneity)", "#1B998B"),
         ("Local effect", "The effect among people whose exposure the variant actually changes (compliers).", "Needs: the variant never lowers the exposure of anyone (monotonicity)", "#F3A712"), ("Average effect in everyone", "The effect of moving everyone's exposure by one unit.", "Needs: no interaction between the effect and the instrument, plus control of confounding", "#E4572E")]
for i, (t, d, m, col) in enumerate(cards):
    x = 0.15 + (i % 2) * 6.1; y = 2.4 - (i // 2) * 2.05; box(ax, x, y, 5.9, 1.85); ax.add_patch(FancyBboxPatch((x, y + 1.85 - 0.1), 5.9, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col))
    ax.text(x + 0.2, y + 1.45, t, fontsize=11.5, fontweight="bold", color=INK, va="center"); ax.text(x + 0.2, y + 0.95, W(d, 62), fontsize=9.6, color=INK, va="center", linespacing=1.4); ax.text(x + 0.2, y + 0.3, W(m, 66), fontsize=8.9, color=GREY, va="center")
fig.savefig("figs11/diagram_estimands.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("diagrams written")
