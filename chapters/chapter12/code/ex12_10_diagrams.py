"""
Chapter 12 - the two drawn diagrams (no data): three release tiers with what each protects and costs, and the protocol flow that connects the hands-on exercises.
Output: figs12/diagram_release_tiers.png, figs12/diagram_protocol_flow.png
"""
import os, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
os.makedirs("figs12", exist_ok=True); INK = "#1a202c"; GREY = "#4a5568"; EDGE = "#c8c5bb"
def box(ax, x, y, w, h, fc="white", ec=EDGE, lw=1.4): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=lw))
W = lambda t, n: "\n".join(textwrap.wrap(t, n))
# ---- release tiers
fig, ax = plt.subplots(figsize=(12.4, 5.5)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 5.5)
ax.text(0.15, 5.28, "Release is graded: three tiers, and what each protects and costs", fontsize=13.5, fontweight="bold", color=INK, va="center")
tiers = [("Open archive", "#1B998B", "Who: anyone, at once", "Protects: maximum reuse and speed", "Costs: little control after release; metadata must carry minimal identifying detail", "Examples: INSDC archives such as the ENA; Pathoplexus open release"),
         ("Registered use", "#F3A712", "Who: registered users who accept the terms of use", "Protects: names the users and sets terms of use", "Costs: an access step, and terms that set what users may do", "Examples: GISAID; Pathoplexus restricted-use release"),
         ("Controlled access", "#E4572E", "Who: approved applicants under an agreement", "Protects: contributors and privacy", "Costs: slower access; can carry benefit-sharing terms", "Examples: agreement-based access to sensitive or community data")]
for i, (t, col, a, b, c, e) in enumerate(tiers):
    x = 0.15 + i * 4.1; box(ax, x, 0.95, 3.9, 3.95); ax.add_patch(FancyBboxPatch((x, 4.8), 3.9, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col))
    ax.text(x + 0.2, 4.45, t, fontsize=12, fontweight="bold", color=INK, va="center")
    for j, line in enumerate((a, b, c, e)): ax.text(x + 0.2, 3.8 - j * 0.72, W(line, 40), fontsize=9.3, color=INK if j < 3 else GREY, va="center", linespacing=1.3, style="italic" if j == 3 else "normal")
ax.annotate("", xy=(12.2, 0.55), xytext=(0.2, 0.55), arrowprops=dict(arrowstyle="->", color=GREY, lw=1.4)); ax.text(6.2, 0.25, "more control over use, slower access, more governance", fontsize=9.4, color=GREY, ha="center", va="center")
fig.savefig("figs12/diagram_release_tiers.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# ---- protocol flow
fig, ax = plt.subplots(figsize=(12.4, 4.6)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.6)
ax.text(0.15, 4.4, "A protocol assembled from the exercises of this chapter", fontsize=13.5, fontweight="bold", color=INK, va="center")
steps = [("Aim", "R1", "Which of the four aims", "Table of aims", "#1B998B"), ("Target and sizing", "R2, R3", "Variant frequency, confidence, sequences, weeks", "Hands-on 5", "#3D7EA6"), ("Allocation", "R4, R5", "Share of capacity by site against population", "Hands-on 6", "#F3A712"),
         ("Standards", "R6, R7", "Turnaround target, minimum metadata, named standard", "Hands-on 1, 7", "#E4572E"), ("Release", "R8, R9", "Tier, date and place resolution, embargo and credit", "Hands-on 3, 8", "#2E294E"), ("Consent and review", "R10", "Legal basis, credit rule, review date", "People decide", "#6b7280")]
for i, (t, rr, d, h, col) in enumerate(steps):
    x = 0.15 + i * 2.05; box(ax, x, 1.1, 1.9, 2.85); ax.add_patch(FancyBboxPatch((x, 3.85), 1.9, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col))
    ax.text(x + 0.12, 3.5, W(t, 16), fontsize=10.5, fontweight="bold", color=INK, va="center", linespacing=1.2); ax.text(x + 0.12, 2.95, f"Rule {rr}", fontsize=8.8, color=GREY, va="center"); ax.text(x + 0.12, 2.3, W(d, 21), fontsize=8.8, color=INK, va="center", linespacing=1.3); ax.text(x + 0.12, 1.4, h, fontsize=8.8, color=col if col != "#F3A712" else "#b7791f", va="center", fontweight="bold")
    if i < 5: ax.annotate("", xy=(x + 2.03, 2.5), xytext=(x + 1.92, 2.5), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3))
ax.text(0.15, 0.6, "The script checks the arithmetic and the completeness (Hands-on 9); the aim, the consent route and the sharing terms are decided by people.", fontsize=9.6, color=GREY, va="center", style="italic")
fig.savefig("figs12/diagram_protocol_flow.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("diagrams written")
