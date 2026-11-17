"""
Chapter 13 - the two drawn diagrams (no data): the capstone workflow with its four deliverables, and the loop from a question to a decision and its evaluation.
Output: figs13/diagram_capstone_workflow.png, figs13/diagram_evidence_to_action.png
"""
import os, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
os.makedirs("figs13", exist_ok=True); INK = "#1a202c"; GREY = "#4a5568"; EDGE = "#c8c5bb"
def box(ax, x, y, w, h, fc="white", ec=EDGE, lw=1.4): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=lw))
W = lambda t, n: "\n".join(textwrap.wrap(t, n))
# ---- capstone workflow and deliverables
fig, ax = plt.subplots(figsize=(12.4, 5.2)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 5.2)
ax.text(0.15, 4.98, "The capstone: five analysis stages and four deliverables", fontsize=13.5, fontweight="bold", color=INK, va="center")
st = [("QC and typing", "Part A", "#1B998B"), ("Clusters", "Session 4", "#3D7EA6"), ("Trees and dates", "Sessions 5 and 6", "#F3A712"), ("Transmission and spread", "Sessions 7 and 8", "#E4572E"), ("Recommendations", "Part D", "#2E294E")]
for i, (t, s, col) in enumerate(st):
    x = 0.15 + i * 2.45; box(ax, x, 3.0, 2.2, 1.55); ax.add_patch(FancyBboxPatch((x, 4.45), 2.2, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col)); ax.text(x + 0.15, 3.95, W(t, 18), fontsize=10.8, fontweight="bold", color=INK, va="center", linespacing=1.2); ax.text(x + 0.15, 3.3, s, fontsize=9.2, color=GREY, va="center")
    if i < 4: ax.annotate("", xy=(x + 2.43, 3.78), xytext=(x + 2.22, 3.78), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3))
ax.text(0.15, 2.72, "Reuse your portfolio scripts; record every command in scripts/, not only in the terminal", fontsize=9.4, color=GREY, style="italic", va="center")
dl = [("Analysis plan", "Tutors", "Session 13 (formative)", "#1B998B"), ("Team presentation", "Peers and tutors", "Session 15 (formative)", "#3D7EA6"), ("Practical portfolio", "Markers", "One week after the course", "#F3A712"), ("Report: briefing and appendix", "Decision makers and specialists", "Two weeks after the course", "#E4572E")]
for i, (t, a, d, col) in enumerate(dl):
    x = 0.15 + i * 3.05; box(ax, x, 0.35, 2.85, 1.95); ax.add_patch(FancyBboxPatch((x, 2.2), 2.85, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col)); ax.text(x + 0.15, 1.8, W(t, 26), fontsize=10.5, fontweight="bold", color=INK, va="center", linespacing=1.2); ax.text(x + 0.15, 1.2, W(f"Audience: {a}", 28), fontsize=9, color=INK, va="center", linespacing=1.2); ax.text(x + 0.15, 0.62, d, fontsize=9, color=GREY, va="center")
    if i < 3: ax.annotate("", xy=(x + 3.03, 1.3), xytext=(x + 2.87, 1.3), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3))
fig.savefig("figs13/diagram_capstone_workflow.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# ---- evidence to action loop
fig, ax = plt.subplots(figsize=(12.4, 4.4)); ax.axis("off"); ax.set_xlim(0, 12.4); ax.set_ylim(0, 4.4)
ax.text(0.15, 4.18, "From evidence to action, and back: the loop that closes with an evaluation", fontsize=13.5, fontweight="bold", color=INK, va="center")
lp = [("Question", "The decision needed", "#1B998B"), ("Design", "Sampling and its frame", "#3D7EA6"), ("Data", "Molecular and epidemiological", "#F3A712"), ("Analysis", "With uncertainty stated", "#E4572E"), ("Decision", "And its evaluation", "#2E294E")]
for i, (t, d, col) in enumerate(lp):
    x = 0.15 + i * 2.45; box(ax, x, 1.8, 2.2, 1.85); ax.add_patch(FancyBboxPatch((x, 3.55), 2.2, 0.1, boxstyle="round,pad=0.0,rounding_size=0.04", fc=col, ec=col)); ax.text(x + 0.15, 3.1, t, fontsize=11.5, fontweight="bold", color=INK, va="center"); ax.text(x + 0.15, 2.4, W(d, 22), fontsize=9.6, color=INK, va="center", linespacing=1.3)
    if i < 4: ax.annotate("", xy=(x + 2.43, 2.7), xytext=(x + 2.22, 2.7), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3))
ax.add_patch(FancyArrowPatch((11.2, 1.75), (1.2, 1.75), connectionstyle="arc3,rad=-0.2", arrowstyle="-|>", mutation_scale=16, color="#c53030", lw=1.8))
ax.text(6.2, 0.05, "Close the loop: evaluate what the action changed, and let the answer reshape the next question", fontsize=10, color="#c53030", ha="center", va="center", style="italic")
fig.savefig("figs13/diagram_evidence_to_action.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("diagrams written")
