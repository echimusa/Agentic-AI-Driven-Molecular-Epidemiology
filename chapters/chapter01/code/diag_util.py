"""Chapter 3 diagrams: the quality-gate pipeline and the agent-and-human workflow, drawn with matplotlib."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
INK, QUIET, EDGE, ACC = "#1a202c", "#4a5568", "#718096", "#2b6cb0"
def box(ax, x, y, w, h, title, lines=(), fill=None, dashed=False, accent=False, tsize=9.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=1.2", fc=(ACC if accent else (fill or "none")), alpha=0.14 if accent else 1, ec=ACC if accent else EDGE, lw=2 if accent else 1.2, ls="--" if dashed else "-"))
    ax.text(x + 1.4, y + h - 3.0, title, fontsize=tsize, fontweight="bold", color=INK, va="top")
    top = 6.3 if title else 2.3
    for i, l in enumerate(lines): ax.text(x + 1.4, y + h - top - 3.0 * i, l, fontsize=8.0, color=QUIET if not accent else INK, va="top")
def arrow(ax, x1, y1, x2, y2): ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, color=EDGE, lw=1.2))
