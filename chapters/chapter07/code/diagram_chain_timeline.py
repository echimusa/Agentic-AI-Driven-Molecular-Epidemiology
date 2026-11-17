"""Chapter 7 diagram: the timeline of two linked cases and the intervals that can be measured."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
INK, MUTED, DK, TEAL, ORG, AMB, BLU = "#1a202c", "#4a5568", "#0b3c49", "#1b998b", "#e4572e", "#f3a712", "#3d7ea6"
fig = plt.figure(figsize=(15.28, 7.2), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1528); ax.set_ylim(720, 0); ax.axis("off")
ax.text(24, 36, "Two linked cases: what happens, and what can be measured", fontsize=16.5, fontweight="bold", color=INK, va="center")
def bar(x0, x1, y, col, label, tc="white"):
    ax.add_patch(Rectangle((x0, y), x1 - x0, 44, fc=col, ec=col)); ax.text((x0 + x1) / 2, y + 22, label, fontsize=11.5, color=tc, ha="center", va="center", fontweight="bold")
# infector (A)
ax.text(24, 150, "Infector A", fontsize=13, fontweight="bold", color=INK, va="center"); ax.text(24, 372, "Infectee B", fontsize=13, fontweight="bold", color=INK, va="center")
bar(210, 410, 128, MUTED, "latent"); bar(410, 760, 128, TEAL, "infectious"); ax.plot([210, 410, 760], [190, 190, 190], alpha=0)
bar(210, 600, 188, BLU, "incubation (infection to onset)"); ax.plot([600, 600], [120, 250], color=ORG, lw=2.2, ls=(0, (4, 3))); ax.text(600, 112, "onset A", fontsize=11.5, color=ORG, ha="center", fontweight="bold")
ax.text(210, 120, "infection A", fontsize=11.5, color=INK, ha="center")
# infectee (B) infected during A's infectious period
bar(560, 760, 350, MUTED, "latent"); bar(760, 1130, 350, TEAL, "infectious"); bar(560, 960, 410, BLU, "incubation (infection to onset)")
ax.plot([960, 960], [340, 470], color=ORG, lw=2.2, ls=(0, (4, 3))); ax.text(960, 332, "onset B", fontsize=11.5, color=ORG, ha="center", fontweight="bold"); ax.text(575, 322, "infection B (transmission A to B)", fontsize=11.5, color=INK, ha="left")
ax.plot([560, 560], [270, 350], color=INK, lw=1.6, ls=(0, (2, 2)))
# interval arrows
def dim(x0, x1, y, col, text):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color=col, lw=2.2)); ax.text((x0 + x1) / 2, y + 22, text, fontsize=11.8, color=col, ha="center", va="center", fontweight="bold")
dim(210, 560, 285, DK, "generation time (never observed)"); dim(600, 960, 505, ORG, "serial interval (onset to onset)")
ax.text(1180, 150, "Incubation period: infection to onset.\nLatent period: infection to infectious.\nInfectious periods can start before onset,\nso a serial interval can be short or negative.", fontsize=12, color=MUTED, va="top", linespacing=1.6)
ax.add_patch(FancyBboxPatch((24, 580), 1480, 100, boxstyle="round,pad=0,rounding_size=10", fc="#fdf0d0", ec=AMB)); ax.text(44, 630, "What a line list gives you is onset dates, so the serial interval stands in for the generation time (the two have similar means and different spreads).", fontsize=13, fontweight="bold", color=INK, va="center")
fig.savefig("diagram_chain_timeline.png", dpi=100, facecolor="white"); print("ok")
