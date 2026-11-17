"""Chapter 6 diagram: from sequences and dates to growth, five steps and the risk at each."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
INK, MUTED = "#1a202c", "#4a5568"
steps = [("1 Sequences and dates", "Alignment, tip dates,\nday-level where possible", "Wrong or missing dates", "Errors mimic rate variation"), ("2 Temporal signal", "Root-to-tip regression,\ndate randomisation, BETS", "No signal", "The prior drives every age"), ("3 Clock and calibration", "Strict, relaxed or local clock;\ntip dates, nodes, rates", "Wrong or circular anchor", "A false date biases every age"), ("4 Time tree", "TreeTime or BEAST 2;\nages with intervals", "Intervals dropped", "A point date over-claims"), ("5 Growth and R", "Skyline or birth-death;\ngeneration interval", "Sampling bias", "A skyline is not incidence")]
fig = plt.figure(figsize=(15.28, 7.4), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1528); ax.set_ylim(740, 0); ax.axis("off")
ax.text(24, 38, "From sequences and dates to growth: five steps, a risk at each", fontsize=16.5, fontweight="bold", color=INK, va="center")
w, h, gap = 270, 170, 31
for i, (t, d, r, rd) in enumerate(steps):
    x = 24 + i * (w + gap); col = "#0b3c49" if i % 2 == 0 else "#1b998b"
    ax.add_patch(FancyBboxPatch((x, 80), w, h, boxstyle="round,pad=0,rounding_size=12", fc=col, ec=col)); ax.text(x + 16, 112, t, fontsize=13.5, fontweight="bold", color="white", va="center"); ax.text(x + 16, 190, d, fontsize=11.6, color="#e6f1ef", va="center", linespacing=1.5)
    ax.plot([x + w / 2] * 2, [250, 300], color=MUTED, lw=1.5, ls=(0, (4, 3)))
    ax.add_patch(FancyBboxPatch((x, 300), w, 150, boxstyle="round,pad=0,rounding_size=12", fc="#ebf8ff", ec="#3d7ea6", lw=1.8, ls=(0, (4, 3)))); ax.text(x + 16, 340, r, fontsize=13, fontweight="bold", color=INK, va="center"); ax.text(x + 16, 395, rd, fontsize=11.6, color=MUTED, va="center")
    if i < 4: ax.annotate("", xy=(x + w + gap - 4, 165), xytext=(x + w + 4, 165), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6))
hands = ["Hands-on 1, 3", "Hands-on 1, 3", "Hands-on 3, 4", "Hands-on 3", "Hands-on 2, 3"]
for i, hh in enumerate(hands): ax.text(24 + i * (w + gap), 478, hh, fontsize=11, color=MUTED)
ax.add_patch(FancyBboxPatch((24, 520), 1480, 70, boxstyle="round,pad=0,rounding_size=10", fc="#edf2f7", ec="#cbd5e0")); ax.text(44, 555, "Each step depends on the one before; report the interval, the clock model and the calibration with every date.", fontsize=13, color=INK, va="center")
ax.add_patch(FancyBboxPatch((24, 620), 1480, 70, boxstyle="round,pad=0,rounding_size=10", fc="#fdf0d0", ec="#f3a712")); ax.text(44, 655, "No temporal signal, no reliable dates. A rate measured over years should not date events over centuries.", fontsize=13, fontweight="bold", color=INK, va="center")
fig.savefig("diagram_ch6_pipeline.png", dpi=100, facecolor="white"); print("ok")
