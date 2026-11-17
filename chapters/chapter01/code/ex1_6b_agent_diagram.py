"""Chapter 1: the agentic AI foundations diagram (loop, components and levels of autonomy)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from diag_util import box, arrow, INK, QUIET, EDGE, ACC
fig, ax = plt.subplots(figsize=(11.6, 7.4)); ax.set_xlim(0, 100); ax.set_ylim(0, 66); ax.axis("off")
ax.text(1, 63.2, "An agent is a model in a loop with tools, memory and guardrails", fontsize=13, fontweight="bold", color=INK, va="center")
names = [("Goal", ["What to find out,", "and the limits"]), ("Plan", ["Choose the", "next step"]), ("Act", ["Call a tool:", "code, data, search"]), ("Observe", ["Read the", "result"]), ("Check", ["Is it done, wrong", "or out of budget?"])]
w, g = 10.6, 2.3; xs = [1 + i * (w + g) for i in range(5)]
for i, (t, l) in enumerate(names): box(ax, xs[i], 40, w, 14, t, l, accent=(i == 2), tsize=10)
for i in range(4): arrow(ax, xs[i] + w, 47, xs[i + 1], 47)
ax.add_patch(FancyArrowPatch((xs[4] + w / 2, 54), (xs[1] + w / 2, 54), connectionstyle="arc3,rad=0.42", arrowstyle="-|>", mutation_scale=10, color=EDGE, lw=1.3)); ax.text((xs[1] + xs[4]) / 2 + 5, 57.0, "not finished: update the state and plan again", fontsize=8.6, color=QUIET, ha="center", bbox=dict(fc="white", ec="none", pad=1.5))
box(ax, xs[3] - 2, 21, 25.6, 10.2, "Report to a person", ["Finished, stuck or uncertain"], fill="#fffaf0", tsize=9.4); arrow(ax, xs[4] + w / 2, 40, xs[4] + w / 2, 31.2)
ax.text(1.2, 33, "The loop repeats until the check says stop.", fontsize=8.8, color=QUIET)
cx, cw = 71, 28; comp = [("Model", "Reasons about the next step"), ("Tools", "Code, databases, pipelines, search"), ("Memory", "State, logs, earlier results"), ("Guardrails", "Step budget, thresholds, human gate")]
ax.text(cx, 56.2, "Parts of an agent", fontsize=9.6, fontweight="bold", color=INK)
for k, (t, d) in enumerate(comp): box(ax, cx, 44 - 9.4 * k, cw, 8.8, t, [d], fill="#ebf8ff", tsize=9.6)
ax.text(1, 14.6, "Levels of autonomy", fontsize=9.6, fontweight="bold", color=INK)
lev = [("1 Assistant", ["Answers a question", "e.g. explains a SNP table"]), ("2 Workflow", ["Fixed steps in code", "e.g. a nightly cluster run"]), ("3 Agent", ["The model picks the steps", "e.g. triage of new genomes"]), ("4 Multi-agent", ["Roles that check each other", "e.g. QC, analysis, review"])]
lw, lg = 22.6, 2.2
for i, (t, l) in enumerate(lev):
    x = 1 + i * (lw + lg); box(ax, x, 1.0, lw, 11.6, t, l, fill="#f7fafc", accent=(i == 2), tsize=9.6)
    if i < 3: arrow(ax, x + lw, 6.8, x + lw + lg, 6.8)
fig.savefig("figs/fig_agent_foundations.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("agent diagram written")
