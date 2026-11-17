"""Chapter 8, Exercise 2 (continued) - Figure: how the date of the introduction into Brazil depends on the clock rate.
Input: rate_sensitivity.csv from ex8_1_dated_tree.py."""
import pandas as pd, numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
d = pd.read_csv("rate_sensitivity.csv"); s = json.load(open("timetree_summary.json"))
fig, ax = plt.subplots(figsize=(8.2, 3.9))
ax.axhspan(2013 + 4 / 12, 2013 + 11 / 12, color="#fbd38d", alpha=0.55, lw=0)
ax.text(4.75e-4, 2013.63, "published estimate\n(Faria et al. 2016):\nMay to December 2013", fontsize=8.5, color="#7b341e", va="center", ha="left")
for r in d.itertuples():
    c = "#2b6cb0" if r.setting == "estimated" else "#4a5568"
    ax.plot([r.clock_rate] * 2, [r.brazil_parent, r.brazil_mrca], lw=7, solid_capstyle="butt", color=c)
    ax.text(r.clock_rate, r.brazil_mrca + 0.09, f"{r.brazil_parent:.1f} to {r.brazil_mrca:.1f}", ha="center", fontsize=8.5)
ax.set_xticks(d.clock_rate); ax.set_xticklabels([f"{x * 1e4:.1f}\n(estimated)" if st == "estimated" else f"{x * 1e4:.1f}\n(fixed)" for x, st in zip(d.clock_rate, d.setting)])
ax.set(xlabel="Clock rate, $10^{-4}$ substitutions per site per year", ylabel="Date of the introduction branch", ylim=(2011.8, 2014.5), xlim=(4.6e-4, 1.3e-3),
       title="The date of the introduction into Brazil moves with the assumed clock rate")
fig.tight_layout(); fig.savefig("figs8/fig4_clock_sensitivity.png", dpi=200); plt.close(fig); print("ok")
