"""
Chapter 8 - tables and figures of the two new stress tests.
Reads dating_reruns.csv (ex8_5_dating_reruns.py), timetree_summary.json (the dating that the chapter's other figures use) and sampling_truth_runs.csv (ex8_6_sampling_truth.py).
Writes table_dating_reruns_summary.csv, table_sampling_truth_summary.csv and figs8/fig5_dating_reruns.png, figs8/fig9_sampling_truth.png.
"""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs8", exist_ok=True)
NAMES = ["A", "B", "C", "D"]
# ------------------------------------------------------------- dating reruns
d = pd.read_csv("dating_reruns.csv"); s = json.load(open("timetree_summary.json")); ev = pd.read_csv("events.csv")
import phylo8
is_am = lambda x: x in phylo8.AMERICAS or x == "Other Americas"
ev["intro"] = [(not is_am(a)) and is_am(b) for a, b in zip(ev.src, ev.dst)]
main = ev[ev.intro].sort_values("child_tips", ascending=False).iloc[0]
ship = dict(clock_rate=s["clock_rate"], root=s["root_date"], brazil_mrca=s["brazil_mrca"], brazil_parent=s["brazil_parent"], n_events=len(ev), n_intro=int(ev.intro.sum()), main_lo=float(main.date_lo), main_hi=float(main.date_hi), main_tips=int(main.child_tips))
rows = []
for c, lab in [("clock_rate", "Clock rate (substitutions per site per year)"), ("root", "Date of the root"), ("brazil_parent", "Parent of the Brazilian common ancestor"), ("brazil_mrca", "Common ancestor of the Brazilian genomes"), ("n_events", "Inferred migration events"), ("n_intro", "Introductions into the Americas"), ("main_tips", "Genomes below the main introduction")]:
    rows.append(dict(measure=lab, key=c, shipped=ship[c], min=d[c].min(), median=d[c].median(), max=d[c].max(), range=d[c].max() - d[c].min()))
S = pd.DataFrame(rows); S.round(6).to_csv("table_dating_reruns_summary.csv", index=False); print(S.round(5).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.5, 1]}); order = list(range(len(d)))
ax[0].axvspan(2013 + 4 / 12, 2013 + 11 / 12, color="#f3a712", alpha=0.25, label="Published estimate, May to December 2013")
for k, r in enumerate(d.itertuples()): ax[0].plot([r.main_lo, r.main_hi], [k, k], color="#1b998b", lw=4, solid_capstyle="butt")
ax[0].plot([ship["main_lo"], ship["main_hi"]], [len(d), len(d)], color="#e4572e", lw=4, solid_capstyle="butt", label="Run used in Hands-on 1 to 3")
ax[0].set_yticks(list(range(len(d))) + [len(d)]); ax[0].set_yticklabels([f"Rerun {k}" for k in d.rerun] + ["Chapter run"]); ax[0].set_xlim(2011.7, 2014.2); ax[0].set_xlabel("Date of the branch that brought the virus into Brazil"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper left"); ax[0].set_title("Same data, ten reruns of the dating step")
ax[1].ticklabel_format(useOffset=False, style="plain", axis="y"); ax[1].scatter(d.clock_rate * 1e4, d.brazil_mrca, color="#1b998b", s=40, label="Reruns"); ax[1].scatter([ship["clock_rate"] * 1e4], [ship["brazil_mrca"]], color="#e4572e", s=60, marker="D", label="Chapter run")
ax[1].set_xlabel("Estimated clock rate (10$^{-4}$ per site per year)"); ax[1].set_ylabel("Date of the Brazilian common ancestor"); ax[1].legend(frameon=False, fontsize=8.5); ax[1].set_title("Faster clock, later introduction")
fig.tight_layout(); fig.savefig("figs8/fig5_dating_reruns.png", dpi=200); plt.close(fig)
# ------------------------------------------------------------- sampling strategies against the truth
t = pd.read_csv("sampling_truth_runs.csv"); rows = []; per = {}
for sc in t.scenario.unique():
    for st, g in t[t.scenario == sc].groupby("strategy", sort=False):
        tr = g[[f"true_out_{x}" for x in NAMES]].to_numpy(float); inf = g[[f"inf_out_{x}" for x in NAMES]].to_numpy(float)
        share_t = tr / np.maximum(tr.sum(1, keepdims=True), 1); share_i = inf / np.maximum(inf.sum(1, keepdims=True), 1); bias = share_i[:, 3] - share_t[:, 3]; n = len(g)
        per[(sc, st)] = dict(bias=bias, root=(g.true_root == g.inf_root).to_numpy(float), top=(inf.argmax(1) == tr.argmax(1)).astype(float))
        rows.append(dict(scenario=sc, strategy=st, epidemics=n, tips=g.n_tips.mean(), share_tips_D=float(np.mean(g.tips_D / g.n_tips)), true_events=g.true_events.mean(), inferred_events=g.inf_events.mean(),
            true_share_D=share_t[:, 3].mean(), inferred_share_D=share_i[:, 3].mean(), bias_D=bias.mean(), bias_D_lo=bias.mean() - 1.96 * bias.std(ddof=1) / np.sqrt(n), bias_D_hi=bias.mean() + 1.96 * bias.std(ddof=1) / np.sqrt(n),
            root_correct=per[(sc, st)]["root"].mean(), top_exporter_correct=per[(sc, st)]["top"].mean()))
T = pd.DataFrame(rows); T.round(6).to_csv("table_sampling_truth_summary.csv", index=False); print(T.round(3).to_string(index=False))
labels = {"as_sequenced": "As sequenced", "as_sequenced_sbc": "As sequenced,\nbias correction 2", "cap_per_place": "Capped at 15\nper place", "proportional": "Proportional\nto infections"}; order = list(labels)
fig, ax = plt.subplots(1, 3, figsize=(13, 4.1)); cols = {"hub": "#e4572e", "no_hub": "#1b998b"}
for j, sc in enumerate(["hub", "no_hub"]):
    sub = T[T.scenario == sc].set_index("strategy").loc[order]; x = np.arange(4) + (j - 0.5) * 0.18
    ax[0].errorbar(x, sub.bias_D, yerr=[sub.bias_D - sub.bias_D_lo, sub.bias_D_hi - sub.bias_D], fmt="o", color=cols[sc], capsize=3, label=f"{sc.replace('_', '-')} scenario")
    ax[1].bar(np.arange(4) + (j - 0.5) * 0.36, sub.root_correct, width=0.34, color=cols[sc], label=f"{sc.replace('_', '-')} scenario"); ax[2].bar(np.arange(4) + (j - 0.5) * 0.36, sub.top_exporter_correct, width=0.34, color=cols[sc])
ax[0].axhline(0, color="#718096", ls=":"); ax[0].set_ylabel("Inferred minus true share of all exports from place D"); ax[0].set_title("Over-sampled place D looks like an exporter")
ax[1].set_ylabel("Share of epidemics with the root place correct"); ax[1].set_title("Root place"); ax[2].set_ylabel("Share with the top exporter correct"); ax[2].set_title("Top exporter")
for a in ax: a.set_xticks(range(4)); a.set_xticklabels([labels[o] for o in order], fontsize=8)
ax[1].set_ylim(0, 1.05); ax[2].set_ylim(0, 1.05); ax[0].legend(frameon=False, fontsize=8.5); fig.tight_layout(); fig.savefig("figs8/fig9_sampling_truth.png", dpi=200); print("figures written")
