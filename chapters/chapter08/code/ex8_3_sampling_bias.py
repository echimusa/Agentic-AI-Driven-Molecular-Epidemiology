"""
Chapter 8, Exercise 3 - Does uneven sampling change the answer?
Cap the number of genomes per country at 10, 20 or 40, redo the ancestral-country reconstruction on the pruned tree 25 times
for each cap, and compare with the full, uneven sample. Results are cached per cap, so the script can be rerun to finish.
Needs the outputs of ex8_1_dated_tree.py and phylo8.py.
"""
import os, copy, time
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phylo8 import *

REPS, CAPS, SEED = 25, (10, 20, 40), 3
meta = load_states(); tree0, root_date = load_tree()
traits_all = {t.name: meta.loc[t.name, "state"] for t in tree0.get_terminals()}
FILE = "bias_replicates.csv"
rows = pd.read_csv(FILE) if os.path.exists(FILE) else pd.DataFrame()

def one_run(tip_names, label, cap, rep):
    tr = copy.deepcopy(tree0)
    for t in [t.name for t in tr.get_terminals() if t.name not in tip_names]: tr.prune(t)
    traits = {n: traits_all[n] for n in tip_names}
    mug, st, post = run_mugration(tr, traits, iterations=3, seed=1)
    dd = clade_dates(mug.tree, root_date); ev = events(mug, st, dd)
    ev["intro"] = [(not is_americas(a)) and is_americas(b) for a, b in zip(ev.src, ev.dst)]
    main = ev[ev.intro].sort_values("child_tips", ascending=False).iloc[0]
    tips = pd.Series(traits).value_counts()
    out = []
    for s in tips.index:
        out.append(dict(cap=cap, rep=rep, state=s, tips=int(tips[s]), out_events=int((ev.src == s).sum()), in_events=int((ev.dst == s).sum()),
                        n_events=len(ev), n_intro=int(ev.intro.sum()), main_src=main.src, main_dst=main.dst, main_date=float(main.date), main_tips=int(main.child_tips)))
    return out

t0 = time.time()
if "full" not in set(rows.get("cap", pd.Series(dtype=str)).astype(str)):
    rows = pd.concat([rows, pd.DataFrame(one_run(set(traits_all), "full", "full", 0))], ignore_index=True); rows.to_csv(FILE, index=False)
for cap in CAPS:
    if str(cap) in set(rows["cap"].astype(str)): continue
    rng = np.random.default_rng(SEED + cap); new = []
    for rep in range(REPS):
        keep = set()
        for s, names in pd.Series(traits_all).groupby(pd.Series(traits_all)):
            names = list(names.index); keep |= set(rng.choice(names, size=min(cap, len(names)), replace=False))
        new += one_run(keep, cap, cap, rep)
    rows = pd.concat([rows, pd.DataFrame(new)], ignore_index=True); rows.to_csv(FILE, index=False)
    print(f"cap {cap}: done ({time.time() - t0:.0f} s)")
    if time.time() - t0 > 200: break
done = set(rows["cap"].astype(str)); print("completed:", sorted(done))
if len(done) < len(CAPS) + 1: raise SystemExit("rerun to finish the remaining caps")

# ---------------------------------------------------------------- summaries
am = lambda df: df[df.state.map(is_americas) & ~df.state.str.startswith("Other")]
full = rows[rows.cap.astype(str) == "full"]; fr = full.set_index("state")
print("\nfull, uneven sample: outgoing events", fr.out_events.sort_values(ascending=False).head(4).to_dict())
def top_share(df):
    """Share of runs in which each American country has the most outgoing events; ties split the credit equally."""
    out = {}
    for _, g in am(df).groupby("rep"):
        tops = g[g.out_events == g.out_events.max()].state.tolist()
        for t_ in tops: out[t_] = out.get(t_, 0) + 1 / len(tops)
    n = am(df).rep.nunique(); return pd.Series({k: v / n for k, v in out.items()})
shares = {"full": top_share(full)}
for cap in CAPS:
    r = rows[rows.cap.astype(str) == str(cap)]
    shares[str(cap)] = top_share(r)
    mainreps = r.drop_duplicates("rep")
    med = {s_: float(r[r.state == s_].out_events.median()) for s_ in ("Nicaragua", "Brazil", "Dominican Republic")}
    print(f"cap {cap}: top exporter among American countries {shares[str(cap)].round(2).to_dict()} | median outgoing events {med} | "
          f"main introduction is Asia-Pacific to Brazil in {(mainreps.main_dst == 'Brazil').mean():.0%} of replicates (median date {mainreps.main_date.median():.2f}) | "
          f"introductions into the Americas per replicate: median {mainreps.n_intro.median():.0f}")
allr = rows[rows.state.map(lambda s_: not s_.startswith("Other"))]
rho, _ = stats.spearmanr(allr.tips, allr.out_events); print(f"across all runs and countries: outgoing events vs genomes sampled, Spearman rho = {rho:.2f} (n = {len(allr)})")
print("full-sample events per genome:", (fr.out_events / fr.tips).sort_values(ascending=False).head(4).round(2).to_dict())

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
scen = ["full", "10", "20", "40"]; labels = ["full,\nuneven", "cap 10", "cap 20", "cap 40"]
names = ["Nicaragua", "Brazil", "Dominican Republic", "Honduras", "Colombia"]
cols = {"Nicaragua": "#9ecae1", "Brazil": "#1f77b4", "Dominican Republic": "#b794d6", "Honduras": "#9467bd", "Colombia": "#ffbb78"}
bottom = np.zeros(len(scen))
for n in names:
    v = np.array([shares[sc].get(n, 0.0) for sc in scen]); ax[0].bar(range(len(scen)), v, bottom=bottom, color=cols[n], label=n, edgecolor="white"); bottom += v
ax[0].set_xticks(range(len(scen))); ax[0].set_xticklabels(labels); ax[0].set(ylabel="Share of runs in which the country\nhas the most outgoing events", title="Which country is the top exporter?", ylim=(0, 1.02))
ax[0].yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(1.0)); ax[0].legend(frameon=False, fontsize=8.3, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
focus = ["Nicaragua", "Brazil", "Dominican Republic", "Honduras", "Colombia"]; pos = np.arange(len(focus)); w = 0.2
pal = {"10": "#2b6cb0", "20": "#38a169", "40": "#d69e2e"}
for k, c in enumerate(["10", "20", "40"]):
    r = rows[rows.cap.astype(str) == c]
    bp = ax[1].boxplot([r[r.state == s_].out_events.values for s_ in focus], positions=pos + (k - 1) * w, widths=w * 0.9, patch_artist=True, showfliers=False)
    for b in bp["boxes"]: b.set(facecolor=pal[c], alpha=0.65, edgecolor="black")
    for m in bp["medians"]: m.set(color="black")
ax[1].scatter(pos + 2 * w, [fr.loc[s_, "out_events"] for s_ in focus], marker="D", color="#c53030", zorder=4)
ax[1].set_xticks(pos + w / 2); ax[1].set_xticklabels(focus, fontsize=8.5, rotation=15)
ax[1].set(ylabel="Inferred outgoing events", title="Outgoing events: subsamples vs full sample")
ax[1].legend(handles=[plt.Line2D([], [], marker="s", ls="", color=pal[c], alpha=0.7, label=f"cap {c}") for c in pal] + [plt.Line2D([], [], marker="D", ls="", color="#c53030", label="full sample")], frameon=False, fontsize=8.3)
os.makedirs("figs8", exist_ok=True); fig.tight_layout(); fig.savefig("figs8/fig3_sampling_bias.png", dpi=200); plt.close(fig)
print("figure written to figs8/fig3_sampling_bias.png")
