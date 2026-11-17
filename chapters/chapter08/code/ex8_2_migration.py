"""
Chapter 8, Exercise 2 - Ancestral countries, migration events and the figures that show them.
Input: the time-scaled tree and metadata from ex8_1_dated_tree.py. Output: events.csv, state_summary.csv and four figures.
Needs: pip install biopython treetime pandas numpy matplotlib shapely
"""
import os, json, urllib.request, collections
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, FancyArrowPatch
from shapely.geometry import shape
from phylo8 import *

os.makedirs("figs8", exist_ok=True)
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})

# ---------------------------------------------------------------- 1. reconstruct ancestral countries
meta = load_states(); tree, root_date = load_tree()
traits = {t.name: meta.loc[t.name, "state"] for t in tree.get_terminals()}
mug, st, post = run_mugration(tree, traits, iterations=5, seed=1)           # discrete trait on a fixed dated tree
dd = clade_dates(mug.tree, root_date)
ev = events(mug, st, dd); ev.to_csv("events.csv", index=False)
ev["intro"] = [(not is_americas(a)) and is_americas(b) for a, b in zip(ev.src, ev.dst)]
tips_per_state = pd.Series(traits).value_counts()
print(f"{len(traits)} tips, {tips_per_state.size} states, {len(ev)} inferred migration events, {int(ev.intro.sum())} introductions into the Americas")
root = mug.tree.root; print(f"root state {st[id(root)]} (posterior {post[id(root)]:.2f}), root date {root_date:.1f}")
print("introductions into the Americas (interval = parent date to child date):")
print(ev[ev.intro].sort_values("date_lo")[["src", "dst", "date_lo", "date_hi", "child_tips"]].round(2).to_string(index=False))
out = ev.groupby("src").size(); inn = ev.groupby("dst").size()
summary = pd.DataFrame({"tips": tips_per_state, "events_out": out, "events_in": inn}).fillna(0).astype(int).sort_values("tips", ascending=False)
summary.to_csv("state_summary.csv"); print(summary.to_string())

# ---------------------------------------------------------------- colours and basemap
region = lambda s: "Americas" if is_americas(s) else ("Asia-Pacific" if "Asia" in s or s in ASIA_PAC else "Europe/Africa")
states = sorted(tips_per_state.index, key=lambda s: (region(s), -tips_per_state[s]))
cmap = plt.get_cmap("tab20"); COL = {s: cmap(i % 20) for i, s in enumerate(states)}

NE = "ne_50m_countries.geojson"
if not os.path.exists(NE):
    urllib.request.urlretrieve("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson", NE)
geo = json.load(open(NE))["features"]
NAMES = {"USA": "United States of America", "Dominican Republic": "Dominican Rep.", "French Polynesia": "Fr. Polynesia", "Viet Nam": "Vietnam", "Virgin Islands": "U.S. Virgin Is."}
polys = {f["properties"]["NAME"]: shape(f["geometry"]) for f in geo}
def centroid(country):
    g = polys.get(NAMES.get(country, country))
    if g is None: return None
    big = max(g.geoms, key=lambda p: p.area) if g.geom_type == "MultiPolygon" else g       # largest piece (e.g. the contiguous USA)
    return np.array([big.centroid.x, big.centroid.y])
xy, is_group = {}, {}
for s_ in states:
    members = meta[meta.state == s_].cty.value_counts()
    pts = [(centroid(c), n) for c, n in members.items() if centroid(c) is not None]
    lons = np.array([p[0] for p, _ in pts]); lats = np.array([p[1] for p, _ in pts]); w = np.array([n for _, n in pts], float)
    if region(s_) == "Asia-Pacific": lons = lons % 360                                      # French Polynesia sits across the antimeridian
    lon = np.average(lons, weights=w); lon = lon - 360 if lon > 180 else lon
    xy[s_] = np.array([lon, np.average(lats, weights=w)]); is_group[s_] = s_.startswith("Other")

# ---------------------------------------------------------------- Figure 1: dated tree, branches coloured by inferred country
order = {}
def lay(c, y=[0]):
    if c.is_terminal(): order[id(c)] = y[0]; y[0] += 1
    else:
        for ch in c.clades: lay(ch)
        order[id(c)] = np.mean([order[id(ch)] for ch in c.clades])
lay(mug.tree.root)
segs, cols = [], []
for c in mug.tree.find_clades(order="level"):
    for ch in c.clades:
        segs.append([(dd[id(c)], order[id(ch)]), (dd[id(ch)], order[id(ch)])]); cols.append(COL[st[id(ch)]])
    if c.clades:
        ys = [order[id(ch)] for ch in c.clades]
        segs.append([(dd[id(c)], min(ys)), (dd[id(c)], max(ys))]); cols.append(COL[st[id(c)]])
fig, ax = plt.subplots(1, 2, figsize=(11.5, 8), gridspec_kw={"width_ratios": [3.4, 1]})
ax[0].add_collection(LineCollection(segs, colors=cols, linewidths=1.0))
ax[0].set_xlim(2005, 2019.4); ax[0].set_ylim(-5, len(traits) + 5); ax[0].set_yticks([])
ax[0].set(xlabel="Year", title=f"Dated tree of {len(traits)} Zika genomes, branches coloured by inferred country"); ax[0].spines["left"].set_visible(False)
ax[0].text(2005.1, len(traits) + 1, f"root at {root_date:.1f}, off the left edge", fontsize=8.5, color="#4a5568")
ax[1].barh(range(len(states)), [tips_per_state[s] for s in states], color=[COL[s] for s in states])
ax[1].set_yticks(range(len(states))); ax[1].set_yticklabels(states, fontsize=8.5); ax[1].invert_yaxis()
ax[1].set(xlabel="Genomes sampled", title="Sampling is uneven")
fig.tight_layout(); fig.savefig("figs8/fig1_dated_tree.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 2: migration matrix
M = pd.crosstab(ev.src, ev.dst).reindex(index=states, columns=states, fill_value=0)
fig, ax = plt.subplots(figsize=(8.6, 7.4))
im = ax.imshow(M.values, cmap="Blues", vmin=0, vmax=max(M.values.max(), 1))
ax.set_xticks(range(len(states))); ax.set_xticklabels(states, rotation=60, ha="right", fontsize=8.5)
ax.set_yticks(range(len(states))); ax.set_yticklabels(states, fontsize=8.5)
for i in range(len(states)):
    for j in range(len(states)):
        if M.values[i, j]: ax.text(j, i, int(M.values[i, j]), ha="center", va="center", fontsize=8, color="white" if M.values[i, j] > M.values.max() * 0.55 else "black")
cuts = [i - 0.5 for i in range(1, len(states)) if region(states[i]) != region(states[i - 1])]
for c_ in cuts: ax.axhline(c_, color="#2d3748", lw=1.2); ax.axvline(c_, color="#2d3748", lw=1.2)
ax.set(xlabel="Destination (child lineage)", ylabel="Source (parent lineage)", title="Inferred migration events between countries (regions separated by lines)")
fig.colorbar(im, ax=ax, fraction=0.04, label="Events"); fig.tight_layout(); fig.savefig("figs8/fig2_migration_matrix.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 3: map with arcs
fig, ax = plt.subplots(figsize=(11.5, 5.0))
for f in geo:
    g = shape(f["geometry"]); parts = g.geoms if g.geom_type == "MultiPolygon" else [g]
    for p in parts: ax.add_patch(Polygon(np.array(p.exterior.coords), closed=True, facecolor="#edf0f3", edgecolor="white", linewidth=0.3))
flows = ev.groupby(["src", "dst"]).agg(n=("date", "size"), tips=("child_tips", "sum")).reset_index()
for r in flows[flows.n >= 3].itertuples():
    ax.add_patch(FancyArrowPatch(xy[r.src], xy[r.dst], connectionstyle="arc3,rad=0.28", arrowstyle="-|>", mutation_scale=9, lw=0.8 + 0.45 * r.n, color=COL[r.src], alpha=0.8, zorder=3))
main = ev[ev.intro].sort_values("child_tips", ascending=False).iloc[0]                    # the introduction that seeded the American epidemic
ax.add_patch(FancyArrowPatch(xy[main.src], xy[main.dst], connectionstyle="arc3,rad=-0.22", arrowstyle="-|>", mutation_scale=16, lw=3.2, color="#c53030", zorder=4))
for s_ in states:
    ax.scatter(*xy[s_], s=22 + 5 * tips_per_state[s_], marker="s" if is_group[s_] else "o", facecolor="white" if is_group[s_] else COL[s_], edgecolor="black", linewidth=1.0 if is_group[s_] else 0.6, zorder=5)
for s_ in states:
    if tips_per_state[s_] >= 29 or s_ in ("Thailand", "Cuba", "Dominican Republic"):
        off = {"Dominican Republic": (-6, 14), "Puerto Rico": (10, -15), "Cuba": (-30, 9), "Nicaragua": (-48, -16)}.get(s_, (5, 6))
        ax.annotate(s_, xy[s_], xytext=off, textcoords="offset points", fontsize=8.5, zorder=6)
ax.annotate(f"{main.src} to {main.dst}: {main.date_lo:.1f} to {main.date_hi:.1f}", xy=(10, -30), fontsize=9, color="#c53030", ha="center")
ax.set_xlim(-135, 178); ax.set_ylim(-42, 52); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Inferred spread of Zika virus. Red: the introduction behind 487 of 591 genomes. Arcs: flows with 3 or more events. Squares: groups of countries.", fontsize=9)
fig.tight_layout(pad=0.3); fig.savefig("figs8/fig5_map.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 4: timeline of introductions and onward spread
FIRST_DETECTION, PUB_LO, PUB_HI = 2015 + 4 / 12, 2013 + 4 / 12, 2013 + 11 / 12          # May 2015; May-Dec 2013 (Faria et al. 2016)
fig, ax = plt.subplots(2, 1, figsize=(8.8, 6.6), gridspec_kw={"height_ratios": [1, 1.2]}, sharex=True)
intro = ev[ev.intro].sort_values("date_lo")
ax[0].axvspan(PUB_LO, PUB_HI, color="#fbd38d", alpha=0.55, lw=0)
ax[0].text((PUB_LO + PUB_HI) / 2, -0.6, "published estimate\n(Faria et al. 2016)", ha="center", va="bottom", fontsize=8, color="#7b341e")
for k, r in enumerate(intro.itertuples()):
    ax[0].plot([r.date_lo, r.date_hi], [k, k], lw=6, solid_capstyle="butt", color=COL[r.dst]); ax[0].scatter([r.date], [k], color="black", s=14, zorder=3)
ax[0].set_yticks(range(len(intro))); ax[0].set_yticklabels([f"{r.src} to {r.dst}\n({r.child_tips} genome{'s' if r.child_tips != 1 else ''})" for r in intro.itertuples()], fontsize=8.3)
ax[0].set_ylim(len(intro) - 0.4, -1.3); ax[0].set_title("Introductions into the Americas: the interval on each branch")
within = ev[(~ev.intro) & ev.src.map(is_americas) & ev.dst.map(is_americas)].sort_values("date")
ax[1].step(within.date, np.arange(1, len(within) + 1), where="post", color="#2b6cb0", lw=2)
for a_ in ax: a_.axvline(FIRST_DETECTION, color="#c53030", ls="--", lw=1.2)
ax[1].text(FIRST_DETECTION + 0.05, 3, "first confirmed in Brazil,\nMay 2015", color="#c53030", fontsize=8.3)
ax[1].set(xlabel="Year", ylabel="Cumulative events", title="Onward spread between American countries"); ax[1].set_xlim(2011.5, 2017.2)
fig.tight_layout(); fig.savefig("figs8/fig6_timeline.png", dpi=200); plt.close(fig)
print("figures written to figs8/")
