"""
Chapter 1 schematic figures, redrawn from the ideas of the CDC genomic epidemiology toolkit (module 1.1) as original drawings:
  fig_epidemic_to_tree.png  - epidemic process, sampling, mutations, phylogeny and nested relatedness, from a small simulated outbreak (epi_sim.py; mutation rate raised for legibility)
  fig_investigation_uses.png - what sequencing adds to an investigation: hotspots, multiple introductions, a supplement to contact tracing
  fig_sequencing_workflow.png - from swab to shared sequence, and the two platform families
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform
from epi_sim import *
from diag_util import box, arrow, INK, QUIET, EDGE, ACC
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
# ---------------------------------------------------------------- figure A: epidemic to tree
for seed in range(2000):
    rng = np.random.default_rng(seed); o = simulate_outbreak(rng, R=1.7, k=None, Tg=12.0, shape=4.0, lam=0.22, n_cap=17, min_cases=15, max_tries=1)
    if o["n"] < 15: continue
    n = o["n"]; kids = {i: [] for i in range(n)}
    for i in range(1, n): kids[o["parent"][i]].append(i)
    dep = depths(o["parent"]); cand = [i for i in range(n) if dep[i] >= 2]
    if len(cand) < 6: continue
    samp = sorted(rng.choice(cand, 6, replace=False), key=lambda i: o["t_inf"][i]); G = [sample_genome(o, i, 6.0) for i in samp]
    D = np.array([[len(a ^ b) for b in G] for a in G]); iu = np.triu_indices(6, 1)
    if len(set(D[iu])) >= 4 and D[iu].max() <= 8 and D[iu].min() >= 1 and len({tuple(sorted(g)) for g in G}) == 6: break
print("toy outbreak seed", seed, "cases", n, "sampled", samp)
order = []; 
def dfs(i):
    order.append(i)
    for c in sorted(kids[i], key=lambda c: o["t_inf"][c]): dfs(c)
dfs(0); ypos = {c: len(order) - k for k, c in enumerate(order)}; letters = {c: "ABCDEF"[k] for k, c in enumerate(samp)}
def draw_network(ax, sampled=False, stars=False):
    for i in range(n):
        t0 = o["t_inf"][i]; ch = [o["t_inf"][c] for c in kids[i]]; t1 = max(ch + [t0 + 7]) + 2; faded = sampled and i not in letters
        ax.plot([t0, t1], [ypos[i]] * 2, color="#a0aec0" if faded else "#2d3748", lw=1.2)
        for c in kids[i]: ax.plot([o["t_inf"][c]] * 2, [ypos[i], ypos[c]], color="#cbd5e0" if (sampled and (i not in letters and c not in letters)) else "#2d3748", lw=1.0)
    for i in range(n):
        s = sampled and i in letters; ax.plot(o["t_inf"][i], ypos[i], "o", ms=9 if s else 7, color="#2b6cb0" if s else ("#fc8181" if not sampled else "#fed7d7"), zorder=3)
        if s: ax.text(max([o["t_inf"][c] for c in kids[i]] + [o["t_inf"][i] + 7]) + 3.2, ypos[i], letters[i], fontsize=11, fontweight="bold", va="center", color="#2b6cb0")
    if stars:
        allm = [m for g in G for m in g]; cnt = {m: allm.count(m) for m in set(allm)}; info = [m for m, c in cnt.items() if 1 <= c < 6]; cols = ["#d53f8c", "#ed8936", "#38b2ac", "#805ad5", "#e53e3e", "#3182ce", "#d69e2e", "#2f855a"]
        for k, m in enumerate(sorted(info)):
            for i in range(n):
                tt, ids = o["events"][i]
                if m in ids: ax.plot(tt[ids.index(m)], ypos[i], "*", ms=11, color=cols[k % len(cols)], mec="black", mew=0.4, zorder=4)
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlabel("Time (days)")
fig = plt.figure(figsize=(14, 7.6)); gs = fig.add_gridspec(2, 6, hspace=0.42, wspace=0.7)
a1 = fig.add_subplot(gs[0, 0:2]); draw_network(a1); a1.set_title("1  The epidemic: every case in the transmission network", fontsize=10, loc="left")
a2 = fig.add_subplot(gs[0, 2:4]); draw_network(a2, sampled=True); a2.set_title("2  Only some cases are sequenced (A to F)", fontsize=10, loc="left")
a3 = fig.add_subplot(gs[0, 4:6]); draw_network(a3, sampled=True, stars=True); a3.set_title("3  Mutations (stars) mark the lineages", fontsize=10, loc="left")
Z = linkage(squareform(D), "average"); a4 = fig.add_subplot(gs[1, 0:3]); dn = dendrogram(Z, labels=[letters[c] for c in samp], orientation="right", ax=a4, color_threshold=0, above_threshold_color="#2d3748", leaf_font_size=12)
a4.set_xlabel("SNP distance between sequenced cases"); a4.set_title("4  The phylogeny reconstructed from the genomes", fontsize=10, loc="left"); a4.spines["left"].set_visible(False)
a5 = fig.add_subplot(gs[1, 3:6]); a5.axis("off"); leaves = dn["leaves"]; pos = {leaf: k for k, leaf in enumerate(leaves)}; root, nodelist = to_tree(Z, rd=True); nodes = {nd.id: nd for nd in nodelist}
def depth_of(nd, d=0, out=None):
    out = {} if out is None else out; out[nd.id] = d
    if not nd.is_leaf(): depth_of(nd.left, d + 1, out); depth_of(nd.right, d + 1, out)
    return out
dd = depth_of(root); maxd = max(dd.values()); lab = {c: letters[c] for c in samp}
def leaves_of(nd): return nd.pre_order()
for nid, nd in sorted(nodes.items(), key=lambda kv: dd[kv[0]]):
    if nd.is_leaf(): continue
    ls = [pos[x] for x in leaves_of(nd)]; d0 = dd[nid]; a5.add_patch(FancyBboxPatch((min(ls) - 0.42 + 0.05 * d0, -1.15 + 0.2 * d0), max(ls) - min(ls) + 0.84 - 0.1 * d0, 2.3 - 0.4 * d0, boxstyle="round,pad=0.02,rounding_size=0.25", fc="none", ec="#4a5568", lw=1.3))
for leaf, k in pos.items(): a5.text(k, 0, "ABCDEF"[leaf], ha="center", va="center", fontsize=12, fontweight="bold", color="#2b6cb0")
a5.set_xlim(-0.8, 5.8); a5.set_ylim(-1.4, 1.5); a5.set_title("5  The relatedness it implies: nested groups", fontsize=10, loc="left")
fig.savefig("figs/fig_epidemic_to_tree.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# ---------------------------------------------------------------- figure B: investigation uses
def building(ax, x, w, h, roof=True, fc="#edf2f7"):
    ax.add_patch(Rectangle((x, 0), w, h, fc=fc, ec="#a0aec0", lw=1.2))
    if roof: ax.add_patch(Polygon([[x - 0.1, h], [x + w / 2, h + 0.45], [x + w + 0.1, h]], fc="#e2e8f0", ec="#a0aec0", lw=1.2))
    for k in range(int(w / 0.5)): ax.add_patch(Rectangle((x + 0.18 + 0.5 * k, h * 0.55), 0.22, 0.28, fc="white", ec="#cbd5e0", lw=0.8))
COL = {"A": "#c53030", "B": "#d69e2e", "C": "#2b6cb0"}
def pins(ax, x, y, labels, dx=0.3):
    for k, l in enumerate(labels): ax.add_patch(Circle((x + dx * k, y), 0.2, fc="white", ec=COL[l], lw=2, zorder=5)); ax.text(x + dx * k, y, l, ha="center", va="center", fontsize=9, fontweight="bold", color=COL[l], zorder=6)
fig, ax = plt.subplots(1, 3, figsize=(14, 3.2))
for a, title in zip(ax, ["Hotspot of transmission", "Multiple introductions: transmission dispelled", "A supplement to contact tracing"]):
    a.set_xlim(0, 10); a.set_ylim(-0.8, 4.0); a.set_aspect("equal"); a.axis("off"); a.set_title(title, fontsize=11, fontweight="bold", color=INK, loc="left")
    building(a, 0.4, 2.6, 1.6); building(a, 3.6, 2.4, 1.5); building(a, 6.6, 3.0, 1.1); a.plot([0, 10], [-0.15, -0.15], color="#cbd5e0", lw=3)
    for x, lab in ((1.7, "Home"), (4.8, "Workplace"), (8.1, "Care home")): a.text(x, -0.5, lab, ha="center", fontsize=8.5, color=QUIET)
a = ax[0]; pins(a, 3.9, 0.9, list("BBBBB"), 0.36); pins(a, 0.9, 1.0, ["A"]); pins(a, 7.2, 0.7, ["A"]); a.add_patch(Circle((4.6, 1.0), 1.35, fc="none", ec="#2c7a7b", lw=2.5)); a.text(4.8, 3.3, "Five identical genomes in one setting:\nlikely shared transmission", ha="center", fontsize=9, color="#2c7a7b")
a = ax[1]; pins(a, 3.9, 0.9, list("ACB"), 0.45); pins(a, 0.9, 1.0, ["A"]); pins(a, 7.2, 0.7, ["A"]); a.add_patch(Circle((4.55, 1.0), 1.25, fc="none", ec="#2c7a7b", lw=2.5)); a.text(4.8, 3.3, "Three unrelated genomes in one setting:\nseparate introductions, no local chain", ha="center", fontsize=9, color="#2c7a7b")
a = ax[2]; pins(a, 0.9, 1.0, ["A"]); pins(a, 4.3, 0.9, ["A"]); pins(a, 7.5, 0.7, ["A"]); a.annotate("", xy=(4.15, 1.3), xytext=(1.1, 1.3), arrowprops=dict(arrowstyle="-", ls="--", color="#2c7a7b", lw=2)); a.annotate("", xy=(7.3, 1.1), xytext=(4.5, 1.2), arrowprops=dict(arrowstyle="-", ls="--", color="#2c7a7b", lw=2)); a.text(4.8, 3.3, "The same genome in settings with no known contact:\na link that interviews alone had missed", ha="center", fontsize=9, color="#2c7a7b")
fig.savefig("figs/fig_investigation_uses.png", dpi=200, bbox_inches="tight"); plt.close(fig)
# ---------------------------------------------------------------- figure C: workflow and platforms
fig, ax = plt.subplots(figsize=(11.2, 6.2)); ax.set_xlim(0, 100); ax.set_ylim(0, 62); ax.axis("off")
ax.text(1, 59, "From swab to shared sequence: RNA whole-genome sequencing", fontsize=12.5, fontweight="bold", color=INK, va="center")
steps = [("1 Sample", ["Swab or", "clinical specimen"]), ("2 Extract", ["Purify the", "viral RNA"]), ("3 Library", ["Make sequencing-", "ready DNA"]), ("4 Sequence", ["Read the", "library"]), ("5 Share", ["Upload, link to", "case data"])]
w, g = 17.6, 2.7; xs = [1 + i * (w + g) for i in range(5)]
for i, (t, l) in enumerate(steps):
    box(ax, xs[i], 43, w, 12.5, t, l, accent=(i == 2))
    if i < 4: arrow(ax, xs[i] + w, 49.2, xs[i] + w + g, 49.2)
box(ax, 1, 13, 55, 24, "Library preparation in detail", ["1. Reverse transcription: viral RNA to cDNA", "2. Multiplex PCR: copy the genome in overlapping pieces", "3. End repair: make the fragment ends uniform", "4. Indexing: add a barcode so samples can be pooled", "5. Adapter ligation (a motor protein for nanopores)", "Each step can fail: a dropped PCR piece leaves a gap (Chapter 3)"], fill="#ebf8ff", dashed=True)
ax.plot([xs[2] + w / 2, xs[2] + w / 2], [43, 37], color=EDGE, lw=1, ls=":")
box(ax, 60, 21, 39, 16, "Short-read platforms", ["Accurate reads of 100 to 300 bases", "High throughput; higher equipment cost", "Most widely used in public health"], fill="#ebf8ff"); box(ax, 60, 2, 39, 16, "Long-read platforms", ["Reads of thousands of bases", "Fast and portable; low capital cost", "Accuracy has risen (Chapter 3)"], fill="#ebf8ff")
ax.plot([xs[3] + w / 2, xs[3] + w / 2], [43, 37], color=EDGE, lw=1, ls=":"); ax.text(xs[3] + w / 2 + 1, 39.6, "step 4 uses one of two families", fontsize=8.2, color=QUIET)
fig.savefig("figs/fig_sequencing_workflow.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("figures written")
