"""Chapter 5 diagram: reading a tree. One toy tree drawn as a cladogram, as a phylogram and as the same phylogram with two nodes rotated."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
INK, MUTED, BLUE, RED = "#1a202c", "#4a5568", "#2b6cb0", "#c53030"
# toy tree: (name or children, branch length in substitutions per site, tip label with sampling date)
T = ("root", 0, [("n1", 0.0030, [("A", 0.0020, "A  May"), ("B", 0.0010, "B  Jun")]), ("n2", 0.0010, [("C", 0.0030, "C  Jun"), ("n3", 0.0015, [("D", 0.0005, "D  Jul"), ("E", 0.0015, "E  Jul")])])])
def rotate(node, names):
    nm, bl, kids = node
    if isinstance(kids, str): return node
    k = [rotate(c, names) for c in kids]
    if nm in names: k = k[::-1]
    return (nm, bl, k)
def layout(node, x, ys, pos, clad):
    nm, bl, kids = node; xx = x + (1.0 if clad else bl)
    if isinstance(kids, str):
        y = len(ys); ys.append(y); pos[nm] = (xx, y, x); return y
    cy = [layout(c, xx, ys, pos, clad) for c in kids]; y = sum(cy) / len(cy); pos[nm] = (xx, y, x); return y
def depth(node, clad, x=0):
    nm, bl, kids = node; xx = x + (1.0 if clad else bl); return xx if isinstance(kids, str) else max(depth(c, clad, xx) for c in kids)
def draw(ax, node, clad, title, sub, note=None):
    pos = {}; layout((node[0], 0, node[2]), 0, [], pos, clad); sc = depth(node, clad); X = lambda v: v / sc * 0.78
    def walk(n, parent=None):
        nm, bl, kids = n; xx, y, x0 = pos[nm]
        if parent is not None: ax.plot([X(x0), X(xx)], [y, y], color=INK, lw=2.2, solid_capstyle="round")
        if not isinstance(kids, str):
            ys = [pos[c[0]][1] for c in kids]; ax.plot([X(xx), X(xx)], [min(ys), max(ys)], color=INK, lw=2.2)
            for c in kids: walk(c, n)
            ax.plot(X(xx), y, "o", color="white", mec=INK, mew=1.6, ms=7, zorder=3)
        else: ax.text(X(xx) + 0.02, y, kids, va="center", ha="left", fontsize=11.5, color=INK, fontweight="bold")
    walk((node[0], 0, node[2])); ax.set_xlim(-0.04, 1.1); ax.set_ylim(-0.7, 5.2); ax.axis("off"); ax.text(0, 5.45, title, fontsize=13.5, fontweight="bold", color=INK, va="bottom"); ax.text(0, 5.05, sub, fontsize=10.2, color=MUTED, va="bottom")
    if note: ax.text(0, -0.55, note, fontsize=10, color=RED, va="top")
    return pos, X
fig = plt.figure(figsize=(15.28, 7.2), dpi=100); fig.suptitle("Reading a tree: what is on the page and what is not", x=0.015, y=0.975, ha="left", fontsize=15.5, fontweight="bold", color=INK)
axs = [fig.add_axes([0.01 + i * 0.335, 0.2, 0.32, 0.66]) for i in range(3)]
draw(axs[0], T, True, "Cladogram", "Topology only: who joins whom")
pos, X = draw(axs[1], T, False, "Phylogram", "Branch length: substitutions per site")
axs[1].plot([X(0), X(0.002)], [-0.35, -0.35], color=BLUE, lw=3); axs[1].text(X(0.001), -0.55, "0.002", ha="center", fontsize=10, color=BLUE)
R = rotate(T, {"n1", "n2", "n3"}); draw(axs[2], R, False, "The same phylogram, rotated", "Nodes n1, n2 and n3 swapped: same tree", None)
ban = fig.add_axes([0.01, 0.02, 0.98, 0.15]); ban.axis("off"); ban.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.0,rounding_size=0.04", fc="#edf2f7", ec="#cbd5e0", transform=ban.transAxes))
ban.text(0.015, 0.74, "Tips are sampled sequences; the circles are inferred ancestors, none of them sampled. The distance between two tips is the sum of the branches on the path between them.", transform=ban.transAxes, fontsize=11.2, color=INK, va="center")
ban.text(0.015, 0.30, "Tip order on the page is arbitrary: neighbouring tips are not necessarily linked cases, and a sampled sequence is never the ancestor of another sampled sequence.", transform=ban.transAxes, fontsize=11.2, color=INK, va="center", fontweight="bold")
fig.savefig("diagram_reading_tree.png", dpi=100, facecolor="white"); print("ok")
