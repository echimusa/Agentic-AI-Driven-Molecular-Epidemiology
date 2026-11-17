"""
Chapter 8, Hands-on 2 (extension) - is the dated tree itself reproducible?

The dating step of ex8_1_dated_tree.py (TreeTime on the cached alignment and maximum-likelihood tree) is repeated N times with nothing changed. For each rerun the script records
the clock rate, the date of the root, the date of the common ancestor of the Brazilian genomes and of its parent, then reruns the place reconstruction of ex8_2_migration.py on that rerun's
tree and records the number of inferred migration events, the number of introductions into the Americas and the interval and size of the main introduction into Brazil.
Usage: python ex8_5_dating_reruns.py [N]      (default 10 reruns, about 100 s each; finished reruns are kept in dating_reruns.csv and skipped on a restart)
Needs the cached files of the chapter's data folder (zika_selected_meta.csv, tt_dates.csv, zika_aln.fasta, zika_ml.nwk) and phylo8.py.
"""
import os, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from Bio import Phylo
from treetime import TreeTime, utils
import phylo8
FILE = "dating_reruns.csv"
meta = pd.read_csv("zika_selected_meta.csv", dtype=str); country = meta.set_index("acc").cty
def date_tree():
    dates = utils.parse_dates("tt_dates.csv", name_col="name", date_col="date")
    tt = TreeTime(dates=dates, tree="zika_ml.nwk", aln="zika_aln.fasta", gtr="JC69", verbose=0)
    tt.clock_filter(reroot="least-squares", n_iqd=3, plot=False)
    bad = [n.name for n in tt.tree.get_terminals() if getattr(n, "bad_branch", False)]
    tt.run(root="least-squares", infer_gtr=True, max_iter=2, resolve_polytomies=True, branch_length_mode="joint", time_marginal=False, verbose=0)
    tt.run(infer_gtr=False, max_iter=1, time_marginal=True, resolve_polytomies=False, verbose=0)
    return tt, bad
def one_rerun(k):
    t0 = time.time(); tt, bad = date_tree()
    bz = [t for t in tt.tree.get_terminals() if country.get(t.name) == "Brazil" and not getattr(t, "bad_branch", False)]; anc = tt.tree.common_ancestor(bz); parent = tt.tree.get_path(anc)[-2]
    row = dict(rerun=k, clock_rate=float(tt.date2dist.clock_rate), root=float(tt.tree.root.numdate), brazil_mrca=float(anc.numdate), brazil_parent=float(parent.numdate), n_outliers=len(bad))
    for c in tt.tree.find_clades():                                    # branch lengths in years, as in ex8_1
        path = tt.tree.get_path(c); c.branch_length = 0.0 if c is tt.tree.root else max(c.numdate - (path[-2].numdate if len(path) > 1 else tt.tree.root.numdate), 0.0)
    tree = tt.tree
    for acc in bad:
        tree.prune(acc)
    root_date = float(tt.tree.root.numdate)
    st = phylo8.load_states(); traits = {t.name: st.state[t.name] for t in tree.get_terminals()}
    mug, states, post = phylo8.run_mugration(tree, traits)
    dates = phylo8.clade_dates(mug.tree, root_date); par = {}
    for c in mug.tree.find_clades():
        for ch in c.clades: par[id(ch)] = c
    americas_states = {s for s in set(traits.values()) if s in phylo8.AMERICAS or s == "Other Americas"}; ev = []
    for c in mug.tree.find_clades():
        p = par.get(id(c))
        if p is not None and states[id(p)] != states[id(c)]: ev.append((states[id(p)], states[id(c)], dates[id(p)], dates[id(c)], c.count_terminals()))
    intro = [e for e in ev if e[1] in americas_states and e[0] not in americas_states]; main = max(intro, key=lambda e: e[4]) if intro else None
    row.update(n_events=len(ev), n_intro=len(intro), main_src=main[0] if main else "", main_dst=main[1] if main else "", main_lo=main[2] if main else np.nan, main_hi=main[3] if main else np.nan, main_tips=main[4] if main else 0, seconds=round(time.time() - t0))
    return row
if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10; done = pd.read_csv(FILE) if os.path.exists(FILE) else pd.DataFrame()
    for k in range(1, N + 1):
        if len(done) and k in set(done.rerun): continue
        r = one_rerun(k); done = pd.concat([done, pd.DataFrame([r])], ignore_index=True); done.to_csv(FILE, index=False)
        print(f"rerun {k}: rate {r['clock_rate']:.3e}, Brazil ancestor {r['brazil_mrca']:.2f}, {r['n_events']} events, {r['n_intro']} introductions, main {r['main_lo']:.2f}-{r['main_hi']:.2f} ({r['main_tips']} tips), {r['seconds']} s", flush=True)
