"""
Chapter 8, Hands-on 5 - does the sampling strategy recover the true source? A simulated epidemic with a known truth.

An epidemic is simulated forward in time in four places (A to D) with migration between them: in the 'hub' scenario place A has the highest transmission rate and the highest
export rate and B, C and D are recipients; in the 'no_hub' scenario no place is a hub. Infected hosts are sequenced when they are removed, with a probability that is four times higher in place D than elsewhere
(a convenience sample). The genealogy of the sequenced hosts is recorded with the true place of every node, so the true migration events are known.
Four strategies are applied to the same sequenced pool and the places are reconstructed on the time tree with TreeTime's discrete-trait model (the function of the chapter's other exercises):
  as_sequenced   all sequenced hosts (the natural imbalance)
  cap_per_place  at most CAP hosts per place (balanced counts)
  proportional   every place thinned to the lowest sequencing probability, so the sample is proportional to the number of infections
  as_sequenced_sbc   as_sequenced with TreeTime's sampling_bias_correction = 2
Usage: python ex8_6_sampling_truth.py [reps]      (default 60 epidemics in each of the two scenarios, about 7 minutes on one core; finished epidemics are kept in sampling_truth_runs.csv and skipped on a restart)
Writes sampling_truth_runs.csv: one row per epidemic and strategy with the tips, the true and the inferred outgoing events of each place and the root.
"""
import os, sys, copy, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from Bio.Phylo.BaseTree import Clade, Tree
from phylo8 import run_mugration
sys.setrecursionlimit(20000)
NAMES = ["A", "B", "C", "D"]; K = 4
GAMMA = 1.0; N_CUM = 4000
SCEN = {"hub": dict(beta=[1.8, 1.5, 1.5, 1.5], mig=np.array([[0.92, 0.0267, 0.0267, 0.0266], [0.02, 0.94, 0.02, 0.02], [0.02, 0.02, 0.94, 0.02], [0.02, 0.02, 0.02, 0.94]])),       # row = place of the infector, column = place of the infectee; A exports 8% of its transmissions, B to D 6%
        "no_hub": dict(beta=[1.5, 1.5, 1.5, 1.5], mig=np.array([[0.94, 0.02, 0.02, 0.02], [0.02, 0.94, 0.02, 0.02], [0.02, 0.02, 0.94, 0.02], [0.02, 0.02, 0.02, 0.94]]))}          # no hub: every place exports 6%
P_SEQ = np.array([0.04, 0.04, 0.04, 0.16]); CAP = 15; MIN_POOL = 120
def simulate(seed, scen):
    """Forward simulation; returns the hosts (parent, birth, place, removal, u) after N_CUM infections. A sequenced host is one removed before the end with u < P_SEQ[place]."""
    rng = np.random.default_rng(seed); BETA = SCEN[scen]["beta"]; MIG = SCEN[scen]["mig"]
    parent = [-1]; birth = [0.0]; place = [0]; removal = [np.nan]; u = [rng.random()]
    inf = [[0], [], [], []]; t = 0.0
    while len(parent) < N_CUM:
        w = np.array([len(inf[d]) * (BETA[d] + GAMMA) for d in range(K)])
        if w.sum() == 0: return None
        t += rng.exponential(1 / w.sum()); d = rng.choice(K, p=w / w.sum()); j = rng.integers(len(inf[d])); i = inf[d][j]
        if rng.random() < BETA[d] / (BETA[d] + GAMMA):
            d2 = rng.choice(K, p=MIG[d]); h = len(parent); parent.append(i); birth.append(t); place.append(int(d2)); removal.append(np.nan); u.append(rng.random()); inf[d2].append(h)
        else:
            removal[i] = t; inf[d][j] = inf[d][-1]; inf[d].pop()
    return dict(parent=parent, birth=birth, place=place, removal=removal, u=u)
def genealogy(sim, keep):
    """Tree of the hosts in keep (sampled at their removal), with the true place of every node and branch lengths in time."""
    n = len(sim["parent"]); kids = [[] for _ in range(n)]
    for h in range(1, n): kids[sim["parent"][h]].append((sim["birth"][h], h))
    def build(i):
        trunk = None
        if i in keep: trunk = Clade(name=f"h{i}"); trunk.state = sim["place"][i]; trunk.time = sim["removal"][i]
        for tt, c in reversed(kids[i]):
            sub = build(c)
            if sub is None: continue
            if trunk is None: trunk = sub
            else:
                node = Clade(); node.state = sim["place"][i]; node.time = tt; node.clades = [trunk, sub]; trunk = node
        return trunk
    root = build(0)
    def setlen(c, par):
        c.branch_length = 0.0 if par is None else max(c.time - par.time, 1e-9)
        for ch in c.clades: setlen(ch, c)
    setlen(root, None); return Tree(root=root, rooted=True)
def events(tree, states=None):
    """Outgoing events: matrix[src, dst] over branches whose parent and child nodes are in different places (true places, or the inferred ones in states)."""
    M = np.zeros((K, K), int); par = {}
    for c in tree.find_clades():
        for ch in c.clades: par[id(ch)] = c
    for c in tree.find_clades():
        p = par.get(id(c))
        if p is None: continue
        lab = lambda x: NAMES.index(x) if isinstance(x, str) else int(x)
        a, b = (p.state, c.state) if states is None else (lab(states[id(p)]), lab(states[id(c)]))
        if a != b: M[a, b] += 1
    return M
def prune_to(tree, keep_names):
    tr = copy.deepcopy(tree)
    for t in [t.name for t in tr.get_terminals() if t.name not in keep_names]: tr.prune(t)
    return tr
def strategy_sets(sim, pool, rng):
    pl = np.array([sim["place"][i] for i in pool]); sets = {"as_sequenced": list(pool), "as_sequenced_sbc": list(pool)}
    cap = []; prop = []
    for d in range(K):
        mem = [i for i, p in zip(pool, pl) if p == d]; cap += list(rng.choice(mem, size=min(CAP, len(mem)), replace=False)) if mem else []
        prop += [i for i in mem if sim["u"][i] < P_SEQ.min()]          # the same uniform draw that decided sequencing: u < lowest probability keeps every place at the same sampling fraction
    sets["cap_per_place"] = cap; sets["proportional"] = prop; return sets
def one_epidemic(rep, scen):
    for attempt in range(200):
        seed = rep * 1000 + attempt + (0 if scen == "hub" else 500); sim = simulate(seed, scen)
        if sim is None: continue
        pool = [i for i in range(len(sim["parent"])) if not np.isnan(sim["removal"][i]) and sim["u"][i] < P_SEQ[sim["place"][i]]]
        if len(pool) >= MIN_POOL and all(sum(1 for i in pool if sim["place"][i] == d) >= 12 for d in range(K)): break
    else: raise RuntimeError("no usable epidemic")
    full = genealogy(sim, set(pool)); rng = np.random.default_rng(seed + 7); rows = []
    for name, ids in strategy_sets(sim, pool, rng).items():
        keep = {f"h{i}" for i in ids}; tr = prune_to(full, keep); traits = {f"h{i}": NAMES[sim["place"][i]] for i in ids}
        mug, st, post = run_mugration(tr, traits, sbc=2 if name.endswith("sbc") else None, iterations=5, seed=1)
        # run_mugration works on a copy of the tree, so the true places are read from that copy
        inf_M = events(mug.tree, st); true_M = events(mug.tree)
        root = mug.tree.root; r = dict(scenario=scen, rep=rep, seed=seed, strategy=name, n_tips=len(ids), true_root=int(root.state), inf_root=int(NAMES.index(st[id(root)])) if isinstance(st[id(root)], str) else int(st[id(root)]))
        for d in range(K): r[f"tips_{NAMES[d]}"] = sum(1 for i in ids if sim["place"][i] == d)
        for d in range(K): r[f"true_out_{NAMES[d]}"] = int(true_M[d].sum()); r[f"inf_out_{NAMES[d]}"] = int(inf_M[d].sum())
        r["true_events"] = int(true_M.sum()); r["inf_events"] = int(inf_M.sum()); rows.append(r)
    return rows
if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 60; FILE = "sampling_truth_runs.csv"; done = pd.read_csv(FILE) if os.path.exists(FILE) else pd.DataFrame()
    for scen in SCEN:
        for rep in range(1, reps + 1):
            if len(done) and ((done.scenario == scen) & (done.rep == rep)).any(): continue
            t0 = time.time(); rows = one_epidemic(rep, scen); done = pd.concat([done, pd.DataFrame(rows)], ignore_index=True); done.to_csv(FILE, index=False)
            print(f"{scen} epidemic {rep}: {rows[0]['n_tips']} sequenced hosts, {time.time() - t0:.0f} s", flush=True)
