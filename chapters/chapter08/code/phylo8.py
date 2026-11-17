import pandas as pd, numpy as np, time, warnings
warnings.filterwarnings("ignore")
from Bio import Phylo
from treetime.wrappers import reconstruct_discrete_traits

AMERICAS = {"Brazil","Nicaragua","Puerto Rico","Colombia","Mexico","USA","Cuba","Honduras","Dominican Republic","Haiti","Panama","Virgin Islands","Suriname","Guatemala","Jamaica","Venezuela","Ecuador","Peru","Trinidad and Tobago","Argentina","El Salvador","French Guiana","Canada","Martinique","Guadeloupe"}
ASIA_PAC = {"Singapore","Thailand","China","French Polynesia","Viet Nam","India","South Korea","Japan","Taiwan","Indonesia","Australia","Philippines","Micronesia","Cambodia","Bangladesh"}
def region_of(c): return "Americas" if c in AMERICAS else ("Asia-Pacific" if c in ASIA_PAC else "Europe/Africa")

def load_states(meta="zika_selected_meta.csv", min_tips=8):
    """One state per tip: the country if it has at least min_tips genomes, otherwise a regional bucket."""
    m = pd.read_csv(meta, dtype=str); m["cty"] = m.cty.replace({"Viet Nam": "Viet Nam"})
    counts = m.cty.value_counts()
    def state(c):
        if counts[c] >= min_tips: return c
        return {"Americas": "Other Americas", "Asia-Pacific": "Other Asia-Pacific", "Europe/Africa": "Other Europe/Africa"}[region_of(c)]
    m["state"] = m.cty.map(state)
    return m.set_index("acc")

def load_tree(drop_outliers=True):
    """Time-scaled tree (branch lengths in years) and the date of the root."""
    t = Phylo.read("zika_timetree.nwk", "newick")
    root_date = float(pd.read_csv("zika_node_dates.csv").numdate.iloc[0])      # first row = root
    if drop_outliers:
        for acc in pd.read_csv("zika_clock_outliers.csv").acc: t.prune(acc)
    return t, root_date

def clade_dates(tree, root_date):
    """Calendar date of every node = root date + distance from the root (branch lengths are years)."""
    return {id(c): root_date + d for c, d in tree.depths().items()}

def run_mugration(tree, traits, sbc=None, iterations=5, seed=1):
    import copy
    tr = copy.deepcopy(tree)
    mug, l2s, rev = reconstruct_discrete_traits(tr, traits, missing_data="?", pc=1.0, sampling_bias_correction=sbc, iterations=iterations, rng_seed=seed, verbose=0)
    states = {}; post = {}
    for n in mug.tree.find_clades():
        prof = n.marginal_profile[0]; k = int(np.argmax(prof))
        states[id(n)] = l2s[mug.gtr.alphabet[k]]; post[id(n)] = float(prof[k])
    return mug, states, post

def events(mug, states, dates_by_clade):
    """A migration event is a branch whose two ends have different inferred states; its date lies between the parent and child dates."""
    ev = []
    parent = {}
    for c in mug.tree.find_clades(order="level"):
        for ch in c.clades: parent[id(ch)] = c
    for c in mug.tree.find_clades():
        p = parent.get(id(c))
        if p is None: continue
        if states[id(p)] != states[id(c)]:
            dp, dc = dates_by_clade[id(p)], dates_by_clade[id(c)]
            ev.append(dict(src=states[id(p)], dst=states[id(c)], date_lo=dp, date_hi=dc, date=(dp + dc) / 2, child_tips=c.count_terminals()))
    return pd.DataFrame(ev)

def is_americas(state): return state in AMERICAS or state == "Other Americas"
