"""
Chapter 8, Exercise 1 - From public Zika genomes to a time-scaled tree.

Steps: (1) select human Zika genomes with a country and a collection date from NCBI GenBank and sample 600 at random;
(2) align with MAFFT; (3) build a maximum-likelihood tree with FastTree; (4) date the tree with TreeTime and report
credible intervals for the nodes that matter. Outputs are cached, so reruns skip finished steps.
Needs: mafft, FastTree on the PATH; pip install biopython treetime pandas numpy
"""
import os, re, json, time, random, subprocess, urllib.parse, urllib.request
import numpy as np, pandas as pd
from Bio import SeqIO, Phylo
from treetime import TreeTime, utils

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
N_SAMPLE, SEED = 600, 8
EXCLUDE = {"Senegal", "Uganda", "Central African Republic", "Nigeria", "Gabon"}     # African lineage: a separate clock
MON = {m: i + 1 for i, m in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split())}

def get(url):
    for k in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r: return r.read().decode()
        except Exception: time.sleep(3 * (k + 1))
    raise RuntimeError(url)

# ---------------------------------------------------------------- 1. select and download
def parse_date(s):                                   # -> (year, month, day, precision)
    if not isinstance(s, str): return None
    for pat, f in [(r"(\d{4})-(\d{2})-(\d{2})", lambda m: (int(m[1]), int(m[2]), int(m[3]), "day")),
                   (r"(\d{2})-([A-Za-z]{3})-(\d{4})", lambda m: (int(m[3]), MON[m[2]], int(m[1]), "day")),
                   (r"(\d{4})-(\d{2})", lambda m: (int(m[1]), int(m[2]), 15, "month")),
                   (r"([A-Za-z]{3})-(\d{4})", lambda m: (int(m[2]), MON[m[1]], 15, "month"))]:
        m = re.fullmatch(pat, s.strip())
        if m: return f(m)
    return None

if not os.path.exists("zika_selected.fasta"):
    ids = json.loads(get(E + "esearch.fcgi?db=nuccore&retmode=json&retmax=5000&term="
                         + urllib.parse.quote('"Zika virus"[Organism] AND 10000:11000[SLEN]')))["esearchresult"]["idlist"]
    rows = []
    for i in range(0, len(ids), 200):
        r = json.loads(get(E + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ids[i:i + 200])))["result"]
        for u in r["uids"]:
            f = dict(zip(r[u].get("subtype", "").split("|"), r[u].get("subname", "").split("|")))
            rows.append(dict(uid=u, acc=r[u]["accessionversion"], country=f.get("country"), date=f.get("collection_date"), host=f.get("host")))
        time.sleep(0.4)
    df = pd.DataFrame(rows); df.to_csv("zika_all_meta.csv", index=False)
    pdt = df.date.map(parse_date)
    df["dprec"] = pdt.map(lambda t: t[3] if t else None)
    for k, name in enumerate(["year", "month", "day"]): df[name] = pdt.map(lambda t: t[k] if t else np.nan)
    df["cty"] = df.country.str.split(":").str[0].str.strip()
    df["human"] = df.host.fillna("").str.startswith("Homo sapiens")
    sel = df[df.cty.notna() & df.dprec.notna() & df.human & df.year.between(2006, 2018) & ~df.cty.isin(EXCLUDE)]
    print(f"{len(df)} Zika records; {len(sel)} human genomes with a country and a collection date (2006-2018)")
    pick = sel.sample(n=min(N_SAMPLE, len(sel)), random_state=SEED).sort_values("acc")
    pick.to_csv("zika_selected_meta.csv", index=False)
    fasta = ""
    for i in range(0, len(pick), 50):
        fasta += get(E + "efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(pick.uid.iloc[i:i + 50])); time.sleep(0.5)
    open("zika_selected.fasta", "w").write(fasta)
meta = pd.read_csv("zika_selected_meta.csv", dtype=str)
print(f"selected {len(meta)} genomes from {meta.cty.nunique()} countries")

# ---------------------------------------------------------------- 2-3. align and build the tree
if not os.path.exists("zika_aln.fasta"):
    with open("zika_aln.fasta", "w") as out:
        subprocess.run(["mafft", "--retree", "2", "--maxiterate", "0", "--thread", "-1", "--quiet", "zika_selected.fasta"], stdout=out, check=True)
if not os.path.exists("zika_ml.nwk"):
    with open("zika_ml.nwk", "w") as out:
        subprocess.run(["FastTree", "-nt", "-gtr", "-quiet", "zika_aln.fasta"], stdout=out, check=True)

# ---------------------------------------------------------------- 4. date the tree
meta["tt_date"] = [f"{int(float(r.year))}-{int(float(r.month)):02d}-{int(float(r.day)):02d}" if r.dprec == "day"
                   else f"{int(float(r.year))}-{int(float(r.month)):02d}-XX" for r in meta.itertuples()]    # XX = day unknown
meta[["acc", "tt_date"]].rename(columns={"acc": "name", "tt_date": "date"}).to_csv("tt_dates.csv", index=False)
country = meta.set_index("acc").cty

def date_tree(fixed_rate=None, marginal=False):
    dates = utils.parse_dates("tt_dates.csv", name_col="name", date_col="date")
    tt = TreeTime(dates=dates, tree="zika_ml.nwk", aln="zika_aln.fasta", gtr="JC69", verbose=0)
    tt.clock_filter(reroot="least-squares", n_iqd=3, plot=False)             # flag tips whose date and branch length disagree
    bad = [n.name for n in tt.tree.get_terminals() if getattr(n, "bad_branch", False)]
    tt.run(root="least-squares", infer_gtr=True, max_iter=2, resolve_polytomies=True, branch_length_mode="joint",
           fixed_clock_rate=fixed_rate, time_marginal=False, verbose=0)
    if marginal:                                                              # posterior node dates, needed for intervals
        tt.run(infer_gtr=False, max_iter=1, time_marginal=True, resolve_polytomies=False, fixed_clock_rate=fixed_rate, verbose=0)
    return tt, bad

def brazil_nodes(tt):
    bz = [t for t in tt.tree.get_terminals() if country.get(t.name) == "Brazil" and not getattr(t, "bad_branch", False)]
    anc = tt.tree.common_ancestor(bz)
    return anc, tt.tree.get_path(anc)[-2]

if not os.path.exists("timetree_summary.json"):
    tt, bad = date_tree(marginal=True)
    anc, parent = brazil_nodes(tt)
    summary = dict(clock_rate=float(tt.date2dist.clock_rate), root_date=float(tt.tree.root.numdate), n_tips=tt.tree.count_terminals(), n_outliers=len(bad),
                   root_ci90=[float(x) for x in tt.get_max_posterior_region(tt.tree.root, 0.9)],
                   brazil_mrca=float(anc.numdate), brazil_mrca_ci90=[float(x) for x in tt.get_max_posterior_region(anc, 0.9)],
                   brazil_parent=float(parent.numdate), brazil_parent_ci90=[float(x) for x in tt.get_max_posterior_region(parent, 0.9)])
    for c in tt.tree.find_clades():                                           # branch lengths in years, so the tree reads as calendar time
        path = tt.tree.get_path(c)
        c.branch_length = 0.0 if c is tt.tree.root else max(c.numdate - (path[-2].numdate if len(path) > 1 else tt.tree.root.numdate), 0.0)
    Phylo.write(tt.tree, "zika_timetree.nwk", "newick")
    pd.DataFrame([(c.name, c.numdate, c.is_terminal()) for c in tt.tree.find_clades()], columns=["node", "numdate", "terminal"]).to_csv("zika_node_dates.csv", index=False)
    pd.DataFrame({"acc": bad}).to_csv("zika_clock_outliers.csv", index=False)
    json.dump(summary, open("timetree_summary.json", "w"), indent=1)
s = json.load(open("timetree_summary.json"))
print(f"clock rate {s['clock_rate']:.2e} substitutions/site/year; {s['n_outliers']} tips flagged as clock outliers")
print(f"root {s['root_date']:.1f} (90% interval {s['root_ci90'][0]:.1f}-{s['root_ci90'][1]:.1f})")
print(f"ancestor of the Brazilian genomes {s['brazil_mrca']:.2f} (90% interval {s['brazil_mrca_ci90'][0]:.2f}-{s['brazil_mrca_ci90'][1]:.2f}); "
      f"its parent {s['brazil_parent']:.2f} (90% interval {s['brazil_parent_ci90'][0]:.2f}-{s['brazil_parent_ci90'][1]:.2f})")

# ---------------------------------------------------------------- 5. how much do the dates depend on the clock rate?
if not os.path.exists("rate_sensitivity.csv"):
    rows = [dict(clock_rate=s["clock_rate"], setting="estimated", root=s["root_date"], brazil_mrca=s["brazil_mrca"], brazil_parent=s["brazil_parent"])]
    for r in (8e-4, 1.0e-3, 1.2e-3):
        tt, _ = date_tree(fixed_rate=r); anc, parent = brazil_nodes(tt)
        rows.append(dict(clock_rate=r, setting="fixed", root=float(tt.tree.root.numdate), brazil_mrca=float(anc.numdate), brazil_parent=float(parent.numdate)))
    pd.DataFrame(rows).to_csv("rate_sensitivity.csv", index=False)
print(pd.read_csv("rate_sensitivity.csv").round(3).to_string(index=False))
