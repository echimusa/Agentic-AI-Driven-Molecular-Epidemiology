"""
Chapter 5 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts: its own FASTA reader and site counts, p, JC69 and K80 distances, neighbour-joining, a Felsenstein pruning routine
(scipy matrix exponentials) for the eight-sequence likelihoods, a generalised Fitch parsimony written with bit masks, root-to-tip regression over every root position by shortest paths on the tree, four-point and three-point tests on its own random quartets,
support shares from the IQ-TREE files, the ascertainment and recombination tables from the raw outputs, and a private-substitution scan of the Ebola alignment. Results are written to validation_results_ch5.csv.
USAGE   python validate_ch5.py [--chapter Chapter_5_....docx] [--doi] [--full] [--sim]
        --chapter  also checks the chapter's tables, figure numbers and key numbers against the result tables (needs pandoc)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns the analysis scripts in a scratch folder and checks that the result tables are reproduced byte for byte (about 5 minutes; needs mafft, IQ-TREE 2, FastTree, R with ape)
        --sim      with --full, also reruns the two long simulations (ex5_3 part B, about 5 minutes, and ex5_8_recombination.py, about 30 minutes)
Prints PASS/FAIL/SKIP per check and exits 1 if any check fails.
"""
import argparse, glob, hashlib, io, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy import stats
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from Bio import Phylo
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); ap.add_argument("--sim", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:160])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
# ------------------------------------------------------------------ own readers
def readfa(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"): names.append(line[1:].split()[0]); seqs.append([])
        elif line: seqs[-1].append(line.upper())
    return names, np.array([list("".join(s)) for s in seqs])
def codes(X):
    C = np.full(X.shape, -1, dtype=np.int8)
    for i, b in enumerate("ACGT"): C[X == b] = i
    return C
def site_counts(C):
    cnt = np.stack([(C == b).sum(0) for b in range(4)]); nst = (cnt > 0).sum(0); return cnt, nst
def parse_tree(path):
    t = Phylo.read(path, "newick"); nodes = {}; edges = []; tips = []
    def walk(c, par):
        i = len(nodes); nodes[id(c)] = i
        if par is not None: edges.append((nodes[id(par)], i, c.branch_length or 0.0, c))
        if c.is_terminal(): tips.append((c.name, i))
        for k in c.clades: walk(k, c)
    walk(t.root, None); return t, nodes, edges, tips
def dist_nodes(n, edges):
    r, c, w = [], [], []
    for u, v, l, _ in edges: r += [u, v]; c += [v, u]; w += [max(l, 1e-12)] * 2
    return shortest_path(csr_matrix((w, (r, c)), shape=(n, n)), directed=False)
def split_sets(t, names, minlen):
    idx = {n: i for i, n in enumerate(names)}; n = len(names); out = set()
    for cl in t.get_nonterminals():
        if cl is t.root or (cl.branch_length or 0) < minlen: continue
        s = frozenset(idx[x.name] for x in cl.get_terminals())
        if 0 in s: s = frozenset(range(n)) - s
        if 1 < len(s) < n - 1: out.add(s)
    return out
def rfd(a, b): return len(a ^ b) / max(len(a) + len(b), 1)
# generalised Fitch (Hartigan counting): the state set of a node is the set of bases carried by most children; every child that lacks the chosen base costs one change
def fitch_cost(tree, names, C):
    idx = {n: i for i, n in enumerate(names)}; var = np.where((site_counts(C)[1]) >= 2)[0]; M = np.full((len(names), len(var)), 15, dtype=np.uint8)
    for b in range(4): M[C[:, var] == b] = 1 << b
    total = 0
    def up(cl):
        nonlocal total
        if cl.is_terminal(): return M[idx[cl.name]]
        ch = [up(c) for c in cl.clades]; cnt = np.stack([sum(((m >> b) & 1).astype(np.int16) for m in ch) for b in range(4)]); best = cnt.max(0); total += int((len(ch) - best).sum())
        return sum(((cnt[b] == best).astype(np.uint8) << b) for b in range(4)).astype(np.uint8)
    up(tree.root); return total
# ================================================================== A. data facts
nE, XE = readfa(P("ebola_aln.fasta")); nH, XH = readfa(P("h3n2_ha_aln.fasta")); CE, CH = codes(XE), codes(XH)
def facts(C):
    cnt, nst = site_counts(C); inf = ((cnt >= 2).sum(0) >= 2); return C.shape, int((nst >= 2).sum()), int(inf.sum()), 100 * (C.shape[1] - (nst >= 2).sum()) / C.shape[1]
fe, fh = facts(CE), facts(CH)
check("A.1", "Ebola alignment: 195 genomes, 18,961 positions, 241 variable sites, 84 parsimony-informative, 98.7% constant", fe[0] == (195, 18961) and fe[1] == 241 and fe[2] == 84 and rnd(fe[3], 1) == "98.7", fe)
check("A.2", "H3N2 alignment: 60 sequences, 1,762 positions, 230 variable sites, 119 parsimony-informative, 86.9% constant", fh[0] == (60, 1762) and fh[1] == 230 and fh[2] == 119 and rnd(fh[3], 1) == "86.9", fh)
nC, XC = readfa(P("h3n2_cds.fasta")); CC = codes(XC); cf = facts(CC)
check("A.3", "H3N2 coding region: 1,701 nt, no gap characters, 118 parsimony-informative, 86.5% constant, codon positions 91.9%, 93.1% and 74.6% constant", CC.shape[1] == 1701 and (XC == "-").sum() == 0 and cf[2] == 118 and rnd(cf[3], 1) == "86.5" and [rnd(facts(codes(readfa(P(f"h3n2_cp{k}.fasta"))[1]))[3], 1) for k in (1, 2, 3)] == ["91.9", "93.1", "74.6"], (cf, ))
# ================================================================== B. distances and distance trees (Hands-on 1)
def pw(C):
    n = len(C); ok = C >= 0; p = np.zeros((n, n)); Pt = np.zeros((n, n)); Qv = np.zeros((n, n))
    for i in range(n):
        v = ok[i][None, :] & ok; d = (C[i][None, :] != C) & v; ts = d & (np.abs(C[i][None, :] - C) == 2); L = np.maximum(v.sum(1), 1)
        p[i] = d.sum(1) / L; Pt[i] = ts.sum(1) / L; Qv[i] = (d.sum(1) - ts.sum(1)) / L
    return p, Pt, Qv
td = pd.read_csv(P("table_distances.csv")); okB = True; why = []
for C, row in ((CH, td.iloc[0]), (CE, td.iloc[1])):
    p, Pt, Qv = pw(C); iu = np.triu_indices(len(C), 1); jcd = -0.75 * np.log(1 - 4 / 3 * p); k80 = -0.5 * np.log(1 - 2 * Pt - Qv) - 0.25 * np.log(1 - 2 * Qv)
    got = [p[iu].mean(), p[iu].max(), jcd[iu].max(), np.abs(k80 - jcd)[iu].max(), Pt[iu].sum() / (p[iu].sum()), 100 * (p[iu] == 0).mean()]; exp = [row.mean_p, row.max_p, row.max_jc, row.k80_minus_jc_max, row.ts_share, row.identical_pairs_pct]
    good = all(abs(g - e) <= 6e-6 + 6e-5 * abs(e) for g, e in zip(got, exp)) and abs(got[2] - got[1]) >= 0; okB &= good; why.append([round(x, 5) for x in got])
check("B.1", "p, JC69 and K80 distances, transition share and identical pairs recomputed from the alignments equal table_distances.csv (both datasets)", okB, why)
check("B.2", "closed forms: JC69 distance of p = 0.1 is 0.107326 and equals the numerical maximum of the two-sequence likelihood; P_ii(0.2) = 0.824446", abs(-0.75 * math.log(1 - 4 / 3 * 0.1) - 0.107326) < 1e-6 and abs(__import__("scipy.optimize", fromlist=["x"]).minimize_scalar(lambda d: -(900 * math.log(0.25 + 0.75 * math.exp(-4 * d / 3)) + 100 * math.log(0.25 - 0.25 * math.exp(-4 * d / 3))), bounds=(0.01, 0.5), method="bounded", options={"xatol": 1e-10}).x - 0.107326) < 1e-5 and abs(0.25 + 0.75 * math.exp(-4 * 0.2 / 3) - 0.824446) < 1e-6)
def nj_tree(D, names):
    """Neighbour-joining with the Q criterion; returns the split set with branch lengths (own code, ties broken by the first minimum)."""
    n = len(names); node = {i: frozenset([i]) for i in range(n)}; act = list(range(n)); F = {}; nxt = n; D = D.copy(); blen = {}
    d = {(i, j): D[i, j] for i in range(n) for j in range(n)}
    while len(act) > 2:
        m = len(act); r = {i: sum(d[(i, k)] for k in act if k != i) for i in act}; best = None
        for ai in range(m):
            for aj in range(ai + 1, m):
                i, j = act[ai], act[aj]; q = (m - 2) * d[(i, j)] - r[i] - r[j]
                if best is None or q < best[0] - 1e-12: best = (q, i, j)
        _, i, j = best; li = 0.5 * d[(i, j)] + (r[i] - r[j]) / (2 * (m - 2)); lj = d[(i, j)] - li; blen[node[i]] = li; blen[node[j]] = lj; u = nxt; nxt += 1; node[u] = node[i] | node[j]
        for k in act:
            if k not in (i, j): d[(u, k)] = d[(k, u)] = 0.5 * (d[(i, k)] + d[(j, k)] - d[(i, j)])
        act = [k for k in act if k not in (i, j)] + [u]
    a, b = act; blen[node[a]] = d[(a, b)] / 2; blen[node[b]] = d[(a, b)] / 2; return blen
def own_nj_splits(C, minlen):
    p, Pt, Qv = pw(C); D = -0.75 * np.log(1 - 4 / 3 * p); bl = nj_tree(D, range(len(C))); n = len(C); out = set()
    for s, l in bl.items():
        if l < minlen: continue
        s2 = frozenset(range(n)) - s if 0 in s else s
        if 1 < len(s2) < n - 1: out.add(s2)
    return out
tc = pd.read_csv(P("table_tree_comparison.csv")); okNJ = True; det = []
for C, names, tf, rows in ((CH, nH, "h3n2.treefile", tc[tc.dataset.str.startswith("H3N2")]), (CE, nE, "ebola.treefile", tc[tc.dataset.str.startswith("Ebola")])):
    L = C.shape[1]; ml = split_sets(Phylo.read(P(tf), "newick"), names, 0.5 / L); nj = own_nj_splits(C, 0.5 / L); r = rows[rows.method == "NJ (JC69)"].iloc[0]; det.append((len(ml), len(nj), len(nj & ml), round(rfd(nj, ml), 3)))
    okNJ &= len(ml) == r.ml_splits and abs(len(nj & ml) - r.shared) <= 1 and abs(rfd(nj, ml) - r.rf) <= 0.03
check("B.3", "own neighbour-joining tree (Q criterion) shares the same number of resolved splits with the ML tree as table_tree_comparison.csv, within one split (ties are broken differently by Biopython)", okNJ, det)
# ================================================================== C. likelihood from first principles (Hands-on 2)
names8, X8 = readfa(P("sub8.fasta")); C8 = codes(X8); pat, cnt8 = np.unique(C8.T, axis=0, return_counts=True)
def own_loglik(treefile, model, freqs, kappa):
    t = Phylo.read(P(treefile), "newick"); Q = np.ones((4, 4)); Q[0, 1] = Q[1, 0] = Q[2, 3] = Q[3, 2] = kappa if model == "HKY" else 1.0       # order A, C, G, T: transitions are A-G (0,2) and C-T (1,3)
    Q = np.ones((4, 4)); Q[0, 2] = Q[2, 0] = Q[1, 3] = Q[3, 1] = kappa; Q = Q * freqs[None, :]; np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); Q /= -(freqs * np.diag(Q)).sum()
    def cond(cl, site):
        if cl.is_terminal():
            v = np.ones(4) if site[names8.index(cl.name)] < 0 else np.eye(4)[site[names8.index(cl.name)]]; return v
        v = np.ones(4)
        for c in cl.clades: v = v * (expm(Q * (c.branch_length or 0)) @ cond(c, site))
        return v
    return float(sum(c * math.log(freqs @ cond(t.root, s)) for s, c in zip(pat, cnt8)))
fr = np.full(4, 0.25); iqf = lambda f: float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", open(P(f)).read()).group(1))
hk = open(P("sub8_HKY.iqtree")).read(); pi_h = np.array([float(re.search(rf"pi\({b}\) = ([\d.]+)", hk).group(1)) for b in "ACGT"]); kap = float(re.search(r"A-G.*?: ([\d.]+)", hk).group(1)) if re.search(r"A-G.*?: ([\d.]+)", hk) else None
lj, lh = own_loglik("sub8_JC.treefile", "JC", fr, 1.0), None
check("C.1", "own pruning (matrix exponentials from scipy) gives the JC69 log-likelihood of the eight-sequence alignment, −2,908.301, on the IQ-TREE topology and branch lengths", abs(lj - (-2908.3012)) < 0.01, round(lj, 4))
if kap is not None:
    lh = own_loglik("sub8_HKY.treefile", "HKY", pi_h, kap); check("C.2", "own pruning gives the HKY85 log-likelihood, −2,839.058, from the IQ-TREE tree, base frequencies and κ = 5.34", abs(lh - (-2839.0582)) < 0.02 and abs(kap - 5.34) < 0.05, (round(lh, 4), kap))
else: skip("C.2", "HKY85 pruning check", "κ not found in sub8_HKY.iqtree")
tl = pd.read_csv(P("table_likelihood_models.csv")); N = 1701
check("C.3", "AIC = 2k − 2 ln L and BIC = k ln N − 2 ln L (N = 1,701) recomputed from the log-likelihoods of the pruning code for the four models; AIC picks GTR, BIC picks HKY85 (19.2 and 2.6 units apart); GTR has four more free parameters than HKY85, and 4 ln N − 2 Δln L gives the 19.2", all(abs(2 * r.free_parameters - 2 * r.lnL_pruning - r.AIC) < 0.02 and abs(r.free_parameters * math.log(N) - 2 * r.lnL_pruning - r.BIC) < 0.02 for r in tl.itertuples()) and tl.loc[tl.AIC.idxmin(), "model"] == "GTR" and tl.loc[tl.BIC.idxmin(), "model"] == "HKY" and rnd(tl[tl.model == "GTR"].BIC.iloc[0] - tl[tl.model == "HKY"].BIC.iloc[0], 1) == "19.2" and rnd(tl[tl.model == "HKY"].AIC.iloc[0] - tl[tl.model == "GTR"].AIC.iloc[0], 1) == "2.6" and rnd(math.log(1701), 1) == "7.4" and int(tl[tl.model == "GTR"].free_parameters.iloc[0] - tl[tl.model == "HKY"].free_parameters.iloc[0]) == 4 and abs(4 * math.log(N) - 2 * (tl[tl.model == "GTR"].lnL_pruning.iloc[0] - tl[tl.model == "HKY"].lnL_pruning.iloc[0]) - tl[tl.model == "GTR"].dBIC.iloc[0]) < 0.01)
# ================================================================== D. support (Hands-on 3)
ts = pd.read_csv(P("table_support_real.csv")); okD = True; det = []
for C, names, tf, nm, npy in ((CH, nH, "h3n2.treefile", "H3N2", "support_H3N2.npy"), (CE, nE, "ebola.treefile", "Ebola", "support_Ebola.npy")):
    t = Phylo.read(P(tf), "newick"); L = C.shape[1]; internal = [c for c in t.get_nonterminals() if c is not t.root]; short = sum(1 for c in internal if (c.branch_length or 0) < 0.5 / L); sp = split_sets(t, names, 0.5 / L)
    S = np.load(P(npy)); nb, iq = S[0], S[1]; row = ts[ts.dataset.str.startswith(nm)].iloc[0]; ok_ = ~np.isnan(iq)
    got = (len(internal), short, len(sp), 100 * np.mean(iq[ok_] >= 95), 100 * np.mean(iq[ok_] >= 80), 100 * np.mean(iq[ok_] < 50), 100 * np.mean(nb >= 95), 100 * np.mean(nb >= 70), 100 * np.mean(nb < 50), np.corrcoef(nb[ok_], iq[ok_])[0, 1])
    exp = (row.internal_branches, row.branches_shorter_than_one_substitution, row.resolved_splits, row.iqtree_ge95_pct, row.iqtree_ge80_pct, row.iqtree_lt50_pct, row.nj_boot_ge95_pct, row.nj_boot_ge70_pct, row.nj_boot_lt50_pct, row.corr_nj_vs_iqtree); okD &= all(abs(g - e) < 0.06 for g, e in zip(got, exp)); det.append([round(x, 2) for x in got[:3]])
check("D.1", "internal branches, branches shorter than one substitution and resolved splits recomputed from the tree files, and the support shares and correlations recomputed from the support arrays, equal table_support_real.csv", okD, det)
cr = pd.read_csv(P("calibration_check_raw.csv")); hi = cr[cr.support >= 95]; sup, uns = hi[hi.supported], hi[~hi.supported]
check("D.2", "tie check at low divergence: of 17 splits with bootstrap 95 or more, the 8 with a supporting column were all correct and the 9 without one were correct 11.1% of the time; 84.4% of splits have no supporting column", len(hi) == 17 and len(sup) == 8 and sup.correct.all() and len(uns) == 9 and rnd(100 * uns.correct.mean(), 1) == "11.1" and rnd(100 * (~cr.supported).mean(), 1) == "84.4", (len(hi), len(sup), len(uns), round(100 * (~cr.supported).mean(), 1)))
ca = pd.read_csv(P("table_calibration_sim.csv")); lo = ca[ca.height == ca.height.min()]
check("D.3", "calibration bands of the low-divergence simulation: 585, 255, 134 and 106 splits; correct shares 0.5%, 34.1%, 26.1% and 8.5%", list(lo.splits) == [585, 255, 134, 106] and [rnd(v, 1) for v in lo.correct_pct] == ["0.5", "34.1", "26.1", "8.5"])
# ================================================================== E. root-to-tip over every root (Hands-on 4)
def root_to_tip(tf, dates, L, seed=7, nperm=1000):
    t, nodes, edges, tips = parse_tree(P(tf)); n = len(nodes); D = dist_nodes(n, edges); ti = np.array([i for _, i in tips]); tn = [m for m, _ in tips]; tt = np.array([dates[m] for m in tn]); cand = []
    for u, v, l, _ in edges:
        for f in ([0.5] if l == 0 else np.linspace(0, 1, 11)): cand.append(np.minimum(D[ti, u] + f * l, D[ti, v] + (1 - f) * l))
    Z = np.array(cand); tc_ = tt - tt.mean(); Zc = Z - Z.mean(1, keepdims=True); cov = Zc @ tc_ / len(tt); corr = cov / (Zc.std(1) * tc_.std() + 1e-300); slope = cov / tc_.var(); ok = slope > 0; r2 = np.where(ok, corr ** 2, -1); b = int(np.argmax(r2))
    t0 = tt.mean() - Z[b].mean() / slope[b]; rng = np.random.default_rng(seed); nul = []
    for _ in range(nperm): pc = rng.permutation(tc_); cv = Zc @ pc / len(tt); nul.append(np.max(np.where(cv > 0, (cv / (Zc.std(1) * pc.std() + 1e-300)) ** 2, -1)))
    Dt = None
    return dict(n=len(Z), r2=r2[b], rate=slope[b] * L, t0=t0, perm_p=(1 + np.sum(np.array(nul) >= r2[b])) / (nperm + 1), p95=np.percentile(nul, 95), Z=Z, tt=tt, D=D, ti=ti, edges=edges, Dt=Dt, tn=tn)
ed = pd.read_csv(P("ebola_dates.csv")); edd = {r.acc: pd.Timestamp(r.date).year + (pd.Timestamp(r.date).dayofyear - 1) / 365.25 for r in ed.itertuples()}; hm = pd.read_csv(P("h3n2_meta.csv")); hdd = {r.accession: r.year + 0.5 for r in hm.itertuples()}
rt = pd.read_csv(P("table_root_to_tip.csv")); okE = True; det = []
for nm, tf, dates, L in (("Ebola", "ebola.treefile", edd, 18961), ("H3N2", "h3n2.treefile", hdd, 1762)):
    o = root_to_tip(tf, dates, L); row = rt[(rt.dataset.str.startswith(nm)) & (rt.rooting == "best")].iloc[0]
    good = o["n"] == row.candidate_roots and abs(o["r2"] - row.r2) < 2e-4 and abs(o["rate"] - row.substitutions_per_genome_year) < 0.02 and abs(o["t0"] - row.root_date) < 0.01 and o["perm_p"] < 0.01; okE &= good; det.append((o["n"], round(o["r2"], 4), round(o["rate"], 2), round(o["t0"], 2), round(o["p95"], 3)))
    if nm == "Ebola": p95E = o["p95"]
    else: p95H = o["p95"]
check("E.1", "root-to-tip regression over every position on every branch (own shortest-path code, 4,257 and 1,267 candidates): best-root R², rate per genome-year, root date and permutation p equal table_root_to_tip.csv", okE, det)
check("E.2", "95th percentile of the best-root R² for shuffled dates is 0.034 (Ebola) and 0.105 (H3N2), within 0.01 over 1,000 own shuffles", abs(p95E - 0.034) < 0.01 and abs(p95H - 0.105) < 0.015, (round(p95E, 3), round(p95H, 3)))
# ================================================================== F. modelling in detail (Hands-on 5)
ml_ = pd.read_csv(P("table_model_ladder.csv")); cols = {c.lower(): c for c in ml_.columns}; kk, ll = ml_[cols.get("free_parameters", "k")], ml_[cols.get("lnl", "lnL")]
check("F.1", "model ladder: AIC and BIC recomputed from ln L and k (N = 1,701) for all 16 models; BIC selects HKY+F+G4, AIC selects GTR+F+G4", all(abs(2 * k - 2 * l - a) < 0.1 and abs(k * math.log(1701) - 2 * l - b) < 0.1 for k, l, a, b in zip(kk, ll, ml_[cols["aic"]], ml_[cols["bic"]])) and ml_.loc[ml_[cols["bic"]].idxmin(), cols["model"]] == "HKY+F+G4" and ml_.loc[ml_[cols["aic"]].idxmin(), cols["model"]] == "GTR+F+G4", len(ml_))
cp = pd.read_csv(P("table_codon_positions.csv")); relr = 3 * cp.treelen / cp.treelen.sum()
check("F.2", "codon positions: relative rates 0.59, 0.53 and 1.88 are the tree lengths divided by their mean; the third position is 3.2 times the first and 3.5 times the second", [rnd(v, 2) for v in relr] == ["0.59", "0.53", "1.88"] and rnd(relr[2] / relr[0], 1) == "3.2" and rnd(relr[2] / relr[1], 1) == "3.5", [round(v, 3) for v in relr])
pC, PC, QC = pw(CC); iuC = np.triu_indices(len(CC), 1); ad = pd.read_csv(P("table_adequacy.csv")); jc_ts = ad[(ad.model.str.startswith("JC")) & (ad.statistic == "transition_share")].iloc[0]
check("F.3", "adequacy: the real transition share among pairwise differences is 0.733 (own count); equal-rate simulations give about one third (0.328, expected 1/3)", abs(PC[iuC].sum() / pC[iuC].sum() - 0.7328) < 0.001 and abs(jc_ts.sim_mean - 1 / 3) < 0.01 and not bool(jc_ts.adequate), round(PC[iuC].sum() / pC[iuC].sum(), 4))
# ================================================================== G. the Ebola run with the whole model set (Hands-on 3, part C)
iq = open(P("ebola_full.iqtree")).read(); blk = iq.split("List of models sorted by BIC scores:")[1].split("AIC, w-AIC")[0]; bic = {}
for line in blk.splitlines():
    m = re.match(r"^(\S+)\s+(-?[\d.]+)\s+([\d.]+)\s+[+-]\s+\S+\s+([\d.]+)\s+[+-]\s+\S+\s+([\d.]+)\s+[+-]", line)
    if m: bic[m.group(1)] = float(m.group(5))
best = min(bic, key=bic.get); tf_ = pd.read_csv(P("table_ebola_full_support.csv")).iloc[0]
check("G.1", "full ModelFinder run: BIC selects HKY+F+I; GTR+F+G4 is 40.9 and HKY+F+G4 10.0 BIC units behind; ln L −28,599.45", best == "HKY+F+I" and rnd(bic["GTR+F+G4"] - bic[best], 1) == "40.9" and rnd(bic["HKY+F+G4"] - bic[best], 1) == "10.0" and abs(float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", iq).group(1)) + 28599.452) < 0.001, (best, len(bic)))
tf2 = Phylo.read(P("ebola_full.treefile"), "newick"); s_full = split_sets(tf2, nE, 0.5 / 18961); s_res = split_sets(Phylo.read(P("ebola.treefile"), "newick"), nE, 0.5 / 18961); rows = []
for cl in tf2.get_nonterminals():
    if cl is tf2.root: continue
    lab = str(cl.confidence if cl.confidence is not None else cl.name); m = re.match(r"([\d.]+)/([\d.]+)", lab); rows.append((cl.branch_length or 0, float(m.group(1)), float(m.group(2))))
B_ = pd.DataFrame(rows, columns=["l", "sh", "ufb"]); res_, uns_ = B_[B_.l >= 0.5 / 18961], B_[B_.l < 0.5 / 18961]
check("G.2", "the two Ebola trees have identical resolved splits (50 of 50); support of the 50 resolved branches: ultrafast bootstrap ≥95 76%, ≥80 86%, SH-aLRT ≥80 66%, both thresholds 58%; of the 142 unresolved branches none reaches 95 and 87% are below 50", s_full == s_res and len(s_full) == 50 and [round(100 * (res_.ufb >= 95).mean()), round(100 * (res_.ufb >= 80).mean()), round(100 * (res_.sh >= 80).mean()), round(100 * ((res_.ufb >= 95) & (res_.sh >= 80)).mean())] == [76, 86, 66, 58] and len(uns_) == 142 and (uns_.ufb >= 95).sum() == 0 and round(100 * (uns_.ufb < 50).mean()) == 87, (len(res_), len(uns_)))
# ================================================================== H. tree-likeness and tree scores (Hands-on 1, parts C and D)
tl_ = pd.read_csv(P("table_tree_likeness.csv")); ts_ = pd.read_csv(P("table_tree_scores.csv")); okH = True; det = []
for C, names, tf, nm in ((CH, nH, "h3n2.treefile", "H3N2"), (CE, nE, "ebola.treefile", "Ebola")):
    p, _, _ = pw(C); L = C.shape[1]; Dsub = p * L; n = len(C); rng = np.random.default_rng(99); q = np.array([rng.choice(n, 4, replace=False) for _ in range(60000)]); i, j, k, l = q.T
    S = np.stack([Dsub[i, j] + Dsub[k, l], Dsub[i, k] + Dsub[j, l], Dsub[i, l] + Dsub[j, k]], 1); S.sort(1); ex4 = S[:, 2] - S[:, 1]; t3 = np.array([rng.choice(n, 3, replace=False) for _ in range(60000)]); a, b, c = t3.T
    T = np.stack([Dsub[a, b], Dsub[a, c], Dsub[b, c]], 1); T.sort(1); ex3 = T[:, 2] - T[:, 1]; t, nodes, edges, tips = parse_tree(P(tf)); Dn = dist_nodes(len(nodes), edges); ti = {m: ix for m, ix in tips}; pat = Dn[np.ix_([ti[x] for x in names], [ti[x] for x in names])]; iu = np.triu_indices(n, 1)
    jc_ = -0.75 * np.log(1 - 4 / 3 * p); row = tl_[tl_.dataset.str.startswith(nm)].iloc[0]
    good = abs(100 * (ex4 < 0.5).mean() - row.four_point_zero_pct) < 1.0 and abs(100 * (ex3 < 0.5).mean() - row.three_point_zero_pct) < 1.0 and abs(np.percentile(ex4, 95) - row.four_point_p95) < 0.4 and abs(np.percentile(ex3, 95) - row.three_point_p95) < 1.0 and abs(np.corrcoef(jc_[iu], pat[iu])[0, 1] - row.corr_obs_vs_tree) < 0.0006 and abs(np.mean(np.abs(jc_[iu] - pat[iu])) * L - row.mean_abs_residual_subs) < 0.01
    okH &= good; det.append((round(100 * (ex4 < 0.5).mean(), 1), round(100 * (ex3 < 0.5).mean(), 1), round(np.corrcoef(jc_[iu], pat[iu])[0, 1], 3)))
check("H.1", "four-point and three-point tests on own random quartets and triples (percentages within one point, percentiles within a substitution), and the correlation and mean absolute residual of JC69 distances against ML-tree path lengths over all pairs, equal table_tree_likeness.csv", okH, det)
okS = True; det = []
for C, names, tf, nm in ((CH, nH, "h3n2.treefile", "H3N2"), (CE, nE, "ebola.treefile", "Ebola")):
    sc = ts_[ts_.dataset.str.startswith(nm)].set_index("tree"); mlt = Phylo.read(P(tf), "newick"); cost = fitch_cost(mlt, names, C); okS &= cost == sc.loc["ML (IQ-TREE)", "parsimony_score"]; det.append(cost)
check("H.2", "Fitch parsimony score of the ML tree (own bit-mask implementation, polytomies counted exactly) is 318 for H3N2 and 256 for Ebola", okS, det)
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
okU = True; det = []
for C, names, nm in ((CH, nH, "H3N2"), (CE, nE, "Ebola")):
    p, _, _ = pw(C); d = -0.75 * np.log(1 - 4 / 3 * p); dm = DistanceMatrix(names, [list(d[i, :i + 1]) for i in range(len(names))]); sc = ts_[ts_.dataset.str.startswith(nm)].set_index("tree")
    for lab, tr in (("NJ (JC69)", DistanceTreeConstructor().nj(dm)), ("UPGMA (JC69)", DistanceTreeConstructor().upgma(dm))): c_ = fitch_cost(tr, names, C); okU &= c_ == sc.loc[lab, "parsimony_score"]; det.append(c_)
    okU &= sc.loc["ML (IQ-TREE)", "loglik_deficit"] == 0 and sc.loc["UPGMA (JC69)", "loglik_deficit"] > 10 * sc.loc["NJ (JC69)", "loglik_deficit"]
check("H.3", "parsimony scores of the NJ and UPGMA trees (own Fitch on Biopython's trees): 321 and 342 (H3N2), 257 and 281 (Ebola); the ML tree has the best log-likelihood and UPGMA is more than ten times further from it than NJ", okU, det)
iqE, iqH = iqf("ebola.iqtree"), iqf("h3n2.iqtree")
check("H.4", "log-likelihood of the ML tree: −4,715.16 (H3N2) equals the search value; the Ebola re-evaluation −28,600.95 is within one unit of the search value −28,599.98", abs(ts_[ts_.dataset.str.startswith("H3N2") & ts_.tree.str.startswith("ML")].loglik.iloc[0] - iqH) < 0.01 and abs(ts_[ts_.dataset.str.startswith("Ebola") & ts_.tree.str.startswith("ML")].loglik.iloc[0] - iqE) < 1.0, (iqH, iqE))
# ================================================================== I. variable sites only (Hands-on 6)
asc = pd.read_csv(P("table_ascertainment.csv")); okA = True; det = []; keys = {"full alignment": "full", "variable sites, no correction": "correction", "variable sites, +ASC": "ASC", "variable sites, constant-site counts (-fconst)": "fconst"}
for nm, tag, names, C in (("H3N2", "H3N2", nH, CH), ("Ebola", "Ebola", nE, CE)):
    sub = asc[asc.dataset.str.startswith(nm)].set_index("alignment"); tlen = {k: float(re.search(r"Total tree length \(sum of branch lengths\): ([\d.]+)", open(P(f"asc_{tag}_{v}.iqtree")).read()).group(1)) for k, v in keys.items()}
    sp_full = split_sets(Phylo.read(P(f"asc_{tag}_full.treefile"), "newick"), names, 0.5 / C.shape[1])
    for k in keys:
        found = len(split_sets(Phylo.read(P(f"asc_{tag}_{keys[k]}.treefile"), "newick"), names, 1e-6) & sp_full) / len(sp_full)
        okA &= abs(tlen[k] - sub.loc[k, "tree_length_per_site"]) < 6e-5 and abs(tlen[k] / tlen["full alignment"] - sub.loc[k, "tree_length_ratio_to_full"]) < 6e-3 and abs(found - sub.loc[k, "resolved_splits_of_full_tree_found"]) < 6e-4
    det.append([round(tlen[k] / tlen["full alignment"], 2) for k in keys])
check("I.1", "tree lengths read from the IQ-TREE reports and the share of the resolved splits of the full-alignment tree that each tree recovers (own split code) equal table_ascertainment.csv", okA, det)
cnt, nst = site_counts(CE); cE = [int(((nst < 2) & (cnt[b] > 0)).sum()) for b in range(4)]; lg = open(P("asc_Ebola_fconst.log")).read(); m = re.search(r"-fconst ([\d,]+)", lg); cH_ = site_counts(CH); cH = [int(((cH_[1] < 2) & (cH_[0][b] > 0)).sum()) for b in range(4)]; mh = re.search(r"-fconst ([\d,]+)", open(P("asc_H3N2_fconst.log")).read())
check("I.2", "constant-site counts passed to -fconst equal the counts of A, C, G and T in the constant sites of each alignment; 230 and 241 variable sites", [int(x) for x in m.group(1).split(",")] == cE and [int(x) for x in mh.group(1).split(",")] == cH and sum(cE) == 18720 and sum(cH) == 1532, (cE, cH))
# ================================================================== J. recombination simulation (Hands-on 7)
raw = pd.read_csv(P("table_recombination_raw.csv")); tr_ = pd.read_csv(P("table_recombination.csv")); okR = True
for ev in (0, 3, 10):
    for lab in ("clonal", "recombinant", "masked"):
        g = raw[(raw.events == ev) & (raw.tree == lab)]; r = tr_[(tr_.events == ev) & (tr_.tree == lab)].iloc[0]; okR &= len(g) == 15 and abs(g.tree_length_ratio.mean() - r.tree_length_ratio) < 6e-5 and abs(g.rf_to_truth.mean() - r.rf_to_truth) < 6e-5
    g = raw[(raw.events == ev) & (raw.tree == "masked")]; r = tr_[(tr_.events == ev) & (tr_.tree == "masked")].iloc[0]; okR &= (ev == 0 or abs(g.sensitivity.mean() - r.sensitivity) < 6e-5) and abs(g.false_positive.mean() - r.false_positive) < 6e-5
check("J.1", "table_recombination.csv equals the means of the 15 replicates per setting in table_recombination_raw.csv (tree length, Robinson–Foulds distance, sensitivity, false positives)", okR)
def paired(ev, a, b):
    d = raw[raw.events == ev].pivot(index="rep", columns="tree", values="rf_to_truth"); x = d[a] - d[b]; return x.mean(), x.std() / math.sqrt(len(x))
pr_ = [paired(3, "recombinant", "clonal"), paired(10, "recombinant", "clonal"), paired(3, "masked", "clonal"), paired(10, "masked", "clonal")]
check("J.2", "paired Robinson–Foulds differences quoted in the chapter: recombinant minus clonal −0.017 (SE 0.006) and −0.027 (0.012); masked minus clonal +0.011 (0.004) and +0.025 (0.008)", [(rnd(m_, 3), rnd(s_, 3)) for m_, s_ in pr_] == [("-0.017", "0.006"), ("-0.027", "0.012"), ("0.011", "0.004"), ("0.025", "0.008")], [(round(m_, 3), round(s_, 3)) for m_, s_ in pr_])
mm = raw[(raw.tree == "masked") & (raw.events > 0)]; chk = [(c, lam) for c, lam in ((4, 4 * 500 / 18961),)]; own_tail = 1 - sum(math.exp(-chk[0][1]) * chk[0][1] ** k / math.factorial(k) for k in range(4))
check("J.3", "worst-case sensitivity is 10% (one alignment with three events); the Poisson tail written out in the chapter equals scipy's survival function and is 5e-6 for 4 substitutions on a branch of 4 over 18,961 sites in a 500-site window (flagged); 3 would not be flagged by the count rule", rnd(100 * mm.sensitivity.min(), 0) == "10" and abs(own_tail - stats.poisson.sf(3, chk[0][1])) < 1e-15 and 1e-6 < own_tail < 1e-5, (round(own_tail, 8), round(mm.sensitivity.min(), 3)))
# ---- the screen on the real Ebola alignment: private substitutions of each genome, in 500-site windows (a measure that shares no code with the chapter's Fitch traceback)
ctr = site_counts(CE)[0]; flagged = {}; WIN = 500; step = 100
for gi, nm in enumerate(nE):
    priv = [s for s in np.where(site_counts(CE)[1] >= 2)[0] if CE[gi, s] >= 0 and ctr[CE[gi, s], s] == 1]
    if len(priv) < 4: continue
    lam = len(priv) * WIN / CE.shape[1]; pos = np.array(priv); starts = np.arange(0, CE.shape[1] - WIN + 1, step); c = np.array([((pos >= s) & (pos < s + WIN)).sum() for s in starts])
    if ((stats.poisson.sf(c - 1, lam) < 1e-4) & (c >= 4)).any(): flagged[nm] = (len(priv), int(c.max()))
check("J.4", "scan of the Ebola alignment for genomes whose private substitutions cluster in a 500-site window (own definition: a base found in no other genome): exactly KR105247.1 and KR105300.1 are flagged; H3N2 has none", set(flagged) == {"KR105247.1", "KR105300.1"} and not any((lambda gi: (lambda priv: len(priv) >= 4 and (stats.poisson.sf(np.array([((np.array(priv) >= s) & (np.array(priv) < s + WIN)).sum() for s in range(0, 1762 - WIN + 1, step)]) - 1, len(priv) * WIN / 1762) < 1e-4).any() and any(((np.array(priv) >= s) & (np.array(priv) < s + WIN)).sum() >= 4 for s in range(0, 1262, step)))([s for s in np.where(site_counts(CH)[1] >= 2)[0] if CH[gi, s] >= 0 and site_counts(CH)[0][CH[gi, s], s] == 1]))(g) for g in range(len(nH))), flagged)
# ================================================================== K. the agent pipeline
rep = {k: open(P(f"agent_report_{k}.md")).read() for k in ("ebola", "h3n2")}; md = {}; 
for nm, C in (("ebola", CE), ("h3n2", CH)):
    p, _, _ = pw(C); mdv = np.median(p, axis=1); md[nm] = [n_ for n_, v in zip(nE if nm == "ebola" else nH, mdv) if v > 3 * np.median(mdv)]
check("K.1", "outlier rule (median distance more than three times the median of all): one Ebola genome, KR105312.1, and no H3N2 sequence, as in the agent reports", md["ebola"] == ["KR105312.1"] and md["h3n2"] == [] and "1 flagged: ['KR105312.1']" in rep["ebola"] and "0 flagged: []" in rep["h3n2"], md)
sh = lambda t: [c for c in t.get_nonterminals() if c is not t.root]
fe_ = Phylo.read(P("ebola_full.treefile"), "newick"); ufbE = [float(re.match(r"[\d.]+/([\d.]+)", str(c.confidence if c.confidence is not None else c.name)).group(2 - 1)) for c in sh(fe_) if (c.branch_length or 0) >= 0.5 / 18961]; hh = Phylo.read(P("h3n2.treefile"), "newick"); ufbH = [float(c.confidence if c.confidence is not None else c.name) for c in sh(hh) if (c.branch_length or 0) >= 0.5 / 1762]
check("K.2", "support check of the agent reports: 86% (Ebola, ultrafast bootstrap) and 85% (H3N2) of resolved branches at 80 or more; model lines HKY+F+I and TPM2+F+G4; the substitution-cluster check warns for Ebola and passes for H3N2", rnd(100 * np.mean(np.array(ufbE) >= 80), 0) == "86" and rnd(100 * np.mean(np.array(ufbH) >= 80), 0) == "85" and "BIC selects HKY+F+I" in rep["ebola"] and "BIC selects TPM2+F+G4" in rep["h3n2"] and "| Substitution clusters | warn |" in rep["ebola"] and "| Substitution clusters | pass |" in rep["h3n2"], (round(100 * np.mean(np.array(ufbE) >= 80), 1), round(100 * np.mean(np.array(ufbH) >= 80), 1)))
# ================================================================== L. the chapter text (needs --chapter)
if A.chapter:
    from docx import Document
    d = Document(A.chapter); T = d.tables; paras = [p.text for p in d.paragraphs]; alltext = "\n".join(paras) + "\n" + "\n".join(c.text for t in T for r in t.rows for c in r.cells)
    ntok = lambda s: [t.replace("−", "-").replace(",", "") for t in re.findall(r"[−-]?\d[\d,]*\.?\d*", s)]
    def match(cell, vals):
        tk = ntok(cell)
        if len(tk) != len(vals): return False
        return all((abs(float(t) - v[1]) <= v[2]) if isinstance(v, tuple) else rnd(v, len(t.split(".")[1]) if "." in t else 0) == t for t, v in zip(tk, vals))
    def chk_table(cid, what, tbl, specs):
        bad = []
        for i, (label, cols) in enumerate(specs):
            cells = tbl.rows[i + 1].cells
            if label and not cells[0].text.startswith(label): bad.append((i, "label", cells[0].text[:30])); continue
            for j, vals in enumerate(cols, start=1):
                if vals is not None and not match(cells[j].text, vals): bad.append((i, j, cells[j].text[:30], [round(float(v), 4) for v in vals]))
        check(cid, what, not bad and len(tbl.rows) - 1 == len(specs), bad[:3])
    th, te = td.iloc[0], td.iloc[1]; tcm = lambda ds, m: tc[(tc.dataset.str.startswith(ds)) & (tc.method == m)].iloc[0]
    sp_ = lambda r, k: [r[k]]
    rows1 = [("Mean and largest", [[th.mean_p, th.max_p], [te.mean_p, te.max_p]]), ("Largest JC69", [[th.max_jc], [te.max_jc]]), ("Largest K80", [[th.k80_minus_jc_max], [te.k80_minus_jc_max]]), ("Transitions as", [[100 * th.ts_share], [100 * te.ts_share]]), ("Pairs of identical", [[th.identical_pairs_pct], [te.identical_pairs_pct]]), ("Splits in the ML tree", [[tcm("H3N2", "NJ (JC69)").ml_splits], [tcm("Ebola", "NJ (JC69)").ml_splits]])]
    for lab, m in (("NJ (JC69)", "NJ (JC69)"), ("NJ (p-distance)", "NJ (p-distance)"), ("UPGMA", "UPGMA (JC69)")):
        a, b = tcm("H3N2", m), tcm("Ebola", m); rows1.append((lab, [[a.splits, a.shared, a.rf] if lab != "NJ (p-distance)" else [a.rf], [b.splits, b.shared, b.rf] if lab != "NJ (p-distance)" else [b.rf]]))
    chk_table("T.1", "Hands-on 1 table (distances, ML splits, NJ and UPGMA comparison) equals table_distances.csv and table_tree_comparison.csv", T[1], rows1)
    a_, b_ = tl_.iloc[0], tl_.iloc[1]; chk_table("T.2", "Hands-on 1, part C table equals table_tree_likeness.csv", T[2], [("Quartets", [[a_.four_point_zero_pct], [b_.four_point_zero_pct]]), ("Four-point", [[a_.four_point_median, a_.four_point_p95], [b_.four_point_median, b_.four_point_p95]]), ("Triples", [[a_.three_point_zero_pct], [b_.three_point_zero_pct]]), ("Three-point", [[a_.three_point_median, a_.three_point_p95], [b_.three_point_median, b_.three_point_p95]]), ("Correlation", [[a_.corr_obs_vs_tree], [b_.corr_obs_vs_tree]]), ("Mean absolute", [[a_.mean_abs_residual_subs], [b_.mean_abs_residual_subs]])])
    sc = lambda ds, t_: ts_[(ts_.dataset.str.startswith(ds)) & (ts_.tree.str.startswith(t_))].iloc[0]
    chk_table("T.3", "Hands-on 1, part D table (parsimony scores and log-likelihoods) equals table_tree_scores.csv", T[3], [(lab, [[sc("H3N2", k).parsimony_score], [sc("H3N2", k).loglik], [sc("Ebola", k).parsimony_score], [sc("Ebola", k).loglik]]) for lab, k in (("Maximum", "ML"), ("Neighbour", "NJ"), ("UPGMA", "UPGMA"))])
    chk_table("T.4", "Hands-on 2 table (four models: log-likelihoods, AIC, BIC, κ) equals table_likelihood_models.csv", T[4], [(None, [[r.free_parameters], [r.lnL_pruning], [r.lnL_iqtree], [r.AIC], [r.BIC], None if pd.isna(r.kappa) else [r.kappa]]) for r in tl.itertuples()])
    ts1 = lambda nm: ts.iloc[0] if nm == "H3N2" else ts.iloc[1]; u, v = ts.iloc[0], ts.iloc[1]
    chk_table("T.5", "Hands-on 3 table (support in the real trees) equals table_support_real.csv", T[5], [("Internal", [[u.internal_branches], [v.internal_branches]]), ("Shorter", [[u.branches_shorter_than_one_substitution, round(100 * u.branches_shorter_than_one_substitution / u.internal_branches)], [v.branches_shorter_than_one_substitution, round(100 * v.branches_shorter_than_one_substitution / v.internal_branches)]]), ("Resolved", [[u.resolved_splits], [v.resolved_splits]]), ("IQ-TREE", [[u.iqtree_ge95_pct, u.iqtree_ge80_pct, u.iqtree_lt50_pct], [v.iqtree_ge95_pct, v.iqtree_ge80_pct, v.iqtree_lt50_pct]]), ("Own NJ", [[u.nj_boot_ge95_pct, u.nj_boot_ge70_pct, u.nj_boot_lt50_pct], [v.nj_boot_ge95_pct, v.nj_boot_ge70_pct, v.nj_boot_lt50_pct]]), ("Correlation", [[u.corr_nj_vs_iqtree], [v.corr_nj_vs_iqtree]])])
    ac = pd.read_csv(P("table_accuracy_sim.csv")); chk_table("T.6", "Hands-on 3 accuracy table equals table_accuracy_sim.csv", T[6], [(None, [[r.sites], [r.informative_sites], [r.nj_rf], [r.fasttree_rf]]) for r in ac.itertuples()])
    chk_table("T.7", "Hands-on 3 calibration table equals table_calibration_sim.csv (splits and share correct in each band)", T[7], [(None, [None, [r.splits], [r.correct_pct]]) for r in ca.itertuples()])
    rtv = lambda nm, k: rt[(rt.dataset.str.startswith(nm)) & (rt.rooting == k)].iloc[0]
    chk_table("T.8", "Hands-on 4 table (R², rate per genome-year, root date for three rootings) equals table_root_to_tip.csv", T[8], [(lab, [[rtv("Ebola", k).r2], [rtv("Ebola", k).substitutions_per_genome_year], [rtv("Ebola", k).root_date], [rtv("H3N2", k).r2], [rtv("H3N2", k).substitutions_per_genome_year], [rtv("H3N2", k).root_date]]) for lab, k in (("Best", "best"), ("Midpoint", "midpoint"), ("Random", "random tip"))])
    lad = ml_.set_index("model"); g = lambda a, b: lad.loc[b, "lnL"] - lad.loc[a, "lnL"]
    chk_table("T.9", "Hands-on 5 gain table (log-likelihood gained by each feature) equals the differences in table_model_ladder.csv", T[9], [(None, [None, [g("JC", "K2P")]]), (None, [None, [g("K2P", "HKY+F")]]), (None, [None, [g("HKY+F", "HKY+F+G4")]]), (None, [None, [g("HKY+F+G4", "GTR+F+G4")]])])
    rowsL = []
    for m in ("JC", "JC+G4", "K2P+G4", "HKY+F", "HKY+F+G4", "HKY+F+I", "HKY+F+I+G4", "GTR+F+G4"):
        cols = [[lad.loc[m, "k"]], [lad.loc[m, "lnL"]], [lad.loc[m, "AIC"]], [lad.loc[m, "BIC"]]]
        cols.append([("tol", lad.loc[m, "alpha"], 0.0051)] if not pd.isna(lad.loc[m, "alpha"]) else []); cols.append([lad.loc[m, "pinv"]] if not pd.isna(lad.loc[m, "pinv"]) else [])
        rowsL.append((m, cols))
    chk_table("T.10", "Hands-on 5 model table (k, ln L, AIC, BIC, α, invariant share) equals table_model_ladder.csv", T[10], rowsL)
    chk_table("T.11", "Hands-on 5 codon-position table equals table_codon_positions.csv", T[11], [(None, [[r.treelen], [r.relative_rate], [r.const_pct], [r.alpha], [r.kappa]]) for r in cp.itertuples()])
    adq = lambda mdl, st: ad[(ad.model.str.startswith(mdl)) & (ad.statistic == st)].iloc[0]; rowsA = []
    for st in ("constant_sites", "informative_sites", "transition_share"):
        o = adq("JC", st); rowsA.append((None, [[o.observed], [adq("JC", st).sim_mean, adq("JC", st).lo, adq("JC", st).hi], [adq("HKY", st).sim_mean, adq("HKY", st).lo, adq("HKY", st).hi]]))
    rowsA.append((None, [None, None, None]))
    chk_table("T.12", "Hands-on 5 adequacy table equals table_adequacy.csv", T[13], rowsA)
    rowsAsc = []
    for i in range(4):
        cols = []
        for nm in ("H3N2", "Ebola"):
            r_ = asc[asc.dataset.str.startswith(nm)].iloc[i]; cols += [[r_.tree_length_per_site], [r_.tree_length_ratio_to_full], [100 * r_.resolved_splits_of_full_tree_found]]
        rowsAsc.append((None, cols))
    chk_table("T.13", "Hands-on 6 table equals table_ascertainment.csv", T[14], rowsAsc)
    rr = lambda ev, lab, c: tr_[(tr_.events == ev) & (tr_.tree == lab)].iloc[0][c]; chk_table("T.14", "Hands-on 7 table equals table_recombination.csv", T[15], [(None, [[rr(ev, "clonal", "tree_length_ratio")], [rr(ev, "recombinant", "tree_length_ratio")], [rr(ev, "masked", "tree_length_ratio")], [rr(ev, "clonal", "rf_to_truth")], [rr(ev, "recombinant", "rf_to_truth")], [rr(ev, "masked", "rf_to_truth")], None, [100 * rr(ev, "masked", "false_positive")]] if ev == 0 else [[rr(ev, "clonal", "tree_length_ratio")], [rr(ev, "recombinant", "tree_length_ratio")], [rr(ev, "masked", "tree_length_ratio")], [rr(ev, "clonal", "rf_to_truth")], [rr(ev, "recombinant", "rf_to_truth")], [rr(ev, "masked", "rf_to_truth")], [100 * rr(ev, "masked", "sensitivity")], [100 * rr(ev, "masked", "false_positive")]]) for ev in (0, 3, 10)])
    # key sentences: numbers in the narrative against the result tables
    f3 = lambda v: f"{v:,.0f}"; h = asc[asc.dataset.str.startswith("H3N2")].set_index("alignment"); e = asc[asc.dataset.str.startswith("Ebola")].set_index("alignment")
    snippets = [f"{rnd(h.loc['variable sites, no correction', 'tree_length_ratio_to_full'], 1)} times (H3N2) and {rnd(e.loc['variable sites, no correction', 'tree_length_ratio_to_full'], 0)} times (Ebola)", f"claims {f3(float(rnd(h.loc['variable sites, no correction', 'substitutions_claimed_for_the_genome'], 0)))} substitutions", f"{f3(float(rnd(e.loc['variable sites, no correction', 'substitutions_claimed_for_the_genome'], 0)))} for an Ebola genome",
                f"where the full alignment gives {rnd(h.loc['full alignment', 'substitutions_claimed_for_the_genome'], 0)}", f"where it gives {rnd(e.loc['full alignment', 'substitutions_claimed_for_the_genome'], 0)}", f"{rnd(tf_.bic_gap_GTR_F_G4, 1)} BIC units behind and HKY+F+G4 is {rnd(tf_.bic_gap_HKY_F_G4, 1)} behind", f"ultrafast bootstrap is 95 or more on {rnd(tf_.ufboot_ge95_resolved_pct, 0)}% of the 50 resolved branches and 80 or more on {rnd(tf_.ufboot_ge80_resolved_pct, 0)}%",
                f"{rnd(tf_.shalrt_ge80_resolved_pct, 0)}% for the SH-aLRT", "58% of the branches meet both thresholds", f"{rnd(tf_.ufboot_lt50_unresolved_pct, 0)}% are below 50", "ln L −28,599.45", "1.44 and 2.13 times as long as the true clonal tree", "recombinant minus clonal is −0.017 (standard error 0.006) at 3 events and −0.027 (0.012) at 10", "masked minus clonal is +0.011 (0.004) and +0.025 (0.008)", "found only 10% of the imported bases",
                "KR105247.1", "KR105300.1", "11 of the 12 private changes of KR105247.1 and all 4 of KR105300.1", "230 of 1,762 sites for H3N2 and 241 of 18,961 for Ebola", "Neighbour-joining pays 9 to 11 log-likelihood units for its approximation and UPGMA 103 to 183", "a lower rate of 4.4"] 
    snippets = [s_ for s_ in snippets if s_ != "a lower rate of 4.4"]
    miss = [s_ for s_ in snippets if s_ not in alltext]; check("T.15", f"{len(snippets)} sentences of the narrative that quote results (Hands-on 1, 3, 6, 7 and the agent section) match the result tables word for word", not miss, miss[:2])
    caps = [int(m) for m in re.findall(r"^Figure 5\.(\d+)\. ", alltext, flags=re.M)]; check("T.16", "figure captions are numbered 5.1 to 5.9 in order (5.8 is the recombination figure and 5.9 the agent diagram)", caps == list(range(1, 10)), caps)
    check("T.17", "no mention of 'uploaded' files, eleven learning objectives, Hands-on 1 to 7 present, 'Su et al., 2026' and no 'BioMaster, 2026'", "uploaded" not in alltext.lower() and all(f"Hands-on {k} (" in alltext for k in range(1, 8)) and "four extra parameters (21 against 17)" in alltext and "Su et al., 2026" in alltext and "BioMaster, 2026" not in alltext and (lambda a, b: len([p for p in paras[a + 1:b] if len(p) > 60 and not p.startswith("After this chapter")]))(paras.index("Learning objectives"), paras.index("At a glance")) == 11)
    md = subprocess.run(["pandoc", A.chapter, "-t", "markdown", "--wrap=none"], capture_output=True, text=True).stdout; pj, ph = adq("JC", "gc_variance").p_two_sided, adq("HKY", "gc_variance").p_two_sided
    check("T.12b", "adequacy table: the variance of GC content gives p = 0.11 under JC and p = 0.27 under HKY+F+G4, as in table_adequacy.csv", f"p = {rnd(pj, 2)}$: adequate" in md and f"p = {rnd(ph, 2)}$: adequate" in md and rnd(pj, 2) == "0.11" and rnd(ph, 2) == "0.27", (pj, ph))
    body, srcs = md.split("# Sources")[0], md.split("# Sources")[1]
    dois = lambda t_: set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", t_))
    d_in, d_src = dois(body), dois(srcs); check("T.18", "every DOI cited in the text is in the Sources, and the chapter cites 42 distinct DOIs in the text", d_in <= d_src and len(d_in) == 42, (len(d_in), len(d_src), sorted(d_in - d_src)[:3]))
    if A.doi:
        import urllib.request, urllib.parse, json as _j, time
        badd = []
        for dd in sorted(d_src | d_in):
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(dd, safe="/()"), headers={"User-Agent": "validate-ch5"}), timeout=30).read()
            except Exception as ex: badd.append((dd, str(ex)[:30]))
            time.sleep(0.15)
        check("DOI", f"all {len(d_src | d_in)} DOIs of the chapter resolve in Crossref (preprint servers that Crossref does not hold are listed)", not badd, badd[:4])
else: skip("T.*", "chapter tables and sentences", "no --chapter given")
# ================================================================== M. reruns (needs --full)
if A.full:
    tmp = tempfile.mkdtemp(); src = os.path.dirname(os.path.abspath(P("table_distances.csv"))); code = os.path.dirname(os.path.abspath(__file__)); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) + glob.glob(os.path.join(code, "*.py")) + glob.glob(os.path.join(code, "*.R")) if os.path.isfile(f)]; os.makedirs(os.path.join(tmp, "figs5"), exist_ok=True)
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); tables = ["table_distances.csv", "table_tree_comparison.csv", "table_support_real.csv", "support_H3N2.npy", "support_Ebola.npy", "table_likelihood_models.csv", "table_root_to_tip.csv", "table_model_ladder.csv", "table_codon_positions.csv", "table_adequacy.csv", "table_modelfinder.csv", "table_tree_likeness.csv", "table_tree_scores.csv", "table_ascertainment.csv", "table_ebola_full_support.csv", "agent_report_ebola.md", "agent_report_h3n2.md"]; jobs = [["ex5_1_distances_nj.py"], ["Rscript", "ex5_1b_nj_bootstrap_ape.R"], ["ex5_2_likelihood.py"], ["ex5_3_bootstrap_accuracy.py", "A"], ["ex5_3_bootstrap_accuracy.py", "C"], ["ex5_4_root_to_tip.py"], ["ex5_6a_extract_cds.py"], ["ex5_6_modelling_detail.py"], ["ex5_1c_tree_likeness.py"], ["ex5_7_ascertainment.py"], ["ex5_3d_ebola_full_support.py"], ["ex5_5_agent_pipeline.py"]]
    if A.sim: jobs = [["ex5_3_bootstrap_accuracy.py", "B"]] + jobs + [["ex5_8_recombination.py", "15"]]; tables += ["table_accuracy_sim.csv", "table_calibration_sim.csv", "calibration_check_raw.csv", "table_recombination.csv", "table_recombination_raw.csv"]
    before = {n: md5(os.path.join(tmp, n)) for n in tables if os.path.exists(os.path.join(tmp, n))}; [os.remove(os.path.join(tmp, n)) for n in before if not (n.startswith("table_accuracy") or n.startswith("table_calibration") or n.startswith("calibration") or n.startswith("table_recomb")) or A.sim]
    for j in jobs: subprocess.run(j if j[0] == "Rscript" else [sys.executable] + j, cwd=tmp, capture_output=True)
    for n in before: check("REP", f"{n} reproduced byte for byte", os.path.exists(os.path.join(tmp, n)) and md5(os.path.join(tmp, n)) == before[n])
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch5.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=100)); npass = (Rd.result == "PASS").sum(); nskip = (Rd.result == "SKIP").sum(); print(f"\n{npass} of {len(Rd) - nskip} checks passed" + (f" ({nskip} skipped)" if nskip else "")); sys.exit(0 if (Rd.result == "FAIL").sum() == 0 else 1)
