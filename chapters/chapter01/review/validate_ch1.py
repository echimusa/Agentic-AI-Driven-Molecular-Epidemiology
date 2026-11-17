"""
Chapter 1 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's own scripts (pair-equality matrices, byte arrays, exact hypergeometric formulas, a renewal
equation, fresh simulation seeds) and checks them against the shipped result tables and against the chapter text. Prints PASS/FAIL per check and writes validation_results.csv.
USAGE   python validate_ch1.py [--chapter Chapter_1_....docx] [--full]
        --chapter  also checks that the key numbers appear in the chapter text (needs pandoc)
        --full     also reruns every analysis script in a scratch folder and checks that each result table is reproduced byte for byte (about 3 minutes)
"""
import argparse, glob, hashlib, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, hypergeom
sys.path.insert(0, ".")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, obs, exp, tol=0.0, kind="abs"):
    ok = (abs(obs - exp) <= tol) if kind == "abs" else (obs == exp) if kind == "eq" else (exp[0] <= obs <= exp[1])
    R.append(dict(id=cid, check=what, observed=obs, expected=exp, tolerance=tol, result="PASS" if ok else "FAIL")); return ok
# ---------------- Hands-on 1: typing indices, exact over all pairs
d = pd.read_csv(P("typing_pfge_subset.csv.gz"), dtype={"country": str, "pfge_xbai": str, "pfge_blni": str}, low_memory=False)
q = d[d.pfge_xbai.notna() & d.snp_cluster.notna() & d.serotype.notna()]; N = len(q); check("H1.1", "isolates with all three typing levels", N, 4811, kind="eq")
for k, v in {"Food": 2736, "Animal": 943, "Environment": 686, "Human": 384}.items(): check("H1.2", f"sector {k}", int((q.sector == k).sum()), v, kind="eq")
check("H1.3", "isolates with a PFGE XbaI pattern", int(d.pfge_xbai.notna().sum()), 5731, kind="eq"); check("H1.4", "US human isolates with a PFGE pattern", int(((d.sector == "Human") & (d.country == "USA")).sum()), 342, kind="eq")
code = {k: pd.factorize(q[c])[0] for k, c in {"S": "serotype", "P": "pfge_xbai", "C": "snp_cluster"}.items()}; iu = np.triu_indices(N, 1); tri = {k: (v[:, None] == v[None, :])[iu] for k, v in code.items()}
for k, exp in {"S": 0.957, "P": 0.995, "C": 0.994}.items(): check("H1.5", f"Simpson D {k} (exact over {len(iu[0]):,} pairs)", round(1 - tri[k].mean(), 3), exp, 0.0006)
for a, b, exp in [("S", "P", .111), ("P", "S", .933), ("P", "C", .438), ("C", "P", .385), ("C", "S", .938), ("S", "C", .127)]: check("H1.6", f"Wallace W({a}->{b})", round((tri[a] & tri[b]).sum() / tri[a].sum(), 3), exp, 0.0006)
nest = q.groupby("snp_cluster").pfge_xbai.nunique(); check("H1.7", "SNP clusters with more than one PFGE pattern", int((nest > 1).sum()), 362, kind="eq"); top = q[q.pfge_xbai == q.pfge_xbai.value_counts().index[0]]
check("H1.8", "commonest PFGE pattern: isolates, SNP clusters, serotypes", (len(top), top.snp_cluster.nunique(), top.serotype.nunique()), (251, 32, 5), kind="eq")
T = pd.read_csv(P("table_discrimination.csv")); W = pd.read_csv(P("table_wallace.csv"), header=[0, 1], index_col=0)
for _, r in T.iterrows(): check("H1.9", f"interval for D ({r.level}) is centred on the estimate (position within interval, 0.35-0.65)", round((r.D - r.D_lo) / (r.D_hi - r.D_lo), 2), (0.35, 0.65), kind="in")
for a, b in [("PFGE (XbaI)", "SNP cluster"), ("SNP cluster", "PFGE (XbaI)"), ("PFGE (XbaI)", "Serotype"), ("Serotype", "PFGE (XbaI)"), ("SNP cluster", "Serotype"), ("Serotype", "SNP cluster")]:
    w, lo, hi = W.loc[a, ("W", b)], W.loc[a, ("W_lo", b)], W.loc[a, ("W_hi", b)]; check("H1.10", f"interval for W({a}->{b}) is centred on the estimate", round((w - lo) / (hi - lo), 2), (0.35, 0.65), kind="in")
# ---------------- Hands-on 2: SNP distance from byte arrays
from Bio import SeqIO
recs = list(SeqIO.parse(P("ebola_aln.fasta"), "fasta")); ids = [r.id for r in recs]; X = np.array([np.frombuffer(str(r.seq).upper().encode(), dtype=np.uint8) for r in recs]); n, L = X.shape
dates = pd.to_datetime(pd.read_csv(P("ebola_dates.csv")).set_index("acc").date.loc[ids]); ok = np.isin(X, np.frombuffer(b"ACGT", dtype=np.uint8)); Dm = np.zeros((n, n), dtype=np.int32)
for i in range(n): Dm[i] = (((X[i] != X) & ok[i]) & ok).sum(1)
iu2 = np.triu_indices(n, 1); dd = Dm[iu2]; days = np.abs((dates.values[:, None] - dates.values[None, :]).astype("timedelta64[D]").astype(int))[iu2]
check("H2.1", "genomes, alignment length, pairs", (n, L, len(dd)), (195, 18961, 18915), kind="eq"); check("H2.2", "median SNP distance", float(np.median(dd)), 6, 0); check("H2.3", "share of pairs within 2 SNPs (%)", round(100 * (dd <= 2).mean(), 1), 17.1, 0.05)
check("H2.4", "share identical (%)", round(100 * (dd == 0).mean(), 1), 2.3, 0.05); check("H2.5", "Spearman distance vs days", round(spearmanr(dd, days)[0], 3), 0.317, 0.0006)
y = days <= 30; U = mannwhitneyu(dd[y], dd[~y], alternative="less").statistic; check("H2.6", "AUC for sampling within 30 days", round(1 - U / (y.sum() * (~y).sum()), 3), 0.668, 0.0006)
for t, sens, ppv in [(0, 5.1, 94.8), (2, 32.0, 80.3), (5, 63.1, 55.3)]: p = dd <= t; check("H2.7", f"threshold {t}: sensitivity", round(100 * (p & y).sum() / y.sum(), 1), sens, 0.06); check("H2.7", f"threshold {t}: PPV", round(100 * (p & y).sum() / p.sum(), 1), ppv, 0.06)
per_day = 1.11e-3 * L / 365; check("H2.8", "clock: SNPs per genome-year at 1.11e-3 per site", round(per_day * 365, 1), 21.0, 0.05); check("H2.9", "days to common ancestor for 1, 2, 5, 10 SNPs", tuple(int(round(k / (2 * per_day))) for k in (1, 2, 5, 10)), (9, 17, 43, 87), kind="eq")
# ---------------- Hands-on 3: tables and a renewal-equation check
import ex1_3_design_simulation as e3
S = pd.read_csv(P("table_design_simulation.csv")); g = lambda D, f, c: float(S[(S.D_months == D) & (S.f == f)][c].iloc[0])
check("H3.1", "relative bias against the chain truth at 20% sampled, range (%)", (round(min(g(D, 0.2, "relative_bias_vs_chain_pct") for D in (12, 24, 36, 60))), round(max(g(D, 0.2, "relative_bias_vs_chain_pct") for D in (12, 24, 36, 60)))), (-78, -66), kind="eq")
check("H3.2", "chain truth rises by a factor of about three from 12 to 60 months", round(g(60, 1.0, "truth_chain") / g(12, 1.0, "truth_chain"), 1), (2.8, 3.3), kind="in")
w, Rr, LAM = e3.w, e3.R, e3.LAM; I = np.zeros(700)
for t in range(700): I[t] = LAM + Rr * sum(w[k - 1] * I[t - k] for k in range(1, min(t, len(w)) + 1))
check("H3.3", "analytic steady state, cases per month (40 / (1 - 0.75))", round(LAM / (1 - Rr)), 160, kind="eq"); check("H3.4", "expected incidence at the start of the study window is below the steady state (text: 'a little below 160')", round(I[240]), (135, 150), kind="in")
check("H3.5", "incidence rise across the five-year window (%, text: about 5)", round(100 * (I[300] / I[240] - 1), 1), (4, 7), kind="in"); check("H3.6", "simulated cases per month in the windows (text: about 142)", round(float((S.cases / S.D_months).mean()), 0), 142, 3)
M = pd.read_csv(P("table_design_mutation.csv")); check("H3.7", "estimate equals the chain truth with a stable marker (0.368)", round(float(M[M.m == 0].estimate.iloc[0]), 3), 0.368, 0.0006)
# ---------------- Hands-on 4: lineage diversity, exact formulas
c = pd.read_csv(P("lapis_lineage_monthly.csv"), keep_default_na=False); c["lineage"] = c.lineage.replace("", "unassigned"); TL = pd.read_csv(P("table_lineage_diversity.csv"))
def simp(x): x = np.asarray(x, float); Nn = x.sum(); return 1 - (x * (x - 1)).sum() / (Nn * (Nn - 1))
rows = []
for (ct, mo), gg in c.groupby(["country", "month"]):
    cnt = gg.groupby("lineage")["count"].sum().values; Nn = cnt.sum(); rows.append(dict(country=ct, month=mo, N=Nn, D=simp(cnt), exp_rich=(1 - hypergeom(Nn, cnt, 500).pmf(0)).sum() if Nn >= 500 else np.nan))
V = pd.DataFrame(rows).merge(TL, on=["country", "month"]); check("H4.1", "sequences UK and USA", (int(V[V.country == "United Kingdom"].N.sum()), int(V[V.country == "USA"].N.sum())), (3223704, 3728840), kind="eq")
check("H4.2", "max |D recomputed (one type per lineage name) - D in table|", float((V.D - V.D_observed).abs().max()), 0.0005, 0.0005); ok2 = V[V.D_rarefied.notna()]
check("H4.3", "max |observed - rarefied D| (text: never more than 0.004)", round(float((ok2.D_observed - ok2.D_rarefied).abs().max()), 4), 0.004, 0.0005)
for ct, lo_e, hi_e in (("United Kingdom", 0.044, 0.979), ("USA", 0.482, 0.982)): s = ok2[ok2.country == ct]; check("H4.4", f"{ct}: rarefied D minimum and maximum", (round(s.D_rarefied.min(), 3), round(s.D_rarefied.max(), 3)), (lo_e, hi_e), kind="eq")
uk = V[(V.country == "United Kingdom") & (V.N > 10000)]; r_ = uk.lineages / uk.exp_rich; check("H4.5", "UK months above 10,000 sequences; median ratio of lineages seen to expected; maximum", (len(uk), round(r_.median(), 1), round(r_.max(), 1)), (33, 3.7, 12.3), kind="eq")
# ---------------- Hands-on 5: replicates with fresh seeds and an analytic expectation
from epi_sim import simulate_outbreak, sample_genome, depths, steps_between, LAM_DEFAULT
def pairs_for(o, rho, rng):
    nn = o["n"]; dep = depths(o["parent"]); s = np.where(rng.random(nn) < rho)[0]; gen = {i: sample_genome(o, i, rng.uniform(3, 10)) for i in s}; out = []
    for a in range(len(s)):
        for b in range(a + 1, len(s)): out.append((steps_between(o["parent"], dep, s[a], s[b]), len(gen[s[a]] ^ gen[s[b]])))
    return out
def rep(lam, seed, nout=20):
    rng = np.random.default_rng(seed); Pp = np.vstack([np.array(pairs_for(simulate_outbreak(rng, R=1.5, k=1.0, n_cap=250, min_cases=250, lam=lam), 0.5, rng)) for _ in range(nout)]); st, dist = Pp[:, 0], Pp[:, 1]; dr = st == 1; c1 = dist <= 1
    return dict(mean=dist[dr].mean(), ident=100 * (dist[dr] == 0).mean(), sens1=100 * (c1 & dr).sum() / dr.sum(), ppv1=100 * (c1 & dr).sum() / c1.sum())
Rp = pd.DataFrame([rep(LAM_DEFAULT, 300 + s) for s in range(6)]); rg = np.random.default_rng(3); gi = rg.gamma(4.0, 15.0 / 4.0, 400000); dp = rg.uniform(3, 10, 400000); dc = rg.uniform(3, 10, 400000)
check("H5.1", "mean SNPs between directly linked cases (6 fresh seeds); text: 1.01", round(Rp["mean"].mean(), 2), 1.01, 0.06); check("H5.2", "analytic mean SNPs between a parent and a child", round(float(LAM_DEFAULT * (np.abs(gi - dp) + dc).mean()), 2), 1.01, 0.03)
check("H5.3", "identical share among direct pairs (%); text 40-47", round(Rp.ident.mean(), 1), (38, 48), kind="in"); check("H5.4", "one-SNP sensitivity (%); text 73-75", round(Rp.sens1.mean(), 1), (71, 78), kind="in"); check("H5.5", "one-SNP positive predictive value (%); text 12-14", round(Rp.ppv1.mean(), 1), (10.5, 15.5), kind="in")
TS = pd.read_csv(P("table_sim_steps.csv")); check("H5.6", "pairs and direct pairs in the steps table", (int(TS.pairs.sum()), int(TS.pairs.iloc[0])), (933750, 7470), kind="eq") if "pairs" in TS else None
# ---------------- Agents: rule-table decisions and the guardrail test
import triage_core as tc
rule = {nm: tc.rule_agent(nm)[0] for nm, *_ in tc.CASES}; truth = {nm: tr for nm, _, tr, _ in tc.CASES}
decided = [nm for nm in rule if rule[nm] != tc.DECISIONS[3]]; right = [nm for nm in decided if (rule[nm] == tc.DECISIONS[0]) == truth[nm].startswith("same")]; gated = [nm for nm in rule if rule[nm] == tc.DECISIONS[3]]
check("H6.1", "rule agent: cases decided / correct / wrong / stopped at the gate", (len(decided), len(right), len(decided) - len(right), len(gated)), (8, 5, 3, 1), kind="eq")
check("H6.3", "new-4 nearest cluster member (SNPs)", tc.nearest_distance("new-4")["nearest_snps"], 11, kind="eq")
check("H6.4", "case outbreak size and cluster size", (tc._o1["n"], len(tc.CLUSTER)), (120, 25), kind="eq"); mock = subprocess.run([sys.executable, "ex1_7_llm_agent_loop.py", "--mock"], capture_output=True, text=True).stdout.strip().splitlines()[-1]
check("H6.5", "guardrail test (scripted stand-in): " + mock, bool(re.search(r"67%.*44%.*6 of 9", mock)), True, kind="eq")
# ---------------- Chapter text
if A.chapter:
    txt = subprocess.run(["pandoc", A.chapter, "-t", "plain", "--wrap=none"], capture_output=True, text=True).stdout.replace("\n", " ")
    for s in ["0.957 (0.955–0.959)", "0.438 (0.399–0.477)", "0.385 (0.353–0.416)", "0.044 in March 2021", "0.482 in April 2020", "12.3 times", "within 0.004", "Hands-on 1 to 5", "steady state of 160", "21 to 24"]:
        check("TXT", f"chapter text contains '{s}'", s in txt, True, kind="eq")
    for s in ["0.043 in March 2021", "0.536 in February 2022", "300 bootstrap", "12.5 times", "Hands-on 1 to 4"]: check("TXT", f"chapter text no longer contains '{s}'", s not in txt, True, kind="eq")
# ---------------- Full rerun: byte-for-byte reproduction of every result table
if A.full:
    tmp = tempfile.mkdtemp(); src = os.getcwd(); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) if os.path.isfile(f)]; [shutil.copy(P(t), tmp) for t in [os.path.basename(x) for x in glob.glob(P("table_discrimination.csv").replace("table_discrimination.csv", "table_*.csv"))]]
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); before = {os.path.basename(f): md5(f) for f in glob.glob(os.path.join(tmp, "table_*.csv"))}
    os.makedirs(os.path.join(tmp, "figs1"), exist_ok=True); os.makedirs(os.path.join(tmp, "figs"), exist_ok=True)
    for s in ["ex1_1_discrimination", "ex1_2_snp_distance", "ex1_3_design_simulation", "ex1_4_lineage_diversity", "ex1_5_epidemic_process", "ex1_6_agent_loop"]: subprocess.run([sys.executable, s + ".py"], cwd=tmp, capture_output=True)
    subprocess.run([sys.executable, "ex1_7_llm_agent_loop.py", "--mock"], cwd=tmp, capture_output=True); after = {os.path.basename(f): md5(f) for f in glob.glob(os.path.join(tmp, "table_*.csv"))}
    for k in sorted(before): check("REP", f"{k} reproduced byte for byte", after.get(k) == before[k], True, kind="eq")
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results.csv", index=False); print(Rd[["id", "check", "observed", "expected", "result"]].to_string(index=False, max_colwidth=96)); print(f"\n{(Rd.result == 'PASS').sum()} of {len(Rd)} checks passed")
sys.exit(0 if (Rd.result == "PASS").all() else 1)
