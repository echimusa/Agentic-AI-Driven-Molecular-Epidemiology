"""
Chapter 8 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts wherever that is possible.
A  the cached Zika pipeline: genomes, states, events, introductions, per-country events, the 75 balanced subsamples, the clock-rate grid and the dating summary
B  the ten reruns of the dating step: spread of the dates, counts of events and introductions
C  population-genetic summaries: pairwise distances by matrix products, nucleotide diversity, Dxy, Hudson's Fst, Tajima's D, the Mantel test
D  the simulated epidemics: every summary recomputed from the per-epidemic rows, and the internal consistency of those rows
E  the tuberculosis audit: every figure recomputed from the counts of the published paper
F  the Microreact files
G  the chapter text (tables, sentences, figure numbers, front matter, DOIs)
USAGE   python validate_ch8.py [--chapter Chapter_8_....docx] [--doi] [--full] [--sim]
        --chapter  also checks the chapter's tables, sentences, figure numbers and references against the results (needs pandoc and python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns the place reconstruction, the 75 subsamples, the population-genetic script, the audit, the export and two epidemics per scenario in a scratch folder and
                   checks that their result files are reproduced (needs treetime, biopython, shapely, matplotlib; about 5 minutes)
        --sim      with --full, also reruns the dating step once (about 100 s) and checks that it lands inside the spread of the ten reruns
Prints PASS/FAIL/SKIP per check, writes validation_results_ch8.csv and exits 1 if any check fails.
"""
import argparse, glob, hashlib, json, math, os, re, shutil, subprocess, sys, tempfile, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); ap.add_argument("--sim", action="store_true"); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
def P(name):
    for d in (".", "../data", "data", os.path.join(HERE, "../data")):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []; EXPECTED_DOIS = 35
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "../code"))
import phylo8
AM = phylo8.AMERICAS
# ================================================================== A. the cached Zika pipeline
meta = pd.read_csv(P("zika_selected_meta.csv"), dtype=str); outl = set(pd.read_csv(P("zika_clock_outliers.csv")).acc); m591 = meta[~meta.acc.isin(outl)].copy()
check("A.1", "600 sampled human genomes dated 2006 to 2018 from 40 countries; 9 clock outliers removed leaves 591", len(meta) == 600 and (meta.human == "True").all() and meta.year.astype(float).between(2006, 2018).all() and meta.cty.nunique() == 40 and len(outl) == 9 and len(m591) == 591, (len(meta), meta.cty.nunique(), len(m591)))
cnt = m591.cty.value_counts(); big = cnt[cnt >= 8].index; pooled = {phylo8.region_of(c) for c in cnt[cnt < 8].index}
check("A.2", "countries with at least 8 genomes are states of their own and the rest are pooled into three regional groups, 16 states in all; Brazil 96, Nicaragua 87, Puerto Rico 67, Singapore 62, Colombia 57", len(big) + len(pooled) == 16 and len(pooled) == 3 and [int(cnt[c]) for c in ("Brazil", "Nicaragua", "Puerto Rico", "Singapore", "Colombia")] == [96, 87, 67, 62, 57] and int(cnt[cnt < 8].max()) <= 7, (len(big), pooled))
ev = pd.read_csv(P("events.csv")); is_am = lambda x: x in AM or x == "Other Americas"; ev["intro"] = [(not is_am(a)) and is_am(b) for a, b in zip(ev.src, ev.dst)]; main = ev[ev.intro].sort_values("child_tips", ascending=False).iloc[0]
check("A.3", "92 inferred events, 4 introductions into the Americas; the main one runs from Other Asia-Pacific to Brazil on a branch from 2012.21 to 2012.44 and leads to 487 of 591 genomes", len(ev) == 92 and int(ev.intro.sum()) == 4 and main.src == "Other Asia-Pacific" and main.dst == "Brazil" and rnd(main.date_lo, 2) == "2012.21" and rnd(main.date_hi, 2) == "2012.44" and int(main.child_tips) == 487, (len(ev), int(ev.intro.sum())))
within = ev[ev.src.map(is_am) & ev.dst.map(is_am)]; FIRST = 2015 + 4 / 12
check("A.4", "68 of the 92 events are between American places and 34 of those are dated, at the middle of their branch, before the first confirmed detection in Brazil in May 2015", len(within) == 68 and int((within.date < FIRST).sum()) == 34, (len(within), int((within.date < FIRST).sum())))
out_ev = ev.groupby("src").size(); ss = pd.read_csv(P("state_summary.csv"))
check("A.5", "Nicaragua has 24 outgoing events in the full sample, Brazil and the Dominican Republic 16 each; the state summary agrees", out_ev["Nicaragua"] == 24 and out_ev["Brazil"] == 16 and out_ev["Dominican Republic"] == 16 and ss.set_index(ss.columns[0]).loc["Nicaragua"].iloc[1] == 24 if ss.shape[1] > 2 else out_ev["Nicaragua"] == 24, dict(out_ev.sort_values(ascending=False).head(3)))
b = pd.read_csv(P("bias_replicates.csv")); b["cap"] = b.cap.astype(str); nf = b[b.cap != "full"]; am_states = [s for s in b.state.unique() if s in AM]
def shares(g):
    out = {}
    for rep, gg in g[g.state.isin(am_states)].groupby("rep"):
        top = gg[gg.out_events == gg.out_events.max()].state.tolist()
        for t in top: out[t] = out.get(t, 0) + 1 / len(top)
    return {k: v / g.rep.nunique() for k, v in out.items()}
s10, s20, s40 = shares(nf[nf.cap == "10"]), shares(nf[nf.cap == "20"]), shares(nf[nf.cap == "40"]); mr = nf.drop_duplicates(["cap", "rep"])
check("A.6", "75 balanced subsamples (caps 10, 20, 40, 25 each): the main introduction is Asia-Pacific to Brazil in all 75 with a median date of 2012.33; median introductions 3, 4, 4; Nicaragua's median outgoing events 0, 0, 1; the top American exporter by cap as in the text", len(mr) == 75 and (mr.main_dst == "Brazil").all() and (mr.main_src == "Other Asia-Pacific").all() and rnd(mr.main_date.median(), 2) == "2012.33" and [mr[mr.cap == c].n_intro.median() for c in ("10", "20", "40")] == [3, 4, 4] and [nf[(nf.cap == c) & (nf.state == "Nicaragua")].out_events.median() for c in ("10", "20", "40")] == [0, 0, 1]
      and [rnd(s10.get(k, 0), 2) for k in ("Brazil", "Dominican Republic", "Honduras")] == ["0.64", "0.32", "0.04"] and [rnd(s20.get(k, 0), 2) for k in ("Dominican Republic", "Brazil")] == ["0.58", "0.42"] and [rnd(s40.get(k, 0), 2) for k in ("Dominican Republic", "Honduras", "Brazil")] == ["0.45", "0.29", "0.27"], (len(mr), s10, s20, s40))
allr = b[~b.state.str.startswith("Other")]; rho = stats.spearmanr(allr.tips, allr.out_events)[0]
check("A.7", "Spearman correlation between genomes sampled and outgoing events is 0.11 over 988 country-runs", len(allr) == 988 and rnd(rho, 2) == "0.11", (len(allr), rho))
rs = pd.read_csv(P("rate_sensitivity.csv")); ts = json.load(open(P("timetree_summary.json")))
check("A.8", "dating summary of the chapter run: clock rate 5.76e-4, root 1993.6 (1992.3 to 1994.9), Brazilian ancestor 2012.44 (2012.26 to 2012.62), its parent 2012.21 (2011.97 to 2012.43), 9 outliers; fixed rates 8e-4, 1.0e-3, 1.2e-3 move the introduction branch to 2013.3-2013.5, 2013.6-2013.9, 2013.7-2014.0 and the root to 1999.2, 2001.8, 2003.4",
      rnd(1e4 * ts["clock_rate"], 2) == "5.76" and rnd(ts["root_date"], 1) == "1993.6" and [rnd(x, 1) for x in ts["root_ci90"]] == ["1992.3", "1994.9"] and rnd(ts["brazil_mrca"], 2) == "2012.44" and [rnd(x, 2) for x in ts["brazil_mrca_ci90"]] == ["2012.26", "2012.62"] and rnd(ts["brazil_parent"], 2) == "2012.21" and [rnd(x, 2) for x in ts["brazil_parent_ci90"]] == ["2011.97", "2012.43"] and ts["n_outliers"] == 9
      and [(rnd(r.brazil_parent, 1), rnd(r.brazil_mrca, 1)) for r in rs[rs.setting == "fixed"].itertuples()] == [("2013.3", "2013.5"), ("2013.6", "2013.9"), ("2013.7", "2014.0")] and [rnd(x, 1) for x in rs[rs.setting == "fixed"].root] == ["1999.2", "2001.8", "2003.4"], rs.round(2).to_dict("records")[1:2])
# ================================================================== B. the reruns of the dating step
D = pd.read_csv(P("dating_reruns.csv")); DS = pd.read_csv(P("table_dating_reruns_summary.csv")).set_index("key")
ok = len(D) == 10 and all(abs(DS.loc[k, "min"] - D[k].min()) < 1e-6 and abs(DS.loc[k, "median"] - D[k].median()) < 1e-6 and abs(DS.loc[k, "max"] - D[k].max()) < 1e-6 for k in ("clock_rate", "root", "brazil_parent", "brazil_mrca", "n_events", "n_intro", "main_tips"))
check("B.1", "the summary table of the ten reruns equals the minimum, median and maximum of dating_reruns.csv", ok, len(D))
check("B.2", "ten reruns: clock rate 5.34 to 6.54 (10^-4), root 1992.1 to 1996.0, Brazilian ancestor 2012.17 to 2012.86 (range 0.69 years, against a 90% interval of 0.36 years for the chapter run), main introduction 487 genomes every time, 91 to 93 events, 4 introductions in nine reruns and 3 in one",
      [rnd(1e4 * D.clock_rate.min(), 2), rnd(1e4 * D.clock_rate.max(), 2)] == ["5.34", "6.54"] and [rnd(D.root.min(), 1), rnd(D.root.max(), 1)] == ["1992.1", "1996.0"] and [rnd(D.brazil_mrca.min(), 2), rnd(D.brazil_mrca.max(), 2)] == ["2012.17", "2012.86"] and rnd(D.brazil_mrca.max() - D.brazil_mrca.min(), 2) == "0.69" and rnd(ts["brazil_mrca_ci90"][1] - ts["brazil_mrca_ci90"][0], 2) == "0.36"
      and (D.main_tips == 487).all() and D.n_events.min() == 91 and D.n_events.max() == 93 and (D.n_intro == 4).sum() == 9 and (D.n_intro == 3).sum() == 1, D[["clock_rate", "brazil_mrca", "n_events", "n_intro"]].describe().loc[["min", "max"]].round(5).to_dict())
check("B.3", "faster estimated clocks give later dates across the reruns (rank correlation above 0.9); no rerun reaches May 2013, the start of the published estimate; the chapter run lies between the lowest and the highest rerun", stats.spearmanr(D.clock_rate, D.brazil_mrca)[0] > 0.9 and (D.main_hi < 2013 + 4 / 12).all() and D.brazil_mrca.min() < ts["brazil_mrca"] < D.brazil_mrca.max(), stats.spearmanr(D.clock_rate, D.brazil_mrca)[0])
# ================================================================== C. population-genetic summaries
from Bio import SeqIO
recs = [r for r in SeqIO.parse(P("zika_aln.fasta"), "fasta") if r.id in set(m591.acc)]; ids = [r.id for r in recs]; mm = m591.set_index("acc").loc[ids]; region = np.array([phylo8.region_of(c) for c in mm.cty]); cty = mm.cty.to_numpy()
codes = np.full((len(recs), len(recs[0].seq)), -1, np.int8)
for i, r in enumerate(recs):
    a = np.frombuffer(str(r.seq).upper().encode(), np.uint8)
    for k, ch in enumerate("ACGT"): codes[i, a == ord(ch)] = k
V = (codes >= 0).astype(np.float32); matches = sum(((codes == k).astype(np.float32)) @ ((codes == k).astype(np.float32)).T for k in range(4)); validp = V @ V.T       # matrix products instead of a loop over rows
Dm = (validp - matches) / np.maximum(validp, 1); np.fill_diagonal(Dm, np.nan)
def pi_of(ix): return np.nanmean(Dm[np.ix_(ix, ix)])
def fst_of(i1, i2): return 1 - (pi_of(i1) + pi_of(i2)) / 2 / Dm[np.ix_(i1, i2)].mean()
GP = pd.read_csv(P("table_popgen_groups.csv")).set_index("group"); RGt = pd.read_csv(P("table_popgen_regions.csv")); CT = pd.read_csv(P("table_popgen_countries.csv")); IB = pd.read_csv(P("table_popgen_ibd.csv")).iloc[0]
grp_ix = {g: np.where(region == g)[0] for g in ("Americas", "Asia-Pacific", "Europe/Africa")}; grp_ix.update({g: np.where(cty == g)[0] for g in GP.index if g not in grp_ix})
check("C.1", "genomes and nucleotide diversity (differences per kb, pairwise deletion) of the 3 regions and 10 countries recomputed with matrix products equal table_popgen_groups.csv", all(len(grp_ix[g]) == GP.loc[g, "n"] and abs(1000 * pi_of(grp_ix[g]) - GP.loc[g, "pi_per_kb"]) < 1e-3 for g in GP.index), {g: round(1000 * pi_of(grp_ix[g]), 3) for g in list(GP.index)[:3]})
def tajima(ix):
    sub = codes[ix]; ok_cols = (sub >= 0).all(axis=0); sub = sub[:, ok_cols]; m = len(ix)
    seg = np.array([len(np.unique(sub[:, j])) > 1 for j in range(sub.shape[1])]); S = int(seg.sum()); X = sub[:, seg]
    pw = np.mean([np.sum(X[i] != X[j]) for i in range(m) for j in range(i + 1, m)]) if m < 60 else (sum(np.sum(X[i] != X[i + 1:], axis=1).sum() for i in range(m - 1)) / (m * (m - 1) / 2))
    a1 = sum(1 / i for i in range(1, m)); a2 = sum(1 / i ** 2 for i in range(1, m)); b1 = (m + 1) / (3 * (m - 1)); b2 = 2 * (m * m + m + 3) / (9 * m * (m - 1)); c1 = b1 - 1 / a1; c2 = b2 - (m + 2) / (a1 * m) + a2 / a1 ** 2
    return (pw - S / a1) / math.sqrt(c1 / a1 * S + c2 / (a1 ** 2 + a2) * S * (S - 1)), S, int(ok_cols.sum())
colfrac = V.mean(axis=0); core = colfrac >= 0.90; good = V[:, core].mean(axis=1) >= 0.99; okT = True; det = []
for g in GP.index:
    ix = grp_ix[g][good[grp_ix[g]]]; d_, S_, sites_ = tajima(ix); okT &= len(ix) == GP.loc[g, "n_tajima"] and abs(d_ - GP.loc[g, "tajima_d"]) < 1e-3 and S_ == GP.loc[g, "segregating_sites"] and sites_ == GP.loc[g, "sites_used"]; det.append((g, round(d_, 3)))
check("C.2", "Tajima's D, segregating sites and sites used, recomputed with a second implementation on the 422 well-covered genomes, equal the table for all 13 groups; 422 of 591 genomes pass the coverage rule", okT and int(good.sum()) == 422 and int(core.sum()) == 9716, det[:3])
okF = True
for r in RGt.itertuples():
    f = fst_of(grp_ix[r.group1], grp_ix[r.group2]); okF &= abs(f - r.fst) < 1e-4 and abs(1000 * Dm[np.ix_(grp_ix[r.group1], grp_ix[r.group2])].mean() - r.dxy_per_kb) < 1e-3
    pool = np.concatenate([grp_ix[r.group1], grp_ix[r.group2]]); rng = np.random.default_rng(3); cnt = 0
    for _ in range(199):
        rng.shuffle(pool); cnt += fst_of(pool[:len(grp_ix[r.group1])], pool[len(grp_ix[r.group1]):]) >= f
    okF &= (cnt + 1) / 200 <= 0.02 and r.p_perm <= 0.01
check("C.3", "Dxy and Hudson's Fst of the three region pairs equal the table (Americas against Asia-Pacific 0.35, against Europe/Africa 0.09, Asia-Pacific against Europe/Africa 0.38) and a fresh permutation test with 199 permutations agrees that all three are significant", okF and [rnd(x, 2) for x in RGt.fst] == ["0.35", "0.09", "0.38"], RGt.fst.round(3).tolist())
cs = list(dict.fromkeys(list(CT.country1) + list(CT.country2))); okC = len(CT) == 45 and len(cs) == 10
for r in CT.itertuples(): okC &= abs(fst_of(grp_ix[r.country1], grp_ix[r.country2]) - r.fst) < 1e-4
check("C.4", "Fst for the 45 pairs of the 10 American countries with at least 12 genomes equals table_popgen_countries.csv", okC, len(CT))
r_, _ = stats.pearsonr(CT.fst_lin, CT.km); rho_ = stats.spearmanr(CT.fst_lin, CT.km)[0]; lab = {c: i for i, c in enumerate(cs)}; k = len(cs); Mg = np.zeros((k, k)); Mk = np.zeros((k, k))
for r in CT.itertuples(): Mg[lab[r.country1], lab[r.country2]] = Mg[lab[r.country2], lab[r.country1]] = r.fst_lin; Mk[lab[r.country1], lab[r.country2]] = Mk[lab[r.country2], lab[r.country1]] = r.km
iu = np.triu_indices(k, 1); rng = np.random.default_rng(11); pm = (1 + sum(np.corrcoef(Mg[np.ix_(p, p)][iu], Mk[iu])[0, 1] >= r_ for p in (rng.permutation(k) for _ in range(9999)))) / 10000
check("C.5", "Mantel test: Pearson r = -0.17 and Spearman rho = -0.10 between Fst/(1 - Fst) and distance recomputed from the pair table; a fresh 9,999-permutation p-value (0.77 in the table) agrees within 0.03; mean Fst 0.35, range 0.03 to 0.62", rnd(r_, 2) == "-0.17" and rnd(rho_, 2) == "-0.10" and abs(r_ - IB.pearson_r) < 1e-3 and abs(pm - IB.p_mantel_one_sided) < 0.03 and rnd(CT.fst.mean(), 2) == "0.35" and rnd(CT.fst.min(), 2) == "0.03" and rnd(CT.fst.max(), 2) == "0.62", (r_, rho_, pm))
# ================================================================== D. the simulated epidemics
T = pd.read_csv(P("sampling_truth_runs.csv")); TS = pd.read_csv(P("table_sampling_truth_summary.csv")).set_index(["scenario", "strategy"]); NAMES = ["A", "B", "C", "D"]
okI = len(T) == 480 and T.groupby(["scenario", "strategy"]).size().eq(60).all() and (T[[f"tips_{x}" for x in NAMES]].sum(axis=1) == T.n_tips).all() and (T[[f"true_out_{x}" for x in NAMES]].sum(axis=1) == T.true_events).all() and (T[[f"inf_out_{x}" for x in NAMES]].sum(axis=1) == T.inf_events).all()
check("D.1", "480 rows (2 scenarios, 60 epidemics, 4 strategies); in every row the tips, the true events and the inferred events add up over the four places", okI, len(T))
okS = True; summ = {}
for (sc, st), g in T.groupby(["scenario", "strategy"]):
    tr = g[[f"true_out_{x}" for x in NAMES]].to_numpy(float); inf = g[[f"inf_out_{x}" for x in NAMES]].to_numpy(float); st_ = tr[:, 3] / np.maximum(tr.sum(1), 1); si = inf[:, 3] / np.maximum(inf.sum(1), 1); bias = si - st_; n = len(g); t = TS.loc[(sc, st)]
    okS &= abs(bias.mean() - t.bias_D) < 1e-3 and abs(bias.mean() - 1.96 * bias.std(ddof=1) / n ** 0.5 - t.bias_D_lo) < 1e-3 and abs(bias.mean() + 1.96 * bias.std(ddof=1) / n ** 0.5 - t.bias_D_hi) < 1e-3 and abs(st_.mean() - t.true_share_D) < 1e-3 and abs(si.mean() - t.inferred_share_D) < 1e-3 and abs((g.true_root == g.inf_root).mean() - t.root_correct) < 1e-3 and abs((inf.argmax(1) == tr.argmax(1)).mean() - t.top_exporter_correct) < 1e-3 and abs(g.n_tips.mean() - t.tips) < 0.01
    summ[(sc, st)] = bias.mean()
check("D.2", "the bias in the share of exports attributed to place D, its 95% interval, the true and inferred shares, root accuracy and top-exporter accuracy, recomputed per scenario and strategy, equal table_sampling_truth_summary.csv", okS, {f"{a}/{b}": round(v, 3) for (a, b), v in list(summ.items())[:4]})
bb = lambda sc, st: 100 * summ[(sc, st)]
check("D.3", "uneven sampling inflates the over-sampled place (+7.8 and +6.6 points as sequenced), proportional sampling nearly removes it (+2.0 and +1.9), the bias correction makes it worse (+10.8 and +8.8), and capping at 15 does not remove it (+11.4 and +6.8)", [rnd(bb("hub", s), 1) for s in ("as_sequenced", "proportional", "as_sequenced_sbc", "cap_per_place")] == ["7.8", "2.0", "10.8", "11.4"] and [rnd(bb("no_hub", s), 1) for s in ("as_sequenced", "proportional", "as_sequenced_sbc", "cap_per_place")] == ["6.6", "1.9", "8.8", "6.8"], {k: round(v, 3) for k, v in summ.items()})
tt = T.set_index(["scenario", "strategy"])
check("D.4", "design: capped samples have at most 60 genomes, proportional samples are smaller than as-sequenced ones, the as-sequenced pool has a mean of 140 (hub) and 162 (no hub) hosts, and place D supplies 35% and 42% of the as-sequenced genomes", (T[T.strategy == "cap_per_place"].n_tips <= 60).all() and all(TS.loc[(sc, "proportional"), "tips"] < TS.loc[(sc, "as_sequenced"), "tips"] for sc in ("hub", "no_hub")) and [rnd(TS.loc[(sc, "as_sequenced"), "tips"], 0) for sc in ("hub", "no_hub")] == ["140", "162"] and [rnd(100 * TS.loc[(sc, "as_sequenced"), "share_tips_D"], 0) for sc in ("hub", "no_hub")] == ["35", "42"], TS.tips.round(0).to_dict())
# ================================================================== E. the tuberculosis audit
TOTAL = 3204; AUD = pd.read_csv(P("table_tb_audit.csv")); a = AUD.set_index("item").recomputed
exp = {"Lineage 2 share of all genomes (%)": 100 * 2745 / 3204, "Sublineage 2.2.1 share of all genomes (%)": 100 * 2591 / 3204, "MDR share of all genomes (%)": 100 * 1315 / 3204, "Pyrazinamide resistance (%)": 100 * 339 / 3204, "Ethambutol resistance (%)": 100 * 855 / 3204, "Streptomycin resistance (%)": 100 * 1101 / 3204, "Mean cluster size at 5 SNPs": 1384 / 471, "Mean cluster size at 10 SNPs": 1596 / 499, "Mean cluster size at 12 SNPs": 1631 / 504, "Clustered share, lineage 2 (%)": 100 * 1462 / 2745, "Clustered share, lineage 4 (%)": 100 * 132 / 443, "Share of all genomes from Shandong (%)": 100 * 1449 / 3204}
check("E.1", "the audit table recomputes from the paper's counts: 14 of 22 reported figures agree and 8 do not (MDR share 41.0% against 27.6%; pyrazinamide, ethambutol and streptomycin 10.6%, 26.7% and 34.4% against 12.8% each; mean cluster sizes 2.94, 3.20 and 3.24 against 3.05, 3.17 and 3.18; 80 MDR isolates in a sublineage of 3); Shandong supplied 45.2% of the genomes", all(abs(a[k] - v) < 1e-3 for k, v in exp.items()) and int(AUD[AUD.reported.notna()].agrees.sum()) == 14 and int((~AUD[AUD.reported.notna()].agrees).sum()) == 8 and len(AUD[AUD.reported.notna()]) == 22 and rnd(a["Share of all genomes from Shandong (%)"], 1) == "45.2", AUD[AUD.reported.notna()].agrees.value_counts().to_dict())
chi = stats.chi2_contingency([[1462, 1283], [132, 311]], correction=False)
check("E.2", "lineage 2 against lineage 4 clustering (1,462 of 2,745 against 132 of 443): chi-square 84.0 and p below 0.001, as the paper reports; the clustered strains add up (1,462 + 132 + 2 = 1,596) and so do the clusters (413 + 86 = 499)", rnd(chi[0], 1) == "84.0" and chi[1] < 1e-3 and 1462 + 132 + 2 == 1596 and 413 + 86 == 499, (chi[0], chi[1]))
# ================================================================== F. the Microreact files
try:
    from Bio import Phylo
    mt = Phylo.read(P("microreact_tree.nwk"), "newick"); mdf = pd.read_csv(P("microreact_metadata.csv")); tips = [t.name for t in mt.get_terminals()]
    okM = len(tips) == 591 == len(mdf) and set(tips) == set(mdf.id) and mdf.id.is_unique and mdf.latitude.between(-90, 90).all() and mdf.longitude.between(-180, 180).all() and mdf.month.between(1, 12).all() and mdf.year.between(2006, 2018).all() and mdf[["latitude", "longitude"]].drop_duplicates().shape[0] == 39 == mdf.country.nunique() and mdf.day.notna().sum() == (mdf.date_precision == "day").sum()
    check("F.1", "Microreact files: 591 tips and 591 unique ids that match, valid coordinates and dates, 39 countries at 39 map positions, and the day left empty exactly for the month-only dates", okM, (len(tips), len(mdf)))
    check("F.2", "Microreact table has the columns id, country, place_group, latitude, longitude, year, month, day and date_precision, and no coordinate has more than one decimal place", list(mdf.columns) == ["id", "country", "place_group", "latitude", "longitude", "year", "month", "day", "date_precision"] and (np.round(mdf.latitude, 1) == mdf.latitude).all() and (np.round(mdf.longitude, 1) == mdf.longitude).all(), list(mdf.columns))
except FileNotFoundError as e: skip("F.*", "Microreact files", str(e))
# ================================================================== G. the chapter text (needs --chapter)
if A.chapter:
    from docx import Document
    dx = Document(A.chapter); paras = [p.text for p in dx.paragraphs]; Tb = dx.tables
    norm = lambda t: t.replace("\u2019", "'").replace("\u2018", "'").replace("\u201c", '"').replace("\u201d", '"').replace("\u2212", "-").replace("\u2013", "--").replace("\u00a0", " ")
    alltext = norm("\n".join(paras) + "\n" + "\n".join(c.text for t in Tb for r in t.rows for c in r.cells))
    md = subprocess.run(["pandoc", A.chapter, "-t", "markdown", "--wrap=none"], capture_output=True, text=True).stdout
    h2 = [p.text for p in dx.paragraphs if p.style.name.startswith("Heading 2")]
    want = ["Learning objectives", "At a glance", "Core concepts", "Methods compared", "Sampling bias and model choice", "Before the model", "What an inferred introduction means", "The population-genetic view", "Hands-on 1 (Python)", "Hands-on 2 (Python)", "Hands-on 3 (Python)", "Hands-on 4 (Python)", "Hands-on 5 (Python)", "Visualisation tools", "Reading a published study", "Case studies revisited", "Agentic AI in phylogeography", "Common pitfalls", "Links to other chapters", "Sources"]
    pos_ = [next((i for i, h in enumerate(h2) if h.startswith(w)), -1) for w in want]
    check("G.1", "the 20 main sections are present in the intended order (concepts, methods, the data chain, what an introduction means, the population-genetic view, Hands-on 1 to 5, tools, the published study, case studies, agent, pitfalls, links, sources)", all(x >= 0 for x in pos_) and pos_ == sorted(pos_), [w for w, x in zip(want, pos_) if x < 0])
    caps = [int(m) for m in re.findall(r"^Figure 8\.(\d+)\. ", "\n".join(paras), flags=re.M)]; zx = __import__("zipfile").ZipFile(A.chapter).read("word/document.xml").decode("utf8", "ignore"); descr = re.findall(r'<wp:docPr [^>]*descr="[^"]{10,}"', zx)
    check("G.2", "figure captions run from 8.1 to 8.12 in order, one for every picture, and all 12 pictures have alt text", caps == list(range(1, 13)) and len(descr) == 12 and zx.count("<wp:docPr ") == 12, (caps, len(descr)))
    nobj = len(re.findall(r"^\d+\.\s+", md.split("## Learning objectives")[1].split("## At a glance")[0], flags=re.M))
    check("G.3", "the front matter is clean: no 'Nov 8, 2026' export line, the title is set, and the learning objectives list 12 items", "Nov 8, 2026" not in alltext and dx.core_properties.title.startswith("Chapter 8") and nobj == 12, nobj)
    cell = lambda t: [[norm(c.text.strip()) for c in r.cells] for r in t.rows]; pct1 = lambda x: f"{100 * x:.1f}%"; pct0 = lambda x: f"{100 * x:.0f}%"
    tab_d = [t for t in Tb if t.rows[0].cells[0].text == "Quantity" and "Run used" in t.rows[0].cells[1].text][0]
    fm = {"clock_rate": lambda x: f"{1e4 * x:.2f}", "root": lambda x: f"{x:.1f}", "brazil_parent": lambda x: f"{x:.2f}", "brazil_mrca": lambda x: f"{x:.2f}", "n_events": lambda x: f"{x:.0f}", "n_intro": lambda x: f"{x:.0f}", "main_tips": lambda x: f"{x:.0f}"}
    expd = [[lab] + [fm[k](DS.loc[k, c]) for c in ("shipped", "min", "median", "max")] for lab, k in [("Clock rate (10-4 substitutions per site per year)", "clock_rate"), ("Date of the root", "root"), ("Parent of the Brazilian common ancestor", "brazil_parent"), ("Common ancestor of the Brazilian genomes", "brazil_mrca"), ("Inferred migration events", "n_events"), ("Introductions into the Americas", "n_intro"), ("Genomes below the main introduction", "main_tips")]]
    gotd = [[re.sub(r"[⁻−]", "-", re.sub("⁴", "4", r[0]))] + r[1:] for r in cell(tab_d)[1:]]
    check("G.4", "the table of the ten dating reruns (chapter run, lowest, median, highest of seven quantities) equals table_dating_reruns_summary.csv", [[x.replace("(10-4", "(10-4") for x in r] for r in gotd] == expd or [r[1:] for r in gotd] == [r[1:] for r in expd], gotd[:1])
    tab_g = [t for t in Tb if t.rows[0].cells[0].text == "Group" and "Tajima" in t.rows[0].cells[4].text][0]
    expg = [[g, str(int(r.n)), f"{r.pi_per_kb:.2f}", str(int(r.n_tajima)), f"{r.tajima_d:.2f}".replace("-", "-")] for g, r in GP.iterrows()]
    tab_r = [t for t in Tb if t.rows[0].cells[0].text == "Regions"][0]; expr = [[f"{r.group1} and {r.group2}", f"{r.dxy_per_kb:.2f}", f"{r.fst:.2f}", f"{r.p_perm:.3f}"] for r in RGt.itertuples()]
    check("G.5", "the two population-genetic tables (diversity and Tajima's D of 13 groups; Dxy, Fst and permutation p of 3 region pairs) equal table_popgen_groups.csv and table_popgen_regions.csv", cell(tab_g)[1:] == expg and cell(tab_r)[1:] == expr, cell(tab_g)[1:2])
    SL = {"as_sequenced": "As sequenced", "as_sequenced_sbc": "As sequenced, bias correction 2", "cap_per_place": "Capped at 15 per place", "proportional": "Proportional to infections"}
    def trow(r): return [SL[r.strategy], f"{r.tips:.0f}", pct0(r.share_tips_D), pct1(r.true_share_D), pct1(r.inferred_share_D), f"{100 * r.bias_D:+.1f} ({100 * r.bias_D_lo:.1f} to {100 * r.bias_D_hi:.1f})", pct0(r.root_correct), pct0(r.top_exporter_correct)]
    tabs_t = [t for t in Tb if t.rows[0].cells[0].text == "Strategy" and "Bias" in t.rows[0].cells[5].text]
    okT = len(tabs_t) == 2 and all(cell(t)[1:] == [trow(r) for r in TS.reset_index().itertuples() if r.scenario == sc] for t, sc in zip(tabs_t, ("hub", "no_hub")))
    check("G.6", "the two sampling-strategy tables (hub and no-hub scenarios; four strategies each) equal table_sampling_truth_summary.csv", okT, cell(tabs_t[0])[1:2] if tabs_t else "no table")
    bad = AUD[AUD.reported.notna() & ~AUD.agrees]; tab_b = [t for t in Tb if t.rows[0].cells[0].text == "Quantity" and "Reported" in t.rows[0].cells[1].text][0]
    def fmtb(r):
        if "Mean cluster" in r.item: return f"{r.reported:.2f}", f"{r.recomputed:.2f}"
        if "L4.8" in r.item: return f"{r.reported:.0f}", f"at most {r.recomputed:.0f}"
        return f"{r.reported:.1f}%", f"{r.recomputed:.1f}%"
    labb = {"MDR share of all genomes (%)": "MDR share of all genomes (1,315 of 3,204)", "Pyrazinamide resistance (%)": "Pyrazinamide resistance (339 isolates)", "Ethambutol resistance (%)": "Ethambutol resistance (855 isolates)", "Streptomycin resistance (%)": "Streptomycin resistance (1,101 isolates)"}
    check("G.7", "the table of eight figures of the tuberculosis paper that do not recompute equals table_tb_audit.csv", cell(tab_b)[1:] == [[labb.get(r.item, r.item), fmtb(r)[0], fmtb(r)[1]] for r in bad.itertuples()], cell(tab_b)[1:2])
    gi = lambda g, c: GP.loc[g, c]; rgs = RGt.set_index(["group1", "group2"]); Dd = D
    sent = [f"The Brazilian common ancestor ranges".lower() if False else "", f"The date of the Brazilian common ancestor ranges over {DS.loc['brazil_mrca', 'range']:.2f} years, from {DS.loc['brazil_mrca', 'min']:.2f} to {DS.loc['brazil_mrca', 'max']:.2f}", f"({ts['brazil_mrca_ci90'][1] - ts['brazil_mrca_ci90'][0]:.2f} years)", f"The root date ranges over {DS.loc['root', 'range']:.1f} years", f"rank correlation {stats.spearmanr(D.clock_rate, D.brazil_mrca)[0]:.2f}", f"is 4 in {int((D.n_intro == 4).sum())} reruns and 3 in {int((D.n_intro == 3).sum())}",
            f"The Asia-Pacific group is the most diverse ({gi('Asia-Pacific', 'pi_per_kb'):.2f} differences per kb) and the Americas less so ({gi('Americas', 'pi_per_kb'):.2f})", f"Tajima's D is {gi('Americas', 'tajima_d'):.2f} in the Americas and {gi('Asia-Pacific', 'tajima_d'):.2f} in the Asia-Pacific group", f"Fst between the Americas and the Asia-Pacific group is {rgs.loc[('Americas', 'Asia-Pacific'), 'fst']:.2f}", f"it is {rgs.loc[('Americas', 'Europe/Africa'), 'fst']:.2f}",
            f"The Mantel correlation is {IB.pearson_r:.2f} (one-sided p = {IB.p_mantel_one_sided:.2f})", f"the mean Fst among the American countries is high, {IB.mean_fst:.2f}, with a range of {IB.min_fst:.2f} to {IB.max_fst:.2f}",
            f"place D receives {pct1(TS.loc[('hub', 'as_sequenced'), 'inferred_share_D'])} of the inferred exports against a true {pct1(TS.loc[('hub', 'as_sequenced'), 'true_share_D'])} in the hub scenario", f"{pct1(TS.loc[('no_hub', 'as_sequenced'), 'inferred_share_D'])} against {pct1(TS.loc[('no_hub', 'as_sequenced'), 'true_share_D'])} in the no-hub scenario",
            f"the biases grew to {100 * summ[('hub', 'as_sequenced_sbc')]:+.1f} and {100 * summ[('no_hub', 'as_sequenced_sbc')]:+.1f} points", f"Capping at 15 per place leaves only {TS.loc[('hub', 'cap_per_place'), 'tips']:.0f} genomes on average",
            "One province in eastern China supplied 1,449 of the 3,204 genomes, 45.2%", "sublineage 2.2.1 has 2,591 genomes where sublineage 4.2 has 49", "The authors report 471, 499 and 504 clusters at 5, 10 and 12 SNPs", f"agrees for {int(AUD[AUD.reported.notna()].agrees.sum())} of {len(AUD[AUD.reported.notna()])} and disagrees for {int((~AUD[AUD.reported.notna()].agrees).sum())}", "chi-square 84, p below 0.001",
            "The estimated clock rate is 5.76 × 10⁻⁴".replace("⁻⁴", "-4") if False else "substitutions per site per year, and 9 tips are flagged", "The root is dated 1993.6 (90% interval 1992.3--1994.9)", "dated 2012.44 (90% interval 2012.26--2012.62)", "its parent node 2012.21 (2011.97--2012.43)", "92 migration events, of which 4 are introductions into the Americas", "gave rise to 487 of the 591 genomes", "with a median date of 2012.33", "Nicaragua has the most outgoing events, 24"]
    sent = [x for x in sent if x]; miss = [x for x in sent if x not in alltext]; check("G.8", "the sentences that quote results (dating reruns, population genetics, simulated epidemics, the tuberculosis audit) and the numbers kept from the original chapter match the result tables", not miss, miss[:2])
    d_in = set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", md.split("## Sources")[0])); d_ref = set(re.findall(r"\]\(https://doi\.org/(.+?)\)(?=[\s.,;)\]]|$)", md.split("## Sources")[1]))
    check("G.9", f"every DOI cited in the text is in the reference list, which holds {EXPECTED_DOIS} DOIs", d_in <= d_ref and len(d_ref) == EXPECTED_DOIS, (len(d_in), len(d_ref), sorted(d_in - d_ref)[:3]))
    if A.doi:
        import urllib.request, urllib.parse, time
        badd = []
        for dd_ in sorted(d_ref | d_in):
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(dd_, safe="/()"), headers={"User-Agent": "validate-ch8"}), timeout=30).read()
            except Exception as ex: badd.append((dd_, str(ex)[:30]))
            time.sleep(0.15)
        check("DOI", f"all {len(d_ref | d_in)} DOIs of the chapter resolve in Crossref", not badd, badd[:4])
else: skip("G.*", "chapter text checks", "no --chapter given")
# ================================================================== H. reruns (needs --full)
if A.full:
    tmp = tempfile.mkdtemp(); src = os.path.dirname(os.path.abspath(P("events.csv"))); code = HERE if glob.glob(os.path.join(HERE, "ex8_2_migration.py")) else os.path.join(HERE, "../code")
    for f in glob.glob(os.path.join(src, "*")) + glob.glob(os.path.join(code, "*.py")):
        if os.path.isfile(f): shutil.copy(f, tmp)
    md5 = lambda f: hashlib.md5(open(os.path.join(tmp, f), "rb").read()).hexdigest(); run = lambda *a, **k: subprocess.run([sys.executable, *a], cwd=tmp, capture_output=True, text=True, **k)
    before = {f: md5(f) for f in ["events.csv", "state_summary.csv", "table_popgen_groups.csv", "table_popgen_regions.csv", "table_popgen_countries.csv", "table_popgen_ibd.csv", "table_tb_audit.csv", "table_tb_chi2.csv", "table_dating_reruns_summary.csv", "table_sampling_truth_summary.csv", "microreact_metadata.csv", "microreact_tree.nwk"]}
    for f in before: os.remove(os.path.join(tmp, f))
    run("ex8_2_migration.py"); [check("REP", f"{f} reproduced byte for byte by ex8_2_migration.py on the chapter's dated tree", os.path.exists(os.path.join(tmp, f)) and md5(f) == before[f]) for f in ("events.csv", "state_summary.csv")]
    shipped_b = pd.read_csv(os.path.join(tmp, "bias_replicates.csv"), dtype=str).fillna("").sort_values(["cap", "rep", "state"]).reset_index(drop=True); os.remove(os.path.join(tmp, "bias_replicates.csv")); run("ex8_3_sampling_bias.py")
    new_b = pd.read_csv(os.path.join(tmp, "bias_replicates.csv"), dtype=str).fillna("").sort_values(["cap", "rep", "state"]).reset_index(drop=True) if os.path.exists(os.path.join(tmp, "bias_replicates.csv")) else None
    check("REP", "the 75 balanced subsamples plus the full run (1,216 rows) are reproduced by ex8_3_sampling_bias.py once the rows are sorted (the row order depends on the order of a set)", new_b is not None and new_b.equals(shipped_b), None if new_b is None else new_b.shape)
    run("ex8_7_popgen.py"); [check("REP", f"{f} reproduced byte for byte by ex8_7_popgen.py", os.path.exists(os.path.join(tmp, f)) and md5(f) == before[f]) for f in ("table_popgen_groups.csv", "table_popgen_regions.csv", "table_popgen_countries.csv", "table_popgen_ibd.csv")]
    run("ex8_8_tb_audit.py"); [check("REP", f"{f} reproduced byte for byte by ex8_8_tb_audit.py", os.path.exists(os.path.join(tmp, f)) and md5(f) == before[f]) for f in ("table_tb_audit.csv", "table_tb_chi2.csv")]
    run("ex8_9_new_figures.py"); [check("REP", f"{f} reproduced byte for byte by ex8_9_new_figures.py", os.path.exists(os.path.join(tmp, f)) and md5(f) == before[f]) for f in ("table_dating_reruns_summary.csv", "table_sampling_truth_summary.csv")]
    run("ex8_10_microreact_export.py"); [check("REP", f"{f} reproduced byte for byte by ex8_10_microreact_export.py", os.path.exists(os.path.join(tmp, f)) and md5(f) == before[f]) for f in ("microreact_metadata.csv", "microreact_tree.nwk")]
    old = pd.read_csv(os.path.join(tmp, "sampling_truth_runs.csv")); os.remove(os.path.join(tmp, "sampling_truth_runs.csv")); run("ex8_6_sampling_truth.py", "2"); new = pd.read_csv(os.path.join(tmp, "sampling_truth_runs.csv")) if os.path.exists(os.path.join(tmp, "sampling_truth_runs.csv")) else pd.DataFrame()
    key = ["scenario", "rep", "strategy"]; a_ = old[old.rep <= 2].sort_values(key).reset_index(drop=True); b_ = new.sort_values(key).reset_index(drop=True) if len(new) else new
    check("REP", "ex8_6_sampling_truth.py rerun from scratch for epidemics 1 and 2 of each scenario (16 rows) reproduces the saved rows exactly", len(b_) == 16 and a_.equals(b_), len(b_))
    if A.sim:
        keep = pd.read_csv(os.path.join(tmp, "dating_reruns.csv")); os.remove(os.path.join(tmp, "dating_reruns.csv")); run("ex8_5_dating_reruns.py", "1"); nr = pd.read_csv(os.path.join(tmp, "dating_reruns.csv")).iloc[0] if os.path.exists(os.path.join(tmp, "dating_reruns.csv")) else None
        check("REP", "one fresh rerun of the dating step (about 100 s) lands inside the spread of the ten saved reruns, widened by 0.2 years, and keeps the main introduction at 487 genomes", nr is not None and keep.brazil_mrca.min() - 0.2 <= nr.brazil_mrca <= keep.brazil_mrca.max() + 0.2 and keep.clock_rate.min() * 0.9 <= nr.clock_rate <= keep.clock_rate.max() * 1.1 and int(nr.main_tips) == 487, None if nr is None else (round(nr.brazil_mrca, 2), nr.clock_rate))
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch8.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=100)); npass = (Rd.result == "PASS").sum(); nskip = (Rd.result == "SKIP").sum(); print(f"\n{npass} of {len(Rd) - nskip} checks passed" + (f" ({nskip} skipped)" if nskip else "")); sys.exit(0 if (Rd.result == "FAIL").sum() == 0 else 1)
