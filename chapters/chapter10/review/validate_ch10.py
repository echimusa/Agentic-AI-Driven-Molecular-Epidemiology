"""
Chapter 10 validation harness. Recomputes the numbers quoted in the chapter with code that does not import the chapter's own scripts (ex10_*.py).
A  the original chapter: loci from the summary statistics, odds ratios and P values, the inflation rescaling, the catalogue tables, the 3p21.31 credible set, the stratification table
B  the new exercises: QC steps and the injected-fault test, Hardy-Weinberg exact tests from genotype counts, Hudson FST, the pattern of LD, power, the inflation scaling, the 3p21.31 table,
   admixture mapping (formula and a fresh simulation) and the audit of the course notes
G  the chapter text (tables, figures, objectives, sentences, references)
H  figures and scripts
I  reruns of the scripts that run offline
USAGE   python validate_ch10.py [--chapter Chapter_10_....docx] [--doi] [--full]
        --chapter  also checks the chapter's tables, figures, objectives, sentences and references against the results (needs python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns the offline scripts in a scratch folder and checks that their result tables are reproduced (byte-identical, except the credible-set table, compared numerically); about 1 minute
Prints PASS/FAIL/SKIP per check, writes validation_results_ch10.csv and exits 1 if any check fails.
"""
import argparse, glob, gzip, hashlib, io, math, os, re, shutil, subprocess, sys, tempfile, warnings
import numpy as np, pandas as pd
from scipy import stats, optimize
from scipy.special import gammaln
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
def P(name):
    for d in (".", "../data", "data", os.path.join(HERE, "../data"), os.path.join(HERE, "data")):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []; EXPECTED_DOIS = 0
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
GROUPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
# ================================================================== I. reruns first, while this process is small
tmp = None
if A.full:
    tmp = tempfile.mkdtemp(prefix="ch10_"); data = os.path.dirname(os.path.abspath(P("chr22_group_geno_counts.tsv.gz"))); src = os.path.join(HERE, "../code") if os.path.isdir(os.path.join(HERE, "../code")) else HERE
    for f_ in glob.glob(os.path.join(src, "*.py")) + glob.glob(os.path.join(src, "*.R")): shutil.copy(f_, tmp)
    for item in os.listdir(data):
        if item.startswith("table_"): continue
        s_, d_ = os.path.join(data, item), os.path.join(tmp, item); shutil.copytree(s_, d_) if os.path.isdir(s_) else shutil.copy(s_, d_)
    os.makedirs(os.path.join(tmp, "figs10"), exist_ok=True)
    shutil.copy(os.path.join(data, "table_credible_set.csv"), os.path.join(tmp, "table_credible_set_shipped.csv")) if os.path.exists(os.path.join(data, "table_credible_set.csv")) else None
    for name, cmd in [("ex10_2b", ["Rscript", "ex10_2b_r_check.R"]), ("ex10_4", ["python3", "ex10_4_finemap.py"]), ("ex10_6", ["python3", "ex10_6_power_lambda.py"]), ("ex10_7", ["python3", "ex10_7_ancestry_3p21.py"]), ("ex10_8", ["python3", "ex10_8_admixture_mapping.py"]), ("ex10_9", ["python3", "ex10_9_notes_audit.py"]), ("ex10_10", ["python3", "ex10_10_diagrams.py"])]:
        r_ = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True); check(f"I.{name}", f"{name} reruns from the shipped data without error", r_.returncode == 0, (r_.stderr or "")[-120:])
    diffs, near = [], True
    for tb in sorted(os.path.basename(f_) for f_ in glob.glob(os.path.join(data, "table_*.csv"))):
        a_, b_ = os.path.join(tmp, tb), os.path.join(data, tb)
        if not os.path.exists(a_): continue                                        # tables of the scripts that need downloads are not rerun
        if tb == "table_credible_set.csv":
            x, y = pd.read_csv(a_), pd.read_csv(b_); near = x.shape == y.shape and (x.rsid == y.rsid).all() and np.allclose(x.select_dtypes("number").values, y.select_dtypes("number").values, rtol=1e-9)
        elif hashlib.md5(open(a_, "rb").read()).hexdigest() != hashlib.md5(open(b_, "rb").read()).hexdigest(): diffs.append(tb)
    check("I.tables", "every result table of the rerun is byte-identical to the shipped one (the credible-set table equal to 9 digits)", not diffs and near, (diffs, near))
    shutil.rmtree(tmp, ignore_errors=True)
else: skip("I.*", "reruns of the offline scripts", "run with --full")
# ================================================================== A. the original chapter
ss = pd.read_csv(P("sumstats_p1e-4.tsv.gz"), sep="\t"); loci_t = pd.read_csv(P("table_loci.csv"))
sig = ss[ss.p < 5e-8].sort_values("p").reset_index(drop=True); taken = np.zeros(len(sig), bool); leads = []
for i in range(len(sig)):
    if taken[i]: continue
    near_ = (sig.chr.values == sig.chr[i]) & (np.abs(sig.pos.values - sig.pos[i]) <= 500_000); taken |= near_; leads.append(i)
L = sig.loc[leads]
check("A.1", "greedy 500 kb clumping of the genome-wide significant variants gives 7 loci with the lead variants of the chapter (rs35081325, rs13050728, rs2109069, rs1894387, rs2660, rs622568, rs920065566)",
      len(L) == 7 and list(L.rsid) == ["rs35081325", "rs13050728", "rs2109069", "rs1894387", "rs2660", "rs622568", "rs920065566"] and list(loci_t.rsid) == list(L.rsid), list(L.rsid))
wald = 2 * stats.norm.sf(np.abs(L.beta / L.se)); dl = np.abs(np.log10(wald) - np.log10(L.p.values)).max()
check("A.2", "the odds ratios of the leads are exp(beta) (rs35081325 1.629, 95% CI 1.532 to 1.733) and the reported P agrees with the Wald P within 0.002 log10 units", [rnd(math.exp(b), 3) for b in L.beta] == [rnd(o, 3) for o in loci_t.OR]
      and rnd(math.exp(L.beta.iloc[0]), 3) == "1.629" and rnd(math.exp(L.beta.iloc[0] - 1.96 * L.se.iloc[0]), 3) == "1.532" and rnd(math.exp(L.beta.iloc[0] + 1.96 * L.se.iloc[0]), 3) == "1.733" and dl < 0.002, dl)
ncase, nctl = 9986, 1877672; lam1000 = 1 + (1.045 - 1) * (1 / ncase + 1 / nctl) / (1 / 1000 + 1 / 1000)
check("A.3", "the inflation factor 1.045 rescaled to 1,000 cases and 1,000 controls is 1.002", rnd(lam1000, 3) == "1.002", lam1000)
cg = pd.read_csv(P("table_catalog_groups.csv")); ca = pd.read_csv(P("table_catalog_ancestry.csv"))
tot = cg[cg.group == "All infection"].iloc[0]; parts = cg[~cg.group.isin(["All infection", "All traits"])]
check("A.4", "the catalogue groups add up to 1,351 infection studies, 7,646 associations and 559 studies with full summary statistics (41.4%), and the ancestry shares add up to about 100%",
      int(parts.studies.sum()) == 1351 == int(tot.studies) and int(parts.associations.sum()) == 7646 and int(parts.full_summary_stats.sum()) == 559 == int(tot.full_summary_stats) and rnd(100 * 559 / 1351, 1) == "41.4"
      and np.allclose(ca.iloc[:, 1:7].sum(axis=1), 100, atol=0.3), (parts.studies.sum(), ca.iloc[:, 1:7].sum(axis=1).round(1).tolist()))
# fine-mapping with own code
CHR, LEAD_POS = 3, 45_848_429; LO, HI = LEAD_POS - 500_000, LEAD_POS + 500_000
r_ = pd.read_csv(P("hgi_chr3_region.tsv.gz"), sep="\t", keep_default_na=False, na_values=[""]); r_ = r_[(r_.chrom == CHR) & (r_.pos >= LO) & (r_.pos <= HI) & np.isfinite(r_.beta) & (r_.se > 0) & (r_.a1.str.len() == 1) & (r_.a2.str.len() == 1)].drop_duplicates("pos", keep=False).reset_index(drop=True)
G = pd.read_csv(P("chr3_region_eur.gt.gz"), sep="\t", header=None, dtype=str); pos_ = G[0].astype(int).values; ref_, alt_ = G[1].values, G[2].values; dos = G.iloc[:, 3:].apply(lambda c: c.str.count("1")).values.astype(np.float32)
m = pd.DataFrame({"pos": pos_, "ref": ref_, "alt": alt_, "row": range(len(pos_))}); J = r_.merge(m, on="pos"); J = J[((J.a1 == J.alt) & (J.a2 == J.ref)) | ((J.a1 == J.ref) & (J.a2 == J.alt))].reset_index(drop=True)
X = dos[J.row.values]; X = np.where((J.a1 == J.alt).values[:, None], X, 2 - X); fq = X.mean(axis=1) / 2; ok_ = (fq > 0.005) & (fq < 0.995); J = J[ok_].reset_index(drop=True); X = X[ok_]
lead = J.index[J.rsid == "rs35081325"][0]; r2 = np.array([np.corrcoef(X[lead], X[i])[0, 1] ** 2 for i in range(len(J))])
def cs(sd):
    V, W = J.se.values ** 2, sd ** 2; rr = W / (V + W); lb = 0.5 * np.log(1 - rr) + 0.5 * rr * (J.beta.values / J.se.values) ** 2; pp = np.exp(lb - lb.max()); pp /= pp.sum(); o = np.argsort(-pp); k = int(np.searchsorted(np.cumsum(pp[o]), 0.95) + 1); return pp, o[:k]
out = {sd: cs(sd) for sd in (0.1, 0.15, 0.2)}; pp, sel = out[0.15]; cs_t = pd.read_csv(P("table_credible_set.csv"))
check("A.5", "an independent fine-mapping gives a 95% credible set of 11 variants at all three priors, lead rs35081325 with posterior probability 0.167 (0.161 to 0.170), r2 above 0.98 with the lead, span 45,818,159 to 45,858,159",
      [len(out[sd][1]) for sd in out] == [11, 11, 11] and J.rsid[sel[0]] == "rs35081325" and [rnd(out[sd][0][lead], 3) for sd in out] == ["0.161", "0.167", "0.170"] and (r2[sel] > 0.98).all() and J.pos[sel].min() == 45818159 and J.pos[sel].max() == 45858159
      and list(J.rsid[sel[:3]]) == ["rs35081325", "rs34326463", "rs34288077"] and set(cs_t.rsid) == set(J.rsid[sel]), ([len(out[sd][1]) for sd in out], [rnd(out[sd][0][lead], 3) for sd in out]))
st_ = pd.read_csv(P("table_stratification.csv")) if os.path.exists(os.path.join(HERE, "../data/table_stratification.csv")) or os.path.exists("table_stratification.csv") or os.path.exists("data/table_stratification.csv") else None
if st_ is not None:
    r0 = st_.iloc[0]; check("A.6", "the stratification table: lambda 71.68, 27,658 of 54,197 variants away from the causal SNPs significant (51.0%) with no covariates, none with 2 or 10 components", rnd(r0.lambda_all, 2) == "71.68" and int(r0.false_positives) == 27658 and int(r0.tested_away) == 54197 and rnd(100 * 27658 / 54197, 1) == "51.0" and (st_.false_positives.iloc[1:] == 0).all(), r0.to_dict())
else: skip("A.6", "stratification table", "table_stratification.csv not found")
# ================================================================== B. the new exercises
qs = pd.read_csv(P("table_qc_steps.csv")); qi = pd.read_csv(P("table_qc_injected.csv")).iloc[0]
check("B.1", "QC steps: --geno and --mind remove nothing, --maf 0.01 removes 83,688 of 100,000 variants, pooled --hwe 1e-6 removes 4,056, all four leave 12,306; the injected faults (10 people, 40 variants) are found exactly, leaving 1,960 variants and 490 people",
      (qs.variants_left + qs.variants_removed == 100000).all() and list(qs.variants_removed) == [0, 0, 83688, 4056, 87694] and qs.people_removed.sum() == 0 and qs.variants_left.iloc[-1] == 12306 and (qi.injected_people, qi.found_people, qi.injected_variants, qi.found_variants, qi.variants_left, qi.people_left) == (10, 10, 40, 40, 1960, 490), qs.variants_removed.tolist())
cnt = pd.read_csv(P("chr22_group_geno_counts.tsv.gz"), sep="\t"); allg = {k: sum(cnt[f"{g}_{k}"] for g in GROUPS) for k in ("n11", "n12", "n22")}
def hwe_p(n11, n12, n22):
    N = n11 + n12 + n22; a = 2 * n11 + n12; b = 2 * n22 + n12
    if N == 0 or a == 0 or b == 0: return 1.0
    hmax = min(a, b); hs = np.arange(a % 2, hmax + 1, 2); lp = gammaln(N + 1) + hs * math.log(2) + gammaln(a + 1) + gammaln(b + 1) - gammaln((a - hs) / 2 + 1) - gammaln(hs + 1) - gammaln((b - hs) / 2 + 1) - gammaln(2 * N + 1)
    lo = lp[(hs == n12).argmax()]; pr = np.exp(lp - lp.max()); pr /= pr.sum(); return float(min(1.0, pr[lp <= lo * (1 - 1e-12) + 1e-12 * abs(lo) + 1e-9].sum()))
def fails(n11, n12, n22, mask, thr=1e-6):
    N = (n11 + n12 + n22).values.astype(float); p = (2 * n11 + n12).values / (2 * np.maximum(N, 1)); e = np.stack([N * p * p, 2 * N * p * (1 - p), N * (1 - p) ** 2]); o = np.stack([n11.values, n12.values, n22.values]).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"): chi = np.nansum(np.where(e > 0, (o - e) ** 2 / e, 0), axis=0)
    pc = stats.chi2.sf(chi, 1); res = np.zeros(len(n11), bool); band = mask & (pc < 1e-2)      # the exact test is run for every variant the chi-square test does not clearly pass: for rare alleles the chi-square P can be far smaller than the exact P
    for i in np.flatnonzero(band): res[i] = hwe_p(int(n11.values[i]), int(n12.values[i]), int(n22.values[i])) < thr
    return res
maf_pool = (2 * allg["n22"] + allg["n12"]) / (2 * (allg["n11"] + allg["n12"] + allg["n22"])); mask = (np.minimum(maf_pool, 1 - maf_pool) >= 0.01).values
fp = fails(allg["n11"], allg["n12"], allg["n22"], mask); fg = {g: fails(cnt[f"{g}_n11"], cnt[f"{g}_n12"], cnt[f"{g}_n22"], mask) for g in GROUPS}; anyg = np.logical_or.reduce(list(fg.values()))
ht = pd.read_csv(P("table_hwe_pooled_within.csv")).set_index("test")
check("B.2", "Hardy-Weinberg exact tests recomputed from the genotype counts: 16,312 variants with frequency of at least 1%; failing at P < 1e-6: pooled 4,006; AFR 227, AMR 140, EAS 161, EUR 140, SAS 163; 3,726 fail pooled but in no group; 28 fail in a group only",
      int(mask.sum()) == 16312 and int(fp.sum()) == 4006 == int(ht.loc["pooled sample", "fail_p_1e6"]) and [int(fg[g].sum()) for g in GROUPS] == [227, 140, 161, 140, 163] and int((fp & ~anyg).sum()) == 3726 and int((~fp & anyg & mask).sum()) == 28, (int(fp.sum()), [int(fg[g].sum()) for g in GROUPS], int((fp & ~anyg).sum())))
com = set(open(P("common.txt")).read().split()); isc = cnt.SNP.isin(com).values; ft = pd.read_csv(P("table_fst_pairs.csv"))
def pf(g): n = cnt[f"{g}_n11"] + cnt[f"{g}_n12"] + cnt[f"{g}_n22"]; return ((2 * cnt[f"{g}_n11"] + cnt[f"{g}_n12"]) / (2 * n)).values, (2 * n).values.astype(float)
def hud(a, b, mk):
    p1, n1 = pf(a); p2, n2 = pf(b); num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1); den = p1 * (1 - p2) + p2 * (1 - p1); g = (den > 0) & mk; return float(num[g].sum() / den[g].sum())
ok_f = all(abs(hud(r.group_a, r.group_b, np.ones(len(cnt), bool)) - r.hudson_fst_all) < 1e-6 and abs(hud(r.group_a, r.group_b, isc) - r.hudson_fst_common) < 1e-6 for r in ft.itertuples()); ee = ft[(ft.group_a == "EAS") & (ft.group_b == "EUR")].iloc[0]
check("B.3", "Hudson FST between all ten pairs of groups recomputed from the counts (all variants and common variants); East Asia against Europe 0.110 on common variants, equal to PLINK's weighted Weir-Cockerham 0.1097; Africa against East Asia 0.180, Europe against the Americas 0.027",
      ok_f and rnd(ee.hudson_fst_common, 3) == "0.110" and abs(ee.hudson_fst_common - ee.plink_wc_weighted_common) < 5e-4 and rnd(ft[(ft.group_a == "AFR") & (ft.group_b == "EAS")].hudson_fst_common.iloc[0], 3) == "0.180" and rnd(ft[(ft.group_a == "AMR") & (ft.group_b == "EUR")].hudson_fst_common.iloc[0], 3) == "0.027", ee.to_dict())
ld = pd.read_csv(P("table_ld_by_distance.csv")); pv = ld.pivot(index="bin_kb", columns="group", values="mean_r2")
check("B.4", "mean r2 falls with distance in every group and is lowest in AFR in every bin (0.182 against 0.262 to 0.293 at 0 to 10 kb)", all(pv.loc[["0-10", "10-25", "25-50", "50-100", "100-250"], g].is_monotonic_decreasing for g in GROUPS) and all(pv.loc[b].idxmin() == "AFR" for b in pv.index) and rnd(pv.loc["0-10", "AFR"], 3) == "0.182" and rnd(pv.loc["0-10"].drop("AFR").min(), 3) == "0.262" and rnd(pv.loc["0-10"].drop("AFR").max(), 3) == "0.293", pv.round(3).to_dict())
# power with own root finding
CRIT = stats.chi2.isf(5e-8, 1)
def pw(nca, nco, p0, orr, crit=CRIT):
    p1 = orr * p0 / (1 - p0 + orr * p0); pb = (nca * p1 + nco * p0) / (nca + nco); return stats.ncx2.sf(crit, 1, (p1 - p0) ** 2 / (pb * (1 - pb) * (1 / (2 * nca) + 1 / (2 * nco))))
need = pd.read_csv(P("table_power_cases_needed.csv")); mino = pd.read_csv(P("table_power_min_or.csv")); ok_n = True
for r in need.itertuples(index=False):
    for j, mf in enumerate((0.01, 0.05, 0.20, 0.50)):
        x = optimize.brentq(lambda n: pw(n, n, mf, r[0]) - 0.8, 50, 1e8, xtol=0.01)
        if abs(x - r[j + 1]) > 1.0: ok_n = False
ok_o = True
for r in mino.itertuples(index=False):
    for (a, b), val in (((1000, 1000), r[1]), ((9986, 1877672), r[2]), ((10000, 10000), r[3])):
        x = optimize.brentq(lambda o: pw(a, b, r[0], o) - 0.8, 1.0001, 30, xtol=1e-6)
        if abs(x - val) > 1e-3: ok_o = False
check("B.5", "power recomputed with an independent root finder: cases needed for 80% power at P < 5e-8 (OR 1.2 at frequency 0.01 needs 110,213; OR 1.1 at 0.5 needs 17,454) and the smallest detectable odds ratio (2.18 at 5% in 1,000 + 1,000; 1.21 in the HGI file)",
      ok_n and ok_o and int(need[need.odds_ratio == 1.2].iloc[0, 1]) == 110213 and int(need[need.odds_ratio == 1.1].iloc[0, 4]) == 17454 and rnd(mino[mino.maf == 0.05].min_or_1000_1000.iloc[0], 2) == "2.18" and rnd(mino[mino.maf == 0.05].min_or_hgi_r5.iloc[0], 2) == "1.21", (ok_n, ok_o))
rng = np.random.default_rng(2026); mc_ok = True; det = []
for nca, nco, p0, orr in [(5000, 5000, 0.2, 1.3), (20000, 20000, 0.3, 1.1)]:
    p1 = orr * p0 / (1 - p0 + orr * p0); ca_ = rng.binomial(2 * nca, p1, 3000); co_ = rng.binomial(2 * nco, p0, 3000); pb = (ca_ + co_) / (2 * nca + 2 * nco); chi = (ca_ / (2 * nca) - co_ / (2 * nco)) ** 2 / (pb * (1 - pb) * (1 / (2 * nca) + 1 / (2 * nco))); emp = float((chi > CRIT).mean()); det.append((round(emp, 3), round(pw(nca, nco, p0, orr), 3)))
    if abs(emp - pw(nca, nco, p0, orr)) > 0.04: mc_ok = False
mct = pd.read_csv(P("table_power_mc.csv")); check("B.6", "a fresh Monte Carlo (3,000 repeats, own seed) agrees with the formula within 0.04, and the shipped Monte Carlo table agrees with the formula within 0.03", mc_ok and (abs(mct.power_formula - mct.power_monte_carlo) < 0.03).all(), det)
ex = pd.read_csv(P("table_example_allelic.csv")).iloc[0]; chi_ex = 2000 * (400 * 700 - 600 * 300) ** 2 / (1000 * 1000 * 700 * 1300)
check("B.7", "the allelic example of the notes (400 and 600 case alleles, 300 and 700 control alleles): chi-square 21.98, odds ratio 1.56, P 2.8e-6, power 0.22 with 500 + 500 at P < 5e-8", rnd(chi_ex, 2) == "21.98" and rnd(stats.chi2.sf(chi_ex, 1) * 1e6, 1) == "2.8" and rnd(400 * 700 / (600 * 300), 2) == "1.56" and rnd(pw(500, 500, 0.3, 400 * 700 / (600 * 300)), 2) == "0.22" and abs(ex.chi_square - chi_ex) < 1e-6, chi_ex)
lt = pd.read_csv(P("table_lambda_scaling.csv")); pe, pa = pf("EUR")[0], pf("EAS")[0]; ok_l = True; Fs = []
for qca, qco in [(0.30, 0.10), (0.11, 0.09)]:
    pca = qca * pa + (1 - qca) * pe; pco = qco * pa + (1 - qco) * pe; pb = (pca + pco) / 2; k = (pb > 0) & (pb < 1) & isc; Fx = float(np.mean((pca[k] - pco[k]) ** 2 / (2 * pb[k] * (1 - pb[k])))); Fs.append(Fx)
    row = lt[(lt.cases_eas_share == qca) & (lt.people == 10000)].iloc[0]
    if abs(row.F - Fx) > 1e-6 * max(1, row.F) * 10 and abs(row.F - Fx) / row.F > 0.01: ok_l = False
k = None
pca = 0.11 * pa + 0.89 * pe; pco = 0.09 * pa + 0.91 * pe; kk = (((pca + pco) / 2) > 0) & (((pca + pco) / 2) < 1); sim = []
for N in (1000, 10000):
    ca_ = rng.binomial(N, pca[kk]); co_ = rng.binomial(N, pco[kk]); q = (ca_ + co_) / (2 * N); okq = (q > 0) & (q < 1); chi = (ca_[okq] / N - co_[okq] / N) ** 2 / (q[okq] * (1 - q[okq]) * (2 / N)); pb_ = (pca[kk] + pco[kk]) / 2; Fn = float(np.mean((pca[kk] - pco[kk]) ** 2 / (2 * pb_ * (1 - pb_)))); sim.append((float(chi.mean()), 1 + N * Fn))
ne = 4 / (1 / 9986 + 1 / 1877672); hg = lt.dropna(subset=["lambda_hgi_rescaled"]).iloc[0]
check("B.8", "inflation scaling: F recomputed from the counts agrees with the table; a fresh simulation gives a mean chi-square within 8% of 1 + N F at N = 1,000 and 10,000; the table shows lambda 1.05, 1.48 and 5.2 for the milder imbalance; HGI: effective N 39,733, implied F 1.1e-6, lambda 1.002 at 1,000 + 1,000",
      ok_l and all(abs(a - b) / b < 0.08 for a, b in sim) and [rnd(x, 2) for x in lt[(lt.cases_eas_share == 0.11)].lambda_median] == ["1.05", "1.48", "5.17"] and int(round(ne)) == 39733 and rnd(hg.F * 1e6, 2) == "1.13" and rnd(hg.lambda_hgi_rescaled, 3) == "1.002", (ok_l, sim))
# 3p21.31
pan = pd.read_csv(P("integrated_call_samples_v3.20130502.ALL.panel"), sep="\t").set_index("sample"); smp = open(P("chr3_region_samples.txt")).read().split(); popn = pan.super_pop.reindex(smp).values
g3 = pd.read_csv(P("chr3_region_all.gt.gz"), sep="\t", header=None, dtype=str); info = g3.iloc[:, :3]; info.columns = ["pos", "ref", "alt"]; info["pos"] = info.pos.astype(int); alt = g3.iloc[:, 3:].apply(lambda c: c.str.count("1")).values.astype(np.int8)
lead_row = int(np.flatnonzero((info.pos == 45848429).values)[0]); ref_l, alt_l = info.ref[lead_row], info.alt[lead_row]; eff_is_alt = (cs_t[cs_t.rsid == "rs35081325"].a1.iloc[0] == alt_l); ld_ = alt[lead_row] if eff_is_alt else 2 - alt[lead_row]
t3 = pd.read_csv(P("table_3p21_by_group.csv")); okt = True; r3 = []
for g in GROUPS:
    mm = popn == g; f = ld_[mm].mean() / 2; c_ = (ld_[mm] >= 1).mean(); row = t3[t3.group == g].iloc[0]; r3.append((g, round(f, 3), round(c_, 3)))
    sd = alt[:, mm].astype(float); sd_l = ld_[mm].astype(float); strong = 0
    for i in range(len(info)):
        if i == lead_row or sd[i].std() == 0 or sd_l.std() == 0: continue
        if np.corrcoef(sd_l, sd[i])[0, 1] ** 2 >= 0.8: strong += 1
    if abs(row.risk_allele_freq - f) > 1e-6 or abs(row.carriers - c_) > 1e-6 or int(row.others_r2_ge_0_8) != strong: okt = False
check("B.9", "3p21.31 by group recomputed from the phased genotypes: carriers of the risk allele 0.8% AFR, 8.4% AMR, 1.0% EAS, 15.9% EUR, 49.7% SAS; frequency 0.004, 0.043, 0.005, 0.080, 0.297; variants in the window with r2 of at least 0.8 with the lead: 18, 23, 14, 24, 23",
      okt and [rnd(100 * x[2], 1) for x in r3] == ["0.8", "8.4", "1.0", "15.9", "49.7"] and list(t3.others_r2_ge_0_8) == [18, 23, 14, 24, 23], r3)
rep_t = pd.read_csv(P("table_3p21_replication.csv")); ok_r = True
for r in rep_t.itertuples():
    if r.risk_allele_freq <= 0: continue
    for a_, col in ((0.05, "cases_p_0_05"), (5e-8, "cases_p_5e8")):
        x = optimize.brentq(lambda n: pw(n, n, r.risk_allele_freq, 1.63, stats.chi2.isf(a_, 1)) - 0.8, 20, 1e9, xtol=0.01)
        if abs(x - getattr(r, col)) > 1.0: ok_r = False
check("B.10", "cases needed to replicate OR 1.63 with 80% power, recomputed: at P < 0.05, SAS 146, EUR 371, AMR 644, EAS 5,283, AFR 6,916", ok_r and [int(round(x)) for x in rep_t.cases_p_0_05] == [6916, 644, 5283, 371, 146], rep_t.cases_p_0_05.round(0).tolist())
rc = pd.read_csv(P("table_3p21_r2_credible.csv")); lowr = rc[(rc.group == "AFR") & (rc.r2_with_lead < 0.95)]
check("B.11", "in the credible set every variant has r2 of at least 0.98 (0.985 at the lowest) with the lead in Europeans, and only rs35731912 (0.09) and rs67959919 (0.83) fall below 0.95 in Africans", (rc[rc.group == "EUR"].r2_with_lead >= 0.98).all() and rnd(rc[rc.group == "EUR"].r2_with_lead.min(), 3) == "0.985" and sorted(lowr.rsid) == ["rs35731912", "rs67959919"] and rnd(lowr.r2_with_lead.min(), 2) == "0.09", lowr.to_dict("records"))
# admixture mapping
ao = pd.read_csv(P("table_admix_omega.csv")); okf = all(abs((1.5 * r.pE + 1 - r.pE) / (1.5 * r.pA + 1 - r.pA) - r.omega_formula) < 1e-9 and abs(r.omega_formula - r.omega_simulated) < 0.03 for r in ao.itertuples())
rng2 = np.random.default_rng(77); n = 1_500_000; sims = []
for pE, pA in ((0.9, 0.1), (0.6, 0.4)):
    l1 = rng2.random(n) < 0.5; l2 = rng2.random(n) < 0.5; k = (rng2.random(n) < np.where(l1, pE, pA)).astype(int) + (rng2.random(n) < np.where(l2, pE, pA)).astype(int); w = 1.5 ** k; case = rng2.random(n) < 0.02 * w / w.mean()
    ec = l1[case].sum() + l2[case].sum(); ac = 2 * case.sum() - ec; eo = l1[~case].sum() + l2[~case].sum(); ao_ = 2 * (~case).sum() - eo; sims.append(((ec / ac) / (eo / ao_), (1.5 * pE + 1 - pE) / (1.5 * pA + 1 - pA)))
check("B.12", "ancestry odds ratios: the six table values (1.50, 0.67, 1.38, 1.22, 1.08, 1.00) equal the formula and the shipped simulation within 0.03; a fresh simulation of two rows agrees with the formula within 0.03", okf and all(abs(a - b) < 0.03 for a, b in sims) and [rnd(x, 2) for x in ao.omega_formula] == ["1.50", "0.67", "1.38", "1.22", "1.08", "1.00"], sims)
apw = pd.read_csv(P("table_admix_power.csv"))
def pc2(x, y, z):
    zc = np.column_stack([np.ones(len(z)), z]); rx = x - zc @ np.linalg.lstsq(zc, x, rcond=None)[0]; ry = y - zc @ np.linalg.lstsq(zc, y, rcond=None)[0]; return (len(y) - 2) * (rx @ ry) ** 2 / ((rx @ rx) * (ry @ ry))
kk_, hits = 1.5 ** 0, 0; crit_am = stats.chi2.isf(1e-5, 1); a0, b0 = 1.566, 5.89; reps = 150
for _ in range(reps):
    th = rng2.beta(a0, b0, 60000); l1 = rng2.random(60000) < th; l2 = rng2.random(60000) < th; k = (rng2.random(60000) < np.where(l1, 0.7, 0.1)).astype(float) + (rng2.random(60000) < np.where(l2, 0.7, 0.1)).astype(float)
    case = rng2.random(60000) < 1 / (1 + np.exp(-(-3.0 + math.log(1.5) * k))); ic = rng2.choice(np.flatnonzero(case), 2000, replace=False); io = rng2.choice(np.flatnonzero(~case), 2000, replace=False); idx = np.concatenate([ic, io]); y = np.r_[np.ones(2000), np.zeros(2000)]
    hits += pc2((l1[idx].astype(float) + l2[idx]), y, th[idx]) > crit_am
emp = hits / reps; p06 = apw[apw.difference == 0.6].power_admixture_mapping.iloc[0]
check("B.13", "admixture mapping power: 0 at a frequency difference of 0.2 or less, 7% at 0.4, 53% at 0.6, 89% at 0.8; a fresh simulation at 0.6 (150 repeats) is within 0.12 of the shipped 0.53; the causal-variant test has 73% to 99% power", (apw[apw.difference <= 0.2].power_admixture_mapping == 0).all() and [rnd(100 * x, 0) for x in apw.power_admixture_mapping[3:]] == ["7", "53", "89"] and abs(emp - p06) < 0.12 and rnd(100 * apw.power_causal_snp.min(), 0) == "73" and rnd(100 * apw.power_causal_snp.max(), 0) == "99", (emp, p06))
au = pd.read_csv(P("table_notes_audit.csv")); su = pd.read_csv(P("table_notes_audit_summary.csv")); chi_ = sum((o - e) ** 2 / e for o, e in ((298, 294.3), (489, 496.4), (213, 209.3)))
check("B.14", "audit of the course sources: 37 figures recomputed, 35 agree; the two that do not are cystic fibrosis (printed 500, 0.25 / 0.0004 = 625) and the Wellcome Trust study total (printed 23, parts add to 24); the Hardy-Weinberg example, the r2 and D', and all 12 admixture figures agree",
      len(au) == 37 and int(au.agrees.sum()) == 35 and sorted(au[~au.agrees].printed) == [23.0, 500.0] and abs(0.25 / 0.0004 - 625) < 1e-9 and 0 + 1 + 9 + 7 + 3 + 4 == 24 and rnd(chi_, 2) == "0.22" and int(su.figures.sum()) == 37 and int(su[su.source.str.contains("admixture")].agree.iloc[0]) == 12, au[~au.agrees].printed.tolist())
# ================================================================== G. the chapter text
def docx_text(path):
    from docx import Document
    d = Document(path); body = [p.text for p in d.paragraphs]; tabs = [[[c.text for c in r.cells] for r in t_.rows] for t_ in d.tables]; return d, "\n".join(body), tabs
if A.chapter:
    d, txt, tabs = docx_text(A.chapter); flat = re.sub(r"\s+", " ", txt + " " + " ".join(c for t_ in tabs for r in t_ for c in r))
    figs = re.findall(r"^Figure 10\.(\d+)\.", txt, flags=re.M); check("G.1", "figures are numbered 10.1 to 10.13 in order, one caption for every image; 13 images (3 diagrams and 10 figures; the R check plot is not in the chapter), 24 tables", figs == [str(i) for i in range(1, 14)] and len(d.inline_shapes) == 13 and len(d.tables) == 24, (figs, len(d.inline_shapes), len(d.tables)))
    objs, on = [], False
    for p in d.paragraphs:
        if p.style.name == "Heading 2": on = (p.text == "Learning objectives")
        elif on and p.text.strip(): objs.append(p.text)
    objs = objs[1:]; check("G.2", "fourteen learning objectives, from what host genetics asks to the connection with Chapter 11", len(objs) == 14 and objs[0].startswith("Explain what host genetic") and objs[-1].startswith("Connect host genetics"), len(objs))
    heads = [p.text for p in d.paragraphs if p.style.name == "Heading 2"]; ho = [h for h in heads if h.startswith("Hands-on")]; check("G.3", "nine hands-on sections numbered 1 to 9 in order", [h.split(" ")[1] for h in ho] == [str(i) for i in range(1, 10)], ho)
    want = ["4,006", "3,726", "87,694", "0.110", "110,213", "39,733", "49.7%", "0.8%", "146 cases", "6,916", "89% power", "35 recompute", "116,204", "7.4-fold", "24 independent signals", "625", "5.81", "2.8 × 10⁻⁶", "1,055,453", "11 variants", "0.167", "λ = 1.045", "27,658"]
    miss = [w for w in want if w not in flat]; check("G.4", "the sentences quote the validated numbers", not miss, miss)
    check("G.5", "no stray backslash, image path, header artefact or placeholder in the text", not re.search(r"\\\\|media/|@Emile|TODO|XXX|Nov 9, 2026", flat), re.findall(r"\\\\|media/|@Emile|TODO|XXX|Nov 9, 2026", flat)[:3])
    def tbl(prefix, second=None): return next((t_ for t_ in tabs if t_[0][0].startswith(prefix) and (second is None or t_[0][1].startswith(second))), None)
    t_h = tbl("Test at P"); ok = t_h is not None and [int(r[1].replace(",", "")) for r in t_h[1:]] == list(ht.fail_p_1e6[:8]); check("G.6", "the Hardy-Weinberg table of the chapter equals the result file", ok, None if t_h is None else t_h[1])
    t_p = tbl("Odds ratio"); ok = t_p is not None and [int(r[3].replace(",", "")) for r in t_p[1:]] == list(need["maf_0.2"]) and [int(r[1].replace(",", "")) for r in t_p[1:]] == list(need["maf_0.01"]); check("G.7", "the table of cases needed equals the result file", ok)
    t_g = tbl("Group", "People"); ok = t_g is not None and [r[0] for r in t_g[1:]] == GROUPS and [r[2] for r in t_g[1:]] == [f"{x:.3f}" for x in t3.risk_allele_freq] and [r[3] for r in t_g[1:]] == [f"{100 * x:.1f}%" for x in t3.carriers]; check("G.8", "the 3p21.31 table of the chapter equals the result file", ok, None if t_g is None else t_g[1])
    t_au = tbl("Source"); ok = t_au is not None and t_au[-1][1:] == ["37", "35", "2"]; check("G.9", "the audit table of the chapter ends with 37 figures, 35 agreeing, 2 not", ok, None if t_au is None else t_au[-1])
    import zipfile; z = zipfile.ZipFile(A.chapter); rels = z.read("word/_rels/document.xml.rels").decode(); dois = {u.rstrip(".,;").replace("&amp;", "&") for u in re.findall(r"https://doi\.org/[^\"'<> ]+", rels)}
    check("G.10", "the chapter links 56 distinct DOIs", len(dois) == 56, len(dois))
    if A.doi:
        import urllib.request, urllib.parse; bad_d = []
        for dd_ in sorted(dois):
            doi = dd_.split("doi.org/")[1]
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "mepi-book/1.0 (mailto:none@example.org)"}), timeout=40).read()
            except Exception as e: bad_d.append((doi, str(e)[:30]))
        check("G.11", "every DOI of the chapter resolves in Crossref", not bad_d, bad_d)
    else: skip("G.11", "every DOI of the chapter resolves in Crossref", "run with --doi")
else:
    for k, w in [("G.1", "figure numbers, tables"), ("G.2", "objectives"), ("G.3", "hands-on sections"), ("G.4", "sentences"), ("G.5", "stray characters"), ("G.6", "HWE table"), ("G.7", "power table"), ("G.8", "3p21.31 table"), ("G.9", "audit table"), ("G.10", "DOI count"), ("G.11", "DOI resolution")]: skip(k, w, "run with --chapter")
# ================================================================== H. figures and scripts
figs_ = ["diagram_pipeline_at_a_glance", "diagram_ladder_of_questions", "diagram_structure_stratification_admixture", "fig1_catalog", "fig2_qc_manhattan", "fig2b_r_leads", "fig3_pca", "fig4_stratification_qq", "fig5_finemap", "fig6_qc_structure", "fig7_power_lambda", "fig8_3p21_ancestry", "fig9_admixture_mapping", "fig10_agentic_workflow_diagram"]
fd = next((d_ for d_ in ("../figures", "figures", os.path.join(HERE, "../figures")) if os.path.isdir(d_)), None)
if fd: miss = [f_ for f_ in figs_ if not os.path.exists(os.path.join(fd, f_ + ".png"))]; check("H.1", "all 14 figure files are present (3 diagrams, 10 figures and the R check plot)", not miss, miss)
else: skip("H.1", "figure files", "no figures folder")
cd = next((d_ for d_ in (".", HERE) if os.path.exists(os.path.join(d_, "ex10_9_notes_audit.py"))), None)
if cd:
    import py_compile; badc = []
    for f_ in sorted(glob.glob(os.path.join(cd, "*.py"))):
        try: py_compile.compile(f_, doraise=True)
        except Exception: badc.append(os.path.basename(f_))
    check("H.2", "every Python script of the package compiles", not badc, badc)
out_ = pd.DataFrame(R); out_.to_csv("validation_results_ch10.csv", index=False)
for r in R: print(f"{r['result']:5} {r['id']:8} {r['check'][:150]}" + (f"   [{r['detail']}]" if r["result"] != "PASS" else ""))
print(f"\n{(out_.result == 'PASS').sum()} PASS, {(out_.result == 'FAIL').sum()} FAIL, {(out_.result == 'SKIP').sum()} SKIP"); sys.exit(1 if (out_.result == "FAIL").any() else 0)
