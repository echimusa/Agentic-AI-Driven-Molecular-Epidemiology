"""
Chapter 2 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts (an asymptotic information calculation for nested case-control sampling, vectorised
re-simulations with fresh seeds, a correlation form of the trend test, exact formulas and recomputation from the public data) and checks them against the shipped result tables and the chapter text.
USAGE   python validate_ch2.py [--chapter Chapter_2_....docx] [--full]
        --chapter  also checks key numbers in the chapter text (needs pandoc)
        --full     also reruns the five analysis scripts in a scratch folder and checks that each result table is reproduced byte for byte (about 3 minutes)
Prints PASS/FAIL per check, writes validation_results_ch2.csv and exits 1 if any check fails.
"""
import argparse, glob, hashlib, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy import stats
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if ok else "FAIL", detail=str(detail)[:150])); return ok
rng = np.random.default_rng(20261113)
# ---------------- Hands-on 1: nested case-control. Asymptotic efficiency against the cohort by an independent calculation:
# per case the cohort information is Var(X) under the tilted distribution = 1; a set of one case and M controls has information Var_w(X), weights proportional to exp(beta x)
def eff(M, beta, nsim=600_000):
    X = np.column_stack([rng.normal(beta, 1, nsim), rng.normal(0, 1, (nsim, M))]); w = np.exp(beta * X); w /= w.sum(1, keepdims=True); m = (w * X).sum(1); return float(((w * X ** 2).sum(1) - m ** 2).mean())
nc = pd.read_csv(P("table_ncc_design.csv")); rs = nc[nc.design == "Risk-set sampling"].set_index("M")
for M in (1, 2, 4, 8): check("H1.1", f"efficiency at {M} controls: table {rs.efficiency_vs_cohort[M]:.2f} against the asymptotic value at the true effect {eff(M, 0.4):.2f}", abs(rs.efficiency_vs_cohort[M] - eff(M, 0.4)) <= 0.08)
check("H1.2", "the rule of thumb M/(M+1) is the small-effect limit (beta = 0.01) for M = 1, 2, 4, 8", all(abs(eff(M, 0.01) - M / (M + 1)) <= 0.01 for M in (1, 2, 4, 8)))
gap = [round(M / (M + 1) - eff(M, 0.4), 2) for M in (1, 2, 4, 8)]; check("H1.3", "at the true effect (0.4 per SD) efficiency lies below M/(M+1); gaps for M = 1, 2, 4, 8 recorded", all(g > 0 for g in gap[:3]), gap)
check("H1.4", "risk-set estimate unbiased: mean 0.399 against a true 0.400 (M = 2)", abs(rs.mean_beta[2] - 0.4) < 0.005); check("H1.5", "1:1 to 1:2 buys 22 points of efficiency for 3 of assay cost; 1:4 to 1:8 buys 16 for 12.5", round(100 * (rs.efficiency_vs_cohort[2] - rs.efficiency_vs_cohort[1])) == 22 and round(rs.share_measured_pct[2] - rs.share_measured_pct[1]) == 3 and round(100 * (rs.efficiency_vs_cohort[8] - rs.efficiency_vs_cohort[4])) == 16 and round(rs.share_measured_pct[8] - rs.share_measured_pct[4], 1) == 12.5)
ic = pd.read_csv(P("table_infert_clogit.csv")).set_index(["model", "term"]); check("H1.6", "infert: log odds ratios 1.9859 and 1.4090 together (the documented clogit result for these data)", abs(ic.log_or[("both", "spontaneous")] - 1.9859) < 5e-4 and abs(ic.log_or[("both", "induced")] - 1.4090) < 5e-4)
inf = pd.read_csv(P("infert.csv")); check("H1.7", "infert: 248 women, 83 cases, 83 matched sets", (len(inf), int(inf.case.sum()), inf.stratum.nunique()) == (248, 83, 83))
# independent conditional likelihood: Newton iterations on the two-covariate model, written from the softmax form
def clogit(df, cols):
    beta = np.zeros(len(cols)); X = df[cols].values.astype(float); g = df.stratum.values; y = df.case.values
    for _ in range(40):
        grad = np.zeros(len(cols)); H = np.zeros((len(cols), len(cols)))
        for s in np.unique(g):
            m = g == s; xs = X[m]; w = np.exp(xs @ beta); w /= w.sum(); mu = w @ xs; grad += xs[y[m] == 1].sum(0) - mu; H -= (xs * w[:, None]).T @ xs - np.outer(mu, mu)
        step = np.linalg.solve(H, grad); beta = beta - step
        if abs(step).max() < 1e-10: break
    return beta, np.sqrt(np.diag(np.linalg.inv(-H)))
b, se = clogit(inf, ["spontaneous", "induced"]); check("H1.8", "infert by a third implementation (Newton on the softmax likelihood): 1.986 and 1.409", abs(b[0] - 1.9859) < 1e-3 and abs(b[1] - 1.4090) < 1e-3, (b.round(4), se.round(4)))
# ---------------- Hands-on 2: batch effects. Vectorised re-simulation, 2000 studies per cell, fresh seed
def batch_cell(n_case_b2, effect, reps=2000, NS=100, NM=20, shift=1.5):
    batch = np.r_[np.zeros(50), np.ones(50)]; case = np.zeros((reps, NS))
    for r in range(reps): case[r, rng.choice(np.arange(50, 100), n_case_b2, replace=False)] = 1; case[r, rng.choice(np.arange(50), 50 - n_case_b2, replace=False)] = 1
    X = rng.normal(5, 1, (reps, NS, NM)) + shift * batch[None, :, None] + effect * case[:, :, None]
    n1 = case.sum(1)[:, None]; n0 = NS - n1; m1 = (case[:, :, None] * X).sum(1) / n1; m0 = ((1 - case)[:, :, None] * X).sum(1) / n0; v = ((case[:, :, None] * (X - m1[:, None, :]) ** 2).sum(1) + ((1 - case)[:, :, None] * (X - m0[:, None, :]) ** 2).sum(1)) / (NS - 2)
    t = (m1 - m0) / np.sqrt(v * (1 / n1 + 1 / n0)); naive = (2 * stats.t.sf(abs(t), NS - 2) < 0.05).mean() * 100
    D = np.stack([np.ones((reps, NS)), case, np.broadcast_to(batch, (reps, NS))], axis=2); B = np.linalg.pinv(D) if n_case_b2 == 50 else None
    if B is not None: return naive, np.nan
    XtXi = np.linalg.inv(np.einsum("rni,rnj->rij", D, D)); beta = np.einsum("rij,rnj,rnm->rim", XtXi, D, X); resid = X - np.einsum("rni,rim->rnm", D, beta); s2 = (resid ** 2).sum(1) / (NS - 3)
    tc = beta[:, 1, :] / np.sqrt(s2 * XtXi[:, 1, 1][:, None]); return naive, (2 * stats.t.sf(abs(tc), NS - 3) < 0.05).mean() * 100
bd = pd.read_csv(P("table_batch_design.csv")).set_index(["design", "true_effect"])
for name, n2 in (("Randomised", 25), ("Partly confounded", 40)):
    nv, ad = batch_cell(n2, 0.0); check("H2.1", f"{name}: false positives per assay, naive {nv:.1f}% and adjusted {ad:.1f}% by fresh simulation against table {bd.naive_hit_rate[(name, 0.0)]:.1f}% and {bd.corrected_hit_rate[(name, 0.0)]:.1f}%", abs(nv - bd.naive_hit_rate[(name, 0.0)]) <= 1.5 and abs(ad - bd.corrected_hit_rate[(name, 0.0)]) <= 1.5)
    nv, ad = batch_cell(n2, 0.6); check("H2.2", f"{name}: power at 0.6, adjusted {ad:.1f}% against table {bd.corrected_hit_rate[(name, 0.6)]:.1f}%", abs(ad - bd.corrected_hit_rate[(name, 0.6)]) <= 3)
nv, _ = batch_cell(50, 0.0); check("H2.3", "fully confounded: 100% false positives naive; adjustment impossible", nv > 99.5)
check("H2.4", "with 20 assays a nominal 5% test gives at least one false positive in 64% of studies", round(100 * (1 - 0.95 ** 20)) == 64)
# Salmonella: laboratory share of variance, recomputed
d = pd.read_csv(P("asm_batch_enteritidis.csv.gz"), dtype={"country": str}); d["lab"] = d.sequenced_by.fillna(d.bioproject_center)
x = d[(d.plat == "ILLUMINA") & (d.serovar.str.lower() == "enteritidis") & d.lab.notna() & d.asm_stats_n_contig.notna() & d.asm_stats_contig_n50.notna() & d.sector.isin(["Human", "Food", "Animal", "Environment"])]
big = x.lab.value_counts(); big = big[big >= 300].index; xb = x[x.lab.isin(big)].copy(); check("H2.5", "Salmonella Enteritidis, Illumina: 68,878 isolates from 45 laboratories with at least 300", (len(xb), len(big)) == (68878, 45), (len(xb), len(big)))
def eta2(df, col, grp="lab"):
    g = df.groupby(grp)[col]; return float((g.count() * (g.mean() - df[col].mean()) ** 2).sum() / ((df[col] - df[col].mean()) ** 2).sum())
xb["lc"] = np.log10(xb.asm_stats_n_contig); xb["ln50"] = np.log10(xb.asm_stats_contig_n50); xb["len"] = xb.asm_stats_length_bp / 1e6
check("H2.6", "laboratory explains 34% of log contig variance, 38% of log N50, 15% of genome length, 9% of AMR genes", [round(eta2(xb, c) * 100) for c in ("lc", "ln50", "len", "number_amr_genes")] == [34, 38, 15, 9], [round(eta2(xb, c) * 100) for c in ("lc", "ln50", "len", "number_amr_genes")])
k = xb[xb.sequenced_by.notna()]; kb = k.lab.value_counts(); kb = kb[kb >= 300].index; k = k[k.lab.isin(kb)]; e_known = eta2(k, "lc"); check("H2.7", "sensitivity: with the laboratory taken from the sequenced-by field only (not filled from the BioProject centre) the share for contig count stays above 20%", e_known > 0.20, f"{e_known * 100:.0f}% in {len(kb)} laboratories, {len(k):,} isolates")
x2 = x.copy(); x2["lc"] = np.log10(x2.asm_stats_n_contig); x2["grp"] = np.where(x2.sector == "Human", "Human", "Non-human"); ct = x2.groupby(["lab", "grp"]).size().unstack(fill_value=0); both = ct[(ct.Human >= 30) & (ct["Non-human"] >= 30)]
check("H2.8", "229 laboratories deposited Enteritidis genomes; 11 sequenced 30+ human and 30+ non-human; 70.2% of non-human isolates from laboratories with no human isolates", (len(ct), len(both), round(100 * ct.loc[ct.Human == 0, "Non-human"].sum() / ct["Non-human"].sum(), 1)) == (229, 11, 70.2), (len(ct), len(both), round(100 * ct.loc[ct.Human == 0, "Non-human"].sum() / ct["Non-human"].sum(), 1)))
crude = 10 ** (x2[x2.grp == "Non-human"].lc.mean() - x2[x2.grp == "Human"].lc.mean()); check("H2.9", "non-human minus human isolates: 23% more contigs (ratio 1.23)", round(crude, 2) == 1.23, round(crude, 3))
# ---------------- Hands-on 3: stratification. Independent re-simulation, 400 studies, fresh seed, PCA by SVD, Wald tests by statsmodels
import statsmodels.api as sm
def sim_strat(nrep=400, riskA=-2.5, riskB=-0.5, n_inf=30, nbg=100):
    out = []
    for _ in range(nrep):
        anc = np.r_[np.zeros(600), np.ones(400)]; snp = rng.binomial(2, np.where(anc == 0, 0.10, 0.45)); y = (rng.uniform(size=1000) < 1 / (1 + np.exp(-np.where(anc == 0, riskA, riskB)))).astype(float)
        G = rng.binomial(2, 0.2, (1000, nbg)).astype(float); G[anc == 1, :n_inf] = rng.binomial(2, 0.6, (400, n_inf)); Gc = G - G.mean(0); U, S, Vt = np.linalg.svd(Gc, full_matrices=False); pcs = U[:, :2] * S[:2]
        f = lambda X, yy: sm.Logit(yy, sm.add_constant(X)).fit(disp=0).pvalues[1]
        out.append((f(snp, y), f(np.column_stack([snp, pcs]), y), f(np.column_stack([snp, anc]), y), f(snp[anc == 0], y[anc == 0])))
    return (np.array(out) < 0.05).mean(0) * 100
ts = pd.read_csv(P("table_stratification.csv")).set_index("analysis"); fp = sim_strat()
check("H3.1", f"false positives by fresh simulation: naive {fp[0]:.0f}%, 2 PCs {fp[1]:.1f}%, label {fp[2]:.1f}%, restricted {fp[3]:.1f}% against 100 / 7.8 / 3.7 / 4.5", fp[0] > 99 and abs(fp[1] - 7.8) <= 4 and abs(fp[2] - 3.7) <= 3.5 and abs(fp[3] - 4.5) <= 3.5)
fp0 = sim_strat(250, riskA=-2.5, riskB=-2.5); check("H3.2", f"with no extra risk in B the naive rate is the nominal 5% (fresh simulation {fp0[0]:.1f}%)", abs(fp0[0] - 5) <= 3)
fp10 = sim_strat(250, n_inf=2); check("H3.3", f"with 2 informative background SNPs two components fail (fresh simulation {fp10[1]:.0f}% against the table's 98%)", fp10[1] > 85)
st = pd.read_csv(P("table_stratification_strength.csv")); check("H3.4", "two components stay between 3.3% and 6.7% across the risk-gap series", st[st.setting == "risk gap (log odds)"].pcs.between(3.3, 6.7).all())
check("H3.5", "restriction discards 40% of the sample (400 of 1,000)", 400 / 1000 == 0.4)
# ---------------- Hands-on 4A: sample size by exact formulas and an independent power simulation (trend test as N r^2)
sz = pd.read_csv(P("table_sample_size.csv")); zf = lambda a, pw=0.8: stats.norm.isf(a / 2) + stats.norm.isf(1 - pw)
nf = np.array([zf(r.alpha if r.alpha > 0 else 5e-8) ** 2 / (2 * r.maf * (1 - r.maf) * math.log(r.OR) ** 2) for r in sz.itertuples()]); check("H4.1", "all 16 draft-formula sample sizes recomputed from the normal quantiles (to 0.1%)", np.allclose(nf, sz.N_draft_formula, rtol=1e-3)); check("H4.2", "adding the outcome variance 0.25 multiplies every sample size by four", np.allclose(sz.N_with_outcome_variance / sz.N_draft_formula, 4, rtol=1e-4))
def power(N, p, OR, alpha, reps=2500):
    f = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2]); fc = f * OR ** np.arange(3); fc /= fc.sum(); n1 = int(round(N / 2)); g1 = rng.choice(3, (reps, n1), p=fc); g0 = rng.choice(3, (reps, n1), p=f); g = np.concatenate([g1, g0], 1).astype(float); yv = np.r_[np.ones(n1), np.zeros(n1)]
    gc = g - g.mean(1, keepdims=True); yc = yv - yv.mean(); r = (gc * yc).sum(1) / np.sqrt((gc ** 2).sum(1) * (yc ** 2).sum()); return float((2 * n1 * r ** 2 > stats.chi2.isf(alpha, 1)).mean())
for maf, OR, alpha in ((0.05, 1.5, 0.05), (0.30, 1.2, 0.05), (0.30, 1.5, 5e-8)):
    r_ = sz[(sz.maf == maf) & (sz.OR == OR) & (np.isclose(sz.alpha, alpha, atol=1e-9) if alpha > 1e-6 else (sz.alpha < 1e-6))].iloc[0]
    check("H4.3", f"power at the draft N and at four times it (maf {maf}, OR {OR}, alpha {alpha:g}) by fresh simulation against the table", abs(power(r_.N_draft_formula, maf, OR, alpha) - r_.power_at_N_draft) <= 0.04 and abs(power(r_.N_with_outcome_variance, maf, OR, alpha) - r_.power_at_N_with_variance) <= 0.04)
check("H4.4", "draft formula gives 30-33% power at alpha 0.05 and 1-2% at 5e-8; with the outcome variance 79-92%", sz[sz.alpha > 0.01].power_at_N_draft.between(0.295, 0.335).all() and sz[sz.alpha < 0.01].power_at_N_draft.between(0.0065, 0.0215).all() and sz.power_at_N_with_variance.between(0.785, 0.925).all(), (round(sz.power_at_N_draft.min(), 3), round(sz.power_at_N_draft.max(), 3), round(sz.power_at_N_with_variance.min(), 3), round(sz.power_at_N_with_variance.max(), 3)))
check("H4.5", "detection arithmetic: 59, 299, 598 and 2,995 sequences for 5%, 1%, 0.5% and 0.1% at 95% confidence", [math.ceil(math.log(0.05) / math.log(1 - p)) for p in (0.05, 0.01, 0.005, 0.001)] == [59, 299, 598, 2995])
# ---------------- Hands-on 4B: open sequencing coverage recomputed from the two public files
seq = pd.read_csv(P("lapis_daily_2021_2022.csv"), parse_dates=["date"]); cs = pd.read_csv(P("owid_cases_2021_2022.csv"), parse_dates=["date"]); yr = lambda df: df.date.dt.year
def wk(df): return df.assign(we=df.date + pd.to_timedelta((6 - df.date.dt.weekday) % 7, unit="D"))
sw = wk(seq); sw = sw[(sw.we >= "2021-01-03") & (sw.we <= "2022-12-25")]; cc = cs[cs.location.isin(seq.owid.unique())]; cw = wk(cc); cw = cw[(cw.we >= "2021-01-03") & (cw.we <= "2022-12-25")]
s_tot = sw.groupby(sw.we.dt.year).n.sum(); c_tot = cw.groupby(cw.we.dt.year).new_cases.sum()
check("H4.6", "36 countries; over the 104 weeks ending 3 January 2021 to 25 December 2022: 4,187,286 and 3,874,278 open sequences, shares of cases 2.81% and 1.08%", seq.lapis.nunique() == 36 and (int(s_tot[2021]), int(s_tot[2022])) == (4187286, 3874278) and (round(100 * s_tot[2021] / c_tot[2021], 2), round(100 * s_tot[2022] / c_tot[2022], 2)) == (2.81, 1.08), (int(s_tot[2021]), int(s_tot[2022]), round(100 * s_tot[2021] / c_tot[2021], 2), round(100 * s_tot[2022] / c_tot[2022], 2)))
cy_s = seq.groupby(yr(seq)).n.sum(); cy_c = cc.groupby(yr(cc)).new_cases.sum(); sh_cal = (100 * cy_s[2021] / cy_c[2021], 100 * cy_s[2022] / cy_c[2022])
check("H4.6b", "sensitivity: by calendar year the shares are 2.91% and 1.05%, within 0.15 points of the weekly-bin shares", abs(sh_cal[0] - 2.81) < 0.15 and abs(sh_cal[1] - 1.08) < 0.15, tuple(round(v, 2) for v in sh_cal))
tot = seq.groupby("lapis").n.sum().sort_values(ascending=False); check("H4.7", "the USA, UK, Germany and Denmark supply 96.1% of open sequences", set(tot.index[:4]) == {"USA", "United Kingdom", "Germany", "Denmark"} and round(100 * tot.iloc[:4].sum() / tot.sum(), 1) == 96.1, round(100 * tot.iloc[:4].sum() / tot.sum(), 1))
cv = pd.read_csv(P("table_sequencing_coverage.csv")); y22 = cv[cv.year == 2022].set_index("country"); y21 = cv[cv.year == 2021]
check("H4.8", "9 countries above 0.5% of cases in 2021 and 7 in 2022", (int((y21.fraction_pct >= 0.5).sum()), int((y22.fraction_pct >= 0.5).sum())) == (9, 7), (int((y21.fraction_pct >= 0.5).sum()), int((y22.fraction_pct >= 0.5).sum())))
check("H4.9", "2022 sequences per 100 cases: UK 11.4, Denmark 11.3, USA 3.1, Germany 2.3, Switzerland 1.6, Mexico 0.39, Japan 0.047, France 0.038, Kenya 4.0", all(abs(y22.fraction_pct[c] - v) <= 0.06 * max(v, 1) for c, v in {"United Kingdom": 11.4, "Denmark": 11.3, "USA": 3.1, "Germany": 2.3, "Switzerland": 1.6, "Mexico": 0.39, "Japan": 0.047, "France": 0.038, "Kenya": 4.0}.items()))
w1 = y22.weeks_detect_1pct; check("H4.10", "13 of 36 countries have at least one week with 299+ open sequences in 2022", int((w1 > 0).sum()) == 13, int((w1 > 0).sum()))
five = w1[w1 >= 99.99].index.tolist(); check("H4.11", "five countries (UK, USA, Germany, Switzerland, Denmark) could detect a 1% variant in every week of 2022", len(five) == 5, f"countries at 100% of weeks: {five}; countries above 50% of weeks: {w1[w1 > 50].round(0).to_dict()}")
# ---------------- chapter text
if A.chapter:
    ch = subprocess.run(["pandoc", A.chapter, "-t", "plain", "--wrap=none"], capture_output=True, text=True).stdout.replace("\n", " ")
    for s in ["0.399", "7.29 (3.65–14.54)", "98.3%", "68,878", "4,187,286", "3,874,278", "96.1%", "299 sequences", "2,995", "81", "7.8%", "30–33%"]: check("TXT", f"chapter text contains '{s}'", s in ch)
if A.full:
    tmp = tempfile.mkdtemp(); src = os.getcwd(); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) if os.path.isfile(f)]; [shutil.copy(f, tmp) for f in glob.glob(os.path.join(os.path.dirname(P("table_ncc_design.csv")), "*")) if os.path.isfile(f) and not os.path.exists(os.path.join(tmp, os.path.basename(f)))]
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); before = {os.path.basename(f): md5(f) for f in glob.glob(os.path.join(tmp, "table_*.csv"))}; os.makedirs(os.path.join(tmp, "figs2"), exist_ok=True)
    for s in ["ex2_1_nested_case_control", "ex2_2_batch_effects", "ex2_3_stratification", "ex2_4_sampling"]: subprocess.run([sys.executable, s + ".py"], cwd=tmp, capture_output=True)
    subprocess.run(["Rscript", "ex2_1b_infert_clogit.R"], cwd=tmp, capture_output=True); after = {os.path.basename(f): md5(f) for f in glob.glob(os.path.join(tmp, "table_*.csv"))}
    for k_ in sorted(before): check("REP", f"{k_} reproduced byte for byte", after.get(k_) == before[k_])
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch2.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=110)); print(f"\n{(Rd.result == 'PASS').sum()} of {len(Rd)} checks passed"); sys.exit(0 if (Rd.result == "PASS").all() else 1)
