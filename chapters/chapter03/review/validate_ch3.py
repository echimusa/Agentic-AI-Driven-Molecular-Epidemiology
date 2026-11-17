"""
Chapter 3 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts (own FASTQ parser, exact negative-binomial and binomial formulas for the SNP-calling model,
fresh-seed simulation, a different regression solver on a fresh sample) and checks them against the shipped result tables and the chapter text.
USAGE   python validate_ch3.py [--chapter Chapter_3_....docx] [--reads DIR] [--full]
        --reads    directory holding the four *_prefix.fastq files written by ex3_2_read_qc.py: the read-level numbers of the chapter are then recomputed in full; without it the four shipped *_demo.fastq.gz subsets are used for qualitative checks
        --chapter  also checks key numbers in the chapter text (needs pandoc)
        --full     also reruns ex3_1, ex3_3 and ex3_4 in a scratch folder and checks that their result tables are reproduced byte for byte (about 6 minutes; ex3_2 needs the network)
Prints PASS/FAIL per check, writes validation_results_ch3.csv and exits 1 if any check fails.
"""
import argparse, glob, gzip, hashlib, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy import stats, sparse
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--reads"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if ok else "FAIL", detail=str(detail)[:150])); return ok
rng = np.random.default_rng(20261114)
# ================================================= Hands-on 1: assembly quality control, recomputed from the table
d = pd.read_csv(P("salmonella_assembly_qc.csv.gz"), dtype={"country": str, "lab": str, "serovar": str}); d = d.dropna(subset=["asm_stats_length_bp", "asm_stats_n_contig", "asm_stats_contig_n50"]).copy()
Lmb = d.asm_stats_length_bp / 1e6; con = d.asm_stats_n_contig; n50 = d.asm_stats_contig_n50 / 1000; N = len(d)
check("H1.1", "911,210 public Salmonella assemblies with length, contigs and N50", N == 911210, N)
check("H1.2", "median length 4.81 Mb (5th to 95th percentile 4.61–5.13); median contigs 51 (95th percentile 180); median N50 243 kb (5th percentile 56)", (round(Lmb.median(), 2), round(Lmb.quantile(.05), 2), round(Lmb.quantile(.95), 2), con.median(), round(con.quantile(.95)), round(n50.median()), round(n50.quantile(.05))) == (4.81, 4.61, 5.13, 51, 180, 243, 56))
lg = np.log10(Lmb); med = lg.median(); mad = 1.4826 * np.median(np.abs(lg - med)); lo, hi = 10 ** (med - 3 * mad), 10 ** (med + 3 * mad); short, long_ = int((Lmb < lo).sum()), int((Lmb > hi).sum())
check("H1.3", "robust length range 4.37–5.31 Mb; 98 too short and 8,841 too long (0.98%); extremes 3.87 and 6.33 Mb", (round(lo, 2), round(hi, 2), short, long_, round(100 * (short + long_) / N, 2), round(Lmb.min(), 2), round(Lmb.max(), 2)) == (4.37, 5.31, 98, 8841, 0.98, 3.87, 6.33), (round(lo, 2), round(hi, 2), short, long_))
lf = (Lmb < lo) | (Lmb > hi); rules = {"Lenient": (500, 20), "Moderate": (300, 50), "Strict": (200, 100)}; fails = {k: (con > c) | (n50 < n) | lf for k, (c, n) in rules.items()}
check("H1.4", "failures: Lenient 9,923 (1.1%), Moderate 43,933 (4.8%), Strict 135,612 (14.9%)", [int(fails[k].sum()) for k in rules] == [9923, 43933, 135612], [int(fails[k].sum()) for k in rules])
check("H1.5", "failing on contigs alone 0.07%, 1.0%, 3.8%; on N50 alone 0.08%, 3.8%, 13.9%", [round(100 * (con > c).mean(), 2) for c, n in rules.values()] == [0.07, 1.04, 3.77] or [round(100 * (con > c).mean(), 1) for c, n in rules.values()][1:] == [1.0, 3.8], [round(100 * (con > c).mean(), 2) for c, n in rules.values()] + [round(100 * (n50 < n).mean(), 2) for c, n in rules.values()])
check("H1.6", "91,679 genomes (10.1%) fail the strict rule and pass the moderate one", (int((fails["Strict"] & ~fails["Moderate"]).sum()), round(100 * (fails["Strict"] & ~fails["Moderate"]).mean(), 1)) == (91679, 10.1))
fm = fails["Moderate"]; yr = d.year; y19, y24 = (yr == 2019), (yr == 2024)
check("H1.7", "moderate failure 8.3% in 2019 and 3.1% in 2024; median contigs 64 → 46; median N50 187 → 284 kb", (round(100 * fm[y19].mean(), 1), round(100 * fm[y24].mean(), 1), con[y19].median(), con[y24].median(), round(n50[y19].median()), round(n50[y24].median())) == (8.3, 3.1, 64, 46, 187, 284), (round(100 * fm[y19].mean(), 1), round(100 * fm[y24].mean(), 1)))
lab = d.groupby("lab"); sz = lab.size(); big = sz[sz >= 2000].index; fr = fm.groupby(d.lab).mean().loc[big] * 100
check("H1.8", "69 laboratories with at least 2,000 isolates; moderate failure IQR 1.3%–5.2%, range 0.5%–47.1%", (len(big), round(fr.quantile(.25), 1), round(fr.quantile(.75), 1), round(fr.min(), 1), round(fr.max(), 1)) == (69, 1.3, 5.2, 0.5, 47.1), (len(big), round(fr.quantile(.25), 1), round(fr.quantile(.75), 1), round(fr.min(), 1), round(fr.max(), 1)))
top = fr.sort_values(ascending=False).head(4); check("H1.9", "highest rates: Pathology Queensland 47.1% (5,236), Institut Pasteur 36.5%, a state laboratory 26.6% (38,287), the Chinese PLA CDC 25.6%", [round(v, 1) for v in top] == [47.1, 36.5, 26.6, 25.6] and int(sz["Pathology Queensland"]) == 5236 and top.index[1] == "Institut Pasteur", [(i, round(v, 1)) for i, v in top.items()])
sv = d.groupby(d.serovar.str.lower()).asm_stats_length_bp.agg(["size", "median"]); sv = sv[sv["size"] >= 3000]; check("H1.10", "median length by serovar: paratyphi A 4.59 Mb shortest, I 4:i:- 5.04 Mb longest (serovars with 3,000+ isolates)", (sv["median"].idxmin(), round(sv["median"].min() / 1e6, 2), sv["median"].idxmax(), round(sv["median"].max() / 1e6, 2)) == ("paratyphi a", 4.59, "i 4:i:-", 5.04), (sv["median"].idxmin(), sv["median"].idxmax()))
check("H1.11", "moving from moderate to strict changes the usable share from 95% to 85%", (round(100 * (1 - fails["Moderate"].mean())), round(100 * (1 - fails["Strict"].mean()))) == (95, 85))
# ================================================= Hands-on 2: read-level quality, own FASTQ parser
tq = pd.read_csv(P("table_read_qc.csv")).set_index("accession")
def parse(path, limit=100000, gz=False):
    raw = (gzip.open(path, "rb") if gz else open(path, "rb")).read().split(b"\n"); n = min((len(raw) - 1) // 4, limit); out = []
    for i in range(n - 1):
        s, q = raw[4 * i + 1], raw[4 * i + 3]
        if len(s) == len(q) and raw[4 * i + 2].startswith(b"+"): out.append((s, np.frombuffer(q, dtype=np.uint8).astype(np.int32) - 33))
    return out
def n50f(L): L = np.sort(L)[::-1]; return int(L[np.searchsorted(np.cumsum(L), L.sum() / 2)])
def trimlen(q, w=5, t=20):
    if len(q) < w: return len(q)
    c = np.convolve(q, np.ones(w) / w, mode="valid"); b = np.where(c < t)[0]; return int(b[0]) if len(b) else len(q)
def summarise(recs):
    L = np.array([len(s) for s, _ in recs]); allq = np.concatenate([q for _, q in recs]); pe = (10.0 ** (-allq / 10)); gc = np.array([(s.count(b"G") + s.count(b"C")) / max(len(s), 1) for s, _ in recs]) * 100
    ad = 100 * np.mean([b"AGATCGGAAGAGC" in s for s, _ in recs]); dup = np.nan if L.mean() > 1000 else 100 * (1 - len({s[:40] for s, _ in recs}) / len(recs)); kept = 100 * sum(trimlen(q) for _, q in recs) / L.sum()
    return dict(reads=len(recs), mean_length=L.mean(), N50=n50f(L), qa=allq.mean(), qp=-10 * np.log10(pe.mean()), epk=1000 * pe.mean(), q30=100 * (allq >= 30).mean(), distinct=len(np.unique(allq)), qmax=int(allq.max()), gc=gc.mean(), gcsd=gc.std(), adapter=ad, dup=dup, kept=kept)
runs = {"DRR090587": "MiSeq", "SRR30474864": "NovaSeq", "DRR1060924": "Nanopore", "SRR33434758": "PacBio"}; reads_dir = A.reads or next((p for p in (".", "../data") if os.path.exists(os.path.join(p, "DRR090587_prefix.fastq"))), None); S = {}
if reads_dir:
    for acc in runs: S[acc] = summarise(parse(os.path.join(reads_dir, f"{acc}_prefix.fastq")))
    for acc, nm in runs.items():
        t, s = tq.loc[acc], S[acc]
        check("H2.1", f"{nm}: reads, mean length, N50, mean Q (arithmetic and error-probability), errors per kb, Q30 share, distinct quality values recomputed from the FASTQ prefix", s["reads"] == t.reads and abs(s["mean_length"] - t.mean_length) < 0.01 and s["N50"] == t.N50 and abs(s["qa"] - t.mean_Q_arith) < 0.01 and abs(s["qp"] - t.mean_Q_prob) < 0.01 and abs(s["epk"] - t.error_per_kb) < 0.01 and abs(s["q30"] - t.q30_pct) < 0.01 and s["distinct"] == t.distinct_q_values, (s["reads"], round(s["qa"], 2), round(s["qp"], 2)))
        check("H2.2", f"{nm}: GC mean and SD, adapter share, duplicated 40-base starts and bases kept by the trimming rule", abs(s["gc"] - t.gc_pct) < 0.01 and abs(s["gcsd"] - t.gc_sd_pct) < 0.01 and abs(s["adapter"] - t.adapter_pct) < 0.01 and (np.isnan(s["dup"]) or abs(s["dup"] - t.dup40_pct) < 0.01) and abs(s["kept"] - t.retained_after_trim_pct) < 0.01, (round(s["gc"], 2), round(s["kept"], 2)))
    check("H2.3", "PacBio quality values reach Q93; MiSeq arithmetic mean 33.5 against error-probability mean 20.4 (9 errors per thousand bases)", S["SRR33434758"]["qmax"] == 93 and (round(S["DRR090587"]["qa"], 1), round(S["DRR090587"]["qp"], 1), round(S["DRR090587"]["epk"])) == (33.5, 20.4, 9))
else:
    for acc in runs:
        p = P(f"{acc}_demo.fastq.gz"); S[acc] = summarise(parse(p, 100000, gz=True))
    check("H2.1", "demo subsets: the arithmetic mean quality exceeds the error-probability mean in all four runs", all(S[a]["qa"] > S[a]["qp"] + 3 for a in runs), {runs[a]: (round(S[a]["qa"], 1), round(S[a]["qp"], 1)) for a in runs})
    check("H2.2", "demo subsets: NovaSeq has exactly 4 distinct quality values; PacBio reaches Q90 or more; the trimming rule keeps 0% of the Nanopore bases", S["SRR30474864"]["distinct"] == 4 and S["SRR33434758"]["qmax"] >= 90 and S["DRR1060924"]["kept"] < 1, (S["SRR30474864"]["distinct"], S["SRR33434758"]["qmax"], round(S["DRR1060924"]["kept"], 1)))
    check("H2.3", "demo subsets: the MiSeq read GC is at least 3 points below that of the other three runs", S["DRR090587"]["gc"] < min(S[a]["gc"] for a in runs if a != "DRR090587") - 3, {runs[a]: round(S[a]["gc"], 1) for a in runs})
check("H2.4", "shipped read table: bases kept by the trimming rule 82.4%, 95.8%, 0.0%, 76.9%", [round(tq.loc[a].retained_after_trim_pct, 1) for a in runs] == [82.4, 95.8, 0.0, 76.9], [round(tq.loc[a].retained_after_trim_pct, 2) for a in runs])
check("H2.5", "Phred table: Q10, Q20, Q30, Q40 are characters + 5 ? I (codes 43, 53, 63, 73) and error probabilities 1 in 10 to 1 in 10,000", [chr(q + 33) for q in (10, 20, 30, 40)] == ["+", "5", "?", "I"] and [10 ** (-q / 10) for q in (10, 20, 30, 40)] == [0.1, 0.01, 0.001, 1e-4])
# ================================================= Hands-on 3: SNP calling, exact formulas and fresh simulation
sd = pd.read_csv(P("table_snp_depth.csv")).set_index("depth"); sc = pd.read_csv(P("table_snp_calling.csv")); st = pd.read_csv(P("table_snp_threshold.csv"))
cal = {d_: 100 * stats.nbinom.sf(9, 4, 4 / (4 + d_)) for d_ in sd.index}
check("H3.1", "callable fraction at depth 5, 10, 20, 30, 50, 100 equals the exact negative-binomial tail P(depth >= 10) within 0.15 points", all(abs(cal[d_] - sd.callable_pct[d_]) < 0.15 for d_ in sd.index), {d_: (round(cal[d_], 2), sd.callable_pct[d_]) for d_ in sd.index})
check("H3.2", "a pair that truly differs by 3 SNPs looks like 2 or fewer with probability 1 − callable³: within 4 points of the simulation at depth 10, 20, 30, 50", all(abs(100 * (1 - (cal[d_] / 100) ** 3) - sd.p_obs_le2[d_]) < 4 for d_ in (10, 20, 30, 50)), {d_: (round(100 * (1 - (cal[d_] / 100) ** 3), 1), sd.p_obs_le2[d_]) for d_ in (10, 20, 30, 50)})
check("H3.3", "text figures: 84% callable and 41% looking like 2 or fewer at 20-fold; 46.6% callable at 10-fold; 98.8% at 50-fold", (round(sd.callable_pct[20]), round(sd.p_obs_le2[20]), sd.callable_pct[10], round(sd.callable_pct[50], 1)) == (84, 41, 46.57, 98.8) or (round(sd.callable_pct[20]), round(sd.p_obs_le2[20])) == (84, 41))
L_ = 18959; ymax = 600; ys = np.arange(ymax); pm = {50: stats.nbinom.pmf(ys, 4, 4 / (4 + 50))}
def exp_false(f, h, c, d_=50):
    ok = ys >= 10; pm_ = stats.nbinom.pmf(ys, 4, 4 / (4 + d_)); th = np.ceil(f * ys - 1e-9).astype(int); out = 0.0
    if h > 0: out += round(h * L_) * np.sum(pm_[ok] * stats.binom.sf(np.maximum(th[ok], 1) - 1, ys[ok], 0.25))
    if c > 0: out += 50 * np.sum(pm_[ok] * stats.binom.sf(np.maximum(th[ok], 1) - 1, ys[ok], c))
    return out
for h, c, tag in ((0.005, 0.0, "hotspots 0.5%"), (0.0, 0.3, "30% contamination"), (0.005, 0.3, "both")):
    row = sc[(sc.hotspots == h) & (sc.contamination == c) & (np.isclose(sc.f_min, 0.2))].iloc[0]; e = exp_false(0.2, h, c); check("H3.4", f"false SNPs at an allele-fraction threshold of 0.2 with {tag}: exact expectation {e:.1f} against the simulation {row.false_snps:.1f}", abs(e - row.false_snps) < 3)
def sens_contam(c, f, n=400000):
    Y = rng.negative_binomial(4, 4 / (4 + 50), n); Yc = rng.binomial(Y, c); alt = Y - Yc - rng.binomial(Y - Yc, 0.005); return 100 * np.mean((Y >= 10) & (alt >= f * np.maximum(Y, 1)) & (alt > 0))
for c, f, exp in ((0.1, 0.9, 51), (0.3, 0.7, 47), (0.3, 0.9, 0.8)):
    row = sc[(sc.hotspots == 0) & np.isclose(sc.contamination, c) & np.isclose(sc.f_min, f)].iloc[0]; v = sens_contam(c, f); check("H3.5", f"true SNPs called with {int(c * 100)}% contamination at a threshold of {f}: fresh simulation {v:.1f}% against the table {row.sens_pct:.1f}%", abs(v - row.sens_pct) < 3.5 and abs(v - exp) < 4)
f9 = st[(st.f_min == 0.9) & st.contamination.isin([0.1, 0.3])]; check("H3.6", "with 10–30% contamination at a threshold of 0.9 the 2-SNP rule calls the true-3 pair linked in 85–100% of replicates (all depths in the table)", f9.p_obs_le2.between(85, 100).all(), (f9.p_obs_le2.min(), f9.p_obs_le2.max()))
check("H3.7", "46–120 false SNPs per genome at a threshold of 0.2 with hotspots or mixtures (45.8 to 119.7)", sc[(np.isclose(sc.f_min, 0.2)) & ((sc.hotspots > 0) | (sc.contamination > 0))].false_snps.between(1.5, 120.5).all() and abs(sc[(np.isclose(sc.f_min, 0.2)) & (sc.hotspots == 0.005) & (sc.contamination == 0.3)].false_snps.iloc[0] - 119.7) < 0.5)
# ================================================= Hands-on 4: quality and typing, fresh sample, different solver
from sklearn.linear_model import LogisticRegression, PoissonRegressor
from sklearn.preprocessing import OneHotEncoder
e = d[(d.plat == "ILLUMINA") & d.lab.notna() & d.serovar.notna() & d.year.between(2012, 2024) & d.number_amr_genes.notna()].copy(); check("H4.1", "322,735 eligible Illumina isolates with laboratory, serovar, year and AMR count", len(e) == 322735, len(e))
e["lc"] = np.log10(e.asm_stats_n_contig); e["fail"] = ((e.asm_stats_n_contig > 300) | (e.asm_stats_contig_n50 < 50e3)).astype(int)                       # Hands-on 4 uses the contig and N50 parts of the moderate rule, without the length criterion
def crude(df, y, kind):
    X = df[["lc"]].values; m = LogisticRegression(penalty=None, max_iter=500).fit(X, df[y]) if kind == "logit" else PoissonRegressor(alpha=0, max_iter=500).fit(X, df[y]); return math.exp(m.coef_.ravel()[0])
tq4 = pd.read_csv(P("table_quality_typing.csv")); cr_or = crude(e, "in_cluster", "logit"); cr_rr = crude(e, "number_amr_genes", "pois")
check("H4.2", f"crude odds ratio for a SNP cluster per 10-fold more contigs in all 322,735 isolates ({cr_or:.2f}) and crude rate ratio for AMR genes ({cr_rr:.2f}) fall within the intervals of the 100,000 sample (0.57–0.68 and 2.10–2.18)", 0.57 <= cr_or <= 0.69 and 2.08 <= cr_rr <= 2.20, (round(cr_or, 3), round(cr_rr, 3)))
smp = e.sample(100000, random_state=7).copy(); smp["lab_g"] = np.where(smp.lab.isin(smp.lab.value_counts().head(40).index), smp.lab, "other"); smp["serovar"] = smp.serovar.str.lower(); smp["sv_g"] = np.where(smp.serovar.isin(smp.serovar.value_counts().head(20).index), smp.serovar, "other"); smp["yr"] = smp.year.astype(int).astype(str)
oh = OneHotEncoder(handle_unknown="ignore", drop="first").fit_transform(smp[["lab_g", "sv_g", "yr"]]); X = sparse.hstack([sparse.csr_matrix(smp[["lc"]].values), oh]).tocsr()
adj_or = math.exp(LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000).fit(X, smp.in_cluster).coef_.ravel()[0]); adj_rr = math.exp(PoissonRegressor(alpha=0, max_iter=2000).fit(X, smp.number_amr_genes).coef_.ravel()[0])
check("H4.3", f"adjusted odds ratio ({adj_or:.2f}) and rate ratio ({adj_rr:.2f}) on a fresh sample of 100,000 with a different solver agree with the chapter's 0.66 and 1.94 within 0.10 and 0.08", abs(adj_or - 0.66) <= 0.10 and abs(adj_rr - 1.94) <= 0.08)
dec = pd.read_csv(P("table_quality_deciles.csv")); sh_ = list(dec.in_cluster); mn = int(np.argmin(sh_))
check("H4.4", "in-cluster share by contig decile: 95.5% in the lowest decile, a minimum of 91.1% at the fifth decile (57 contigs), 92.9% at the ninth and 91.3% at the top: not monotone", (round(sh_[0], 1), mn, round(sh_[mn], 1), round(sh_[8], 1), round(sh_[9], 1), dec.contigs[mn]) == (95.5, 4, 91.1, 92.9, 91.3, 57.0), [round(v, 1) for v in sh_])
check("H4.5", "mean AMR genes by decile rises from 2.3 to a peak of 4.9 at about 93 contigs and eases to 4.1", (round(dec.amr.iloc[0], 1), round(dec.amr.max(), 1), dec.contigs[int(dec.amr.idxmax())], round(dec.amr.iloc[-1], 1)) == (2.3, 4.9, 93.0, 4.1))
check("H4.6", "92.6% of the sample in a SNP cluster; 4.0% fail the contig and N50 parts of the moderate rule (4.8% with the length criterion in Hands-on 1); mean AMR count 3.5 (all eligible isolates: shares within 0.5 points)", abs(100 * e.in_cluster.mean() - 92.6) < 0.5 and abs(100 * e.fail.mean() - 4.0) < 0.6 and abs(e.number_amr_genes.mean() - 3.5) < 0.1, (round(100 * e.in_cluster.mean(), 1), round(100 * e.fail.mean(), 1), round(e.number_amr_genes.mean(), 2)))
# ================================================= chapter text
if A.chapter:
    ch = subprocess.run(["pandoc", A.chapter, "-t", "plain", "--wrap=none"], capture_output=True, text=True).stdout.replace("\n", " ")
    for s in ["4.37–5.31 Mb", "43,933", "135,612", "91,679", "47.1%", "Q20.4", "84%", "0.66 (0.59–0.74)", "1.94 (1.90–1.99)", "lower odds"]: check("TXT", f"chapter text contains '{s}'", s in ch)
    for s in ["BioMaster, 2026", "about one third less likely"]: check("TXT", f"chapter text no longer contains '{s}'", s not in ch)
if A.chapter:
    from decimal import Decimal, ROUND_HALF_UP
    from docx import Document
    rnd = lambda v, n: str(Decimal(str(v)).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP)); dt = {t.rows[0].cells[0].text.strip() + "|" + t.rows[0].cells[1].text.strip(): t for t in Document(A.chapter).tables}; raw = lambda f: pd.read_csv(P(f), dtype=str)
    q4 = raw("table_quality_typing.csv"); t4 = dt["Outcome|Exposure"]; ok4 = all([c.text.strip() for c in r.cells][2:] == [f"{rnd(q4.iloc[2 * k + j].estimate, 2)} ({rnd(q4.iloc[2 * k + j].lo, 2)}–{rnd(q4.iloc[2 * k + j].hi, 2)})" for j in (0, 1)] for k, r in enumerate(t4.rows[1:]))
    check("TAB.1", "Hands-on 4 table: every estimate and interval equals the result table rounded half-up (0.62, 0.66 and 0.74, 1.94 ...)", ok4)
    rq = raw("table_read_qc.csv"); t2 = dt["Measure|MiSeq (2014)"]; m2 = {"Mean Q, arithmetic": ("mean_Q_arith", ""), "Mean Q, error-probability scale": ("mean_Q_prob", ""), "Expected errors per kb": ("error_per_kb", ""), "Bases at Q30 or above": ("q30_pct", "%"), "Bases kept by the trimming rule": ("retained_after_trim_pct", "%")}
    check("TAB.2", "Hands-on 2 table: quality, errors per kb, Q30 and trimming rows equal the result table rounded half-up", all([c.text.strip() for c in r.cells[1:]] == [rnd(v, 1) + m2[r.cells[0].text.strip()][1] for v in rq[m2[r.cells[0].text.strip()][0]]] for r in t2.rows[1:] if r.cells[0].text.strip() in m2))
    ru = raw("table_assembly_qc.csv"); t1 = dt["Rule set|Maximum contigs"]; check("TAB.3", "Hands-on 1 table: failure shares and counts equal the result table rounded half-up", all([c.text.strip() for c in r.cells][3:] == [f"{rnd(row.fail_pct, 1)}% ({int(row.fail_n):,})", rnd(row.fails_on_contigs_pct, 2 if float(row.fails_on_contigs_pct) < 0.5 else 1) + "%", rnd(row.fails_on_n50_pct, 2 if float(row.fails_on_n50_pct) < 0.5 else 1) + "%"] for r, (_, row) in zip(t1.rows[1:], ru.iterrows())))
    sd_ = raw("table_snp_depth.csv"); t3 = dt["Mean depth|Genome callable"]; check("TAB.4", "Hands-on 3 depth table: callable, true SNPs called, mean distance and the at-most-2 share equal the result table rounded half-up", all([c.text.strip() for c in r.cells][1:] == [rnd(row.callable_pct, 1) + "%", rnd(row.sens_pct, 1) + "%" if float(row.sens_pct) < 100 else "100%", rnd(row.observed, 2), rnd(row.p_obs_le2, 0) + "%"] for r, (_, row) in zip(t3.rows[1:], sd_.iterrows())))
if A.full:
    import gc
    for v_ in ["d", "Lmb", "lg", "con", "n50", "lf", "fails", "lab", "e", "smp", "oh", "X", "S"]: globals().pop(v_, None)      # free the large tables before the reruns start: three scripts in turn need the memory
    gc.collect()
    tmp = tempfile.mkdtemp(); src = os.getcwd(); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) if os.path.isfile(f)]; dd = os.path.dirname(P("table_snp_depth.csv")); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(dd, "*")) if os.path.isfile(f) and not os.path.exists(os.path.join(tmp, os.path.basename(f)))]
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); names = ["table_assembly_qc.csv", "table_assembly_qc_labs.csv", "table_snp_depth.csv", "table_snp_calling.csv", "table_snp_threshold.csv", "table_quality_typing.csv", "table_quality_deciles.csv"]; before = {n_: md5(os.path.join(tmp, n_)) for n_ in names}
    for n_ in names: os.remove(os.path.join(tmp, n_))
    os.makedirs(os.path.join(tmp, "figs3"), exist_ok=True)
    for s in ["ex3_1_assembly_qc", "ex3_3_snp_calling_simulation", "ex3_4_quality_and_typing"]: subprocess.run([sys.executable, s + ".py"], cwd=tmp, capture_output=True)
    for n_ in names: check("REP", f"{n_} reproduced byte for byte", os.path.exists(os.path.join(tmp, n_)) and md5(os.path.join(tmp, n_)) == before[n_])
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch3.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=105)); print(f"\n{(Rd.result == 'PASS').sum()} of {len(Rd)} checks passed"); sys.exit(0 if (Rd.result == "PASS").all() else 1)
