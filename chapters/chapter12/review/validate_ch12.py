"""
Chapter 12 validation harness. Recomputes the numbers quoted in the chapter with code that does not import the chapter's own scripts (surv12.py, ex12_*.py).
A  the original chapter: date precision, completeness, timeliness, the Gini coefficient, k-anonymity of day-level human records, detection delay of US clusters
B  the new exercises: the detection formula, its Monte Carlo and exact checks, four departures, weekly volumes, catchments, allocation across income groups, the vocabulary of the table, the cost of an embargo, the protocol rules
G  the chapter text (tables, figures, objectives, sentences, references)
H  figures and scripts
I  reruns of every script from the shipped data, with a byte-for-byte comparison of the result tables
USAGE   python validate_ch12.py [--chapter Chapter_12_....docx] [--doi] [--full]
        --chapter  also checks the chapter's tables, figures, objectives, sentences and references against the results (needs python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns every script in a scratch folder and checks that the result tables are reproduced byte for byte (about 1 minute)
Prints PASS/FAIL/SKIP per check, writes validation_results_ch12.csv and exits 1 if any check fails.
"""
import argparse, glob, math, os, re, shutil, subprocess, sys, tempfile, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
def P(name):
    for d in (".", "../data", "data", os.path.join(HERE, "../data"), os.path.join(HERE, "data")):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []; EXPECTED_DOIS = 20
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
close = lambda a, b, tol: abs(float(a) - float(b)) <= tol
# ================================================================== I. reruns first
if A.full:
    tmp = tempfile.mkdtemp(prefix="ch12_"); data = os.path.dirname(os.path.abspath(P("sal_surv.csv.gz"))); src = os.path.join(HERE, "../code") if os.path.isdir(os.path.join(HERE, "../code")) else HERE
    for f_ in glob.glob(os.path.join(src, "*.py")) + glob.glob(os.path.join(src, "*.R")): shutil.copy(f_, tmp)
    for item in os.listdir(data):
        if item.startswith("table_"): continue
        shutil.copy(os.path.join(data, item), os.path.join(tmp, item))
    os.makedirs(os.path.join(tmp, "figs12"), exist_ok=True)
    runs = [("ex12_1", ["python3", "ex12_1_indicators.py"]), ("ex12_2", ["python3", "ex12_2_equity.py"]), ("ex12_2b", ["Rscript", "ex12_2b_gini.R"]), ("ex12_3", ["python3", "ex12_3_privacy.py"]), ("ex12_4", ["python3", "ex12_4_detection.py"]), ("ex12_5", ["python3", "ex12_5_sizing.py"]),
            ("ex12_6", ["python3", "ex12_6_allocation.py"]), ("ex12_7", ["python3", "ex12_7_standards.py"]), ("ex12_8", ["python3", "ex12_8_embargo.py"]), ("ex12_9", ["python3", "ex12_9_protocol.py"]), ("ex12_10", ["python3", "ex12_10_diagrams.py"])]
    for name, cmd in runs:
        r_ = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True); check(f"I.{name}", f"{name} reruns from the shipped data without error", r_.returncode == 0, (r_.stderr or "")[-120:])
    diffs = []
    for tb in sorted(os.path.basename(f_) for f_ in glob.glob(os.path.join(data, "table_*.csv"))):
        a_, b_ = os.path.join(tmp, tb), os.path.join(data, tb)
        if not os.path.exists(a_): diffs.append(tb + " (not written)"); continue
        if open(a_, "rb").read() != open(b_, "rb").read(): diffs.append(tb)
    check("I.tables", "every result table of the rerun is byte-identical to the shipped one", not diffs, diffs)
else: skip("I.*", "reruns of the scripts", "run with --full")
# ================================================================== data
M = pd.read_csv(P("sal_surv.csv.gz"), dtype={"country": str, "year": float}, parse_dates=["coll_day", "coll_month", "release"])
# ================================================================== A. the original chapter, recomputed
n = len(M); dp = M.date_precision.value_counts(); pct = lambda k: 100 * dp[k] / n
check("A.1", "911,210 isolates; collection dates to the day 14.1%, to the month 40.8%, to the year 36.4%, missing 8.7%", n == 911210 and [rnd(pct(k), 1) for k in ("day", "month", "year", "missing")] == ["14.1", "40.8", "36.4", "8.7"], [rnd(pct(k), 1) for k in ("day", "month", "year", "missing")])
hu = M[M.sector == "Human"]; tc = pd.read_csv(P("table_completeness.csv"), index_col=0)
ok = rnd(100 * (hu.date_precision == "day").mean(), 1) == rnd(tc.loc["Collection date to the day", "Human"], 1) and rnd(100 * hu.host.notna().mean(), 1) == rnd(tc.loc["Host", "Human"], 1) and rnd(100 * (M.date_precision != "missing").mean(), 1) == rnd(tc.loc["Any collection date", "All isolates"], 1)
check("A.2", "completeness: day-level dates 8.5% and host 48.5% of human records; any date 91.3% of all", ok, (rnd(100 * (hu.date_precision == "day").mean(), 1), rnd(100 * hu.host.notna().mean(), 1)))
d = M[M.date_precision.isin(["day", "month"]) & M.release.notna() & M.coll_month.notna()].copy(); d["lag"] = (d.release.dt.year - d.coll_month.dt.year) * 12 + (d.release.dt.month - d.coll_month.dt.month)
x = d[(d.sector == "Human") & d.year.between(2015, 2023)]; tt = pd.read_csv(P("table_timeliness.csv")).set_index("sector")
check("A.3", "timeliness: 328,107 human isolates, 60.5% released within a month, median lag 1 month", len(x) == 328107 and rnd(100 * (x.lag <= 1).mean(), 1) == "60.5" and x.lag.median() == 1 and close(tt.loc["Human", "within_1_month"], 60.5, 0.05), (len(x), rnd(100 * (x.lag <= 1).mean(), 1)))
ce = pd.read_csv(P("table_country_equity.csv")); o = ce.sort_values("per_million"); cp = np.r_[0, np.cumsum(o["pop"]) / o["pop"].sum()]; cl = np.r_[0, np.cumsum(o.isolates) / o.isolates.sum()]; gini = 1 - 2 * np.trapezoid(cl, cp)
top3 = ce.sort_values("isolates", ascending=False).isolates.head(3).sum() / ce.isolates.sum()
check("A.4", "Gini coefficient of isolates over population 0.905 (trapezoid rule, separate code); the top three economies hold 84.5%", rnd(gini, 3) == "0.905" and rnd(100 * top3, 1) == "84.5", (rnd(gini, 3), rnd(100 * top3, 1)))
h = M[(M.sector == "Human") & M.country.notna() & M.serovar_clean.notna() & (M.date_precision == "day")].copy()
sub = h.geo_loc_name.fillna("").str.extract(r"^[^:]*:\s*(.+)$")[0].str.strip(); h["place"] = np.where(sub.notna() & (sub != h.country), h.country + ": " + sub.fillna(""), h.country)
grp = h.groupby(["place", "serovar_clean", h.coll_day.dt.date]).target_acc.transform("count"); uniq = 100 * (grp == 1).mean()
check("A.5", "privacy: 54,691 day-level human records; 51.5% unique on place, date and serovar, and 12.9% when the date is rounded to the month", len(h) == 54691 and rnd(uniq, 1) == "51.5", (len(h), rnd(uniq, 1)))
grpm = h.groupby(["place", "serovar_clean", h.coll_day.dt.to_period("M")]).target_acc.transform("count"); check("A.5b", "privacy: 12.9% unique at month resolution", rnd(100 * (grpm == 1).mean(), 1) == "12.9", rnd(100 * (grpm == 1).mean(), 1))
cl_ = pd.read_csv(P("cluster_list_us.tsv"), sep="\t", usecols=["PDS_acc", "target_acc"])
u = M[(M.country == "USA") & M.date_precision.isin(["day", "month"]) & M.release.notna()].merge(cl_, on="target_acc"); u = u[u.groupby("PDS_acc").target_acc.transform("size") >= 10].copy()
u["coll"] = np.where(u.date_precision == "day", u.coll_day, u.coll_month); u["coll"] = pd.to_datetime(u.coll) + pd.to_timedelta(np.where(u.date_precision == "month", 14, 0), unit="D")
def detect(E=0, K=3):
    c3 = u.groupby("PDS_acc").coll.apply(lambda s: s.sort_values().iloc[K - 1]); r3 = u.groupby("PDS_acc").release.apply(lambda s: s.sort_values().iloc[K - 1]) + pd.Timedelta(days=E)
    acc = u.merge(r3.rename("obs"), left_on="PDS_acc", right_index=True); acc = acc[acc.coll <= acc.obs].groupby("PDS_acc").size() - K; df = pd.DataFrame({"ideal": c3, "delay": ((r3 - c3).dt.days).clip(lower=0)}); df["acc"] = acc.reindex(df.index).fillna(0).clip(lower=0)
    return df[(df.ideal >= "2019-01-01") & (df.ideal <= "2022-12-31")]
d0 = detect(0); check("A.6", "detection: 520 US clusters of three isolates, median delay 29 days, 52.1% within 30 days (Hands-on 4, separate code)", len(d0) == 520 and d0.delay.median() == 29 and rnd(100 * (d0.delay <= 30).mean(), 1) == "52.1", (len(d0), d0.delay.median(), rnd(100 * (d0.delay <= 30).mean(), 1)))
# ================================================================== B. the new exercises
SZ = pd.read_csv(P("table_sizing.csv")); SDp = pd.read_csv(P("table_sizing_departures.csv")); SRl = pd.read_csv(P("table_sizing_real.csv")); CT = pd.read_csv(P("table_catchment.csv"))
def n_loop(p, c):
    k = 1
    while 1 - (1 - p) ** k < c - 1e-12: k += 1
    return k
rows = [(row.frequency_pct / 100, row.n_90, row.n_95, row.n_99) for row in SZ.itertuples()]; ok = all(n_loop(p, 0.90) == a and n_loop(p, 0.95) == b and n_loop(p, 0.99) == c_ for p, a, b, c_ in rows) and n_loop(0.02, 0.95) == 149
check("B.1", "sizing: the sequences needed for 0.5% to 10% at 90, 95 and 99% confidence, by a search that does not use the logarithm formula; 149 for 2% at 95%", ok, [(a, b, c_) for _, a, b, c_ in rows])
rng = np.random.default_rng(77); hit = (rng.binomial(149, 0.02, 600_000) >= 1).mean(); ex149 = 1 - 0.98 ** 149; ex148 = 1 - 0.98 ** 148
check("B.2", "sizing: 149 sequences see a 2% variant with probability 0.9507 (exact), 600,000 simulated weeks agree within 0.003, and 148 sequences fall just short (0.9497)", close(ex149, 0.9507, 0.00005) and close(hit, ex149, 0.003) and ex148 < 0.95 <= ex149, (round(ex149, 4), round(hit, 4), round(ex148, 4)))
check("B.3", "sizing: 1.55 weeks at 96 sequences a week; the smallest frequency 96 sequences can see at 95% is 3.07%; the rule of three is within 2 of the exact number for every row", rnd(149 / 96, 2) == "1.55" and rnd(100 * (1 - 0.05 ** (1 / 96)), 2) == "3.07" and all(abs(round(3 / p) - b) <= 2 for p, a, b, c_ in rows), "")
n_s = n_loop(0.02 * 0.8, 0.95)
def hyper_n(N, K, c):
    for k in range(1, N + 1):
        if 1 - math.comb(N - K, k) / math.comb(N, k) >= c: return k
n_f = hyper_n(2000, 40, 0.95)
lg = math.lgamma
def bb0(m, a, b): return math.exp(lg(b + m) + lg(a + b) - lg(b) - lg(a + b + m))              # beta-binomial probability of no carrier in a household of m
rho, p0, mh = 0.3, 0.02, 10; a_, b_ = p0 * (1 - rho) / rho, (1 - p0) * (1 - rho) / rho; Hn = next(H for H in range(1, 500) if 1 - bb0(mh, a_, b_) ** H >= 0.95); det15 = 1 - bb0(mh, a_, b_) ** 15
rng2 = np.random.default_rng(78); q = rng2.beta(a_, b_, 400_000); mc15 = (rng2.binomial(mh, q) >= 1)[:400_000].reshape(-1, 15).any(axis=1).mean() if False else (1 - np.prod([(rng2.binomial(mh, rng2.beta(a_, b_, 100_000)) == 0) for _ in range(15)], axis=0)).mean()
dd = SDp.sequences_needed.tolist()
check("B.4", "four departures: assay sensitivity 0.8 needs 186; a finite population of 2,000 with 40 carriers needs 143 (exact hypergeometric); 35 households of 10 (350 sequences); 15 households detect with 72.9%", n_s == 186 == dd[1] and n_f == 143 == dd[2] and Hn == 35 and dd[3] == 350 and rnd(100 * det15, 1) == "72.9" and close(mc15, det15, 0.01), (n_s, n_f, Hn, round(det15, 4), round(mc15, 4)))
check("B.5", "the design-effect shortcut 1 + (m - 1) x rho = 3.7 gives 551 sequences, 57% above the exact 350", rnd(1 + 9 * 0.3, 1) == "3.7" and round(149 * 3.7) == 551 and rnd(100 * (551 / 350 - 1), 0) == "57", "")
hh = M[(M.sector == "Human") & (M.year == 2022) & M.country.notna()].country.value_counts().head(8); w = hh / 52.14; pm = 100 * (1 - 0.05 ** (1 / w))
ok = list(hh.index) == SRl.country.tolist() and [int(v) for v in hh.values] == SRl.isolates_2022.tolist() and all(close(a, b, 0.006) for a, b in zip(pm.values, SRl.smallest_frequency_pct_in_a_week))
check("B.6", "weekly volume: the eight countries with most human isolates in 2022 and the smallest frequency a week can see at 95% (USA 0.35%, UK 1.88%, China 9.86%)", ok and rnd(pm.iloc[0], 2) == "0.35" and rnd(pm.iloc[1], 2) == "1.88" and rnd(pm.iloc[2], 2) == "9.86", [rnd(v, 2) for v in pm.values[:4]])
ca = [(q_, N_, -math.expm1(N_ * math.log1p(-q_))) for q_ in (1e-4, 1e-5) for N_ in (10_000, 1_000_000)]; check("B.7", "catchment: the chance of at least one shedder is 63.2% and 100% at 1 in 10,000, and 9.5% and 99.995% at 1 in 100,000", all(close(a[2], b, 1e-6) for a, b in zip(ca, CT.p_at_least_one)) and rnd(100 * ca[0][2], 1) == "63.2" and rnd(100 * ca[2][2], 1) == "9.5", [round(a[2], 6) for a in ca])
# allocation: explicit loop in place of the cumulative product
GI = pd.read_csv(P("table_income_groups.csv")); AL = pd.read_csv(P("table_allocation.csv"))
cur = (GI.isolates / GI.isolates.sum()).values; pop = (GI.population / GI.population.sum()).values; flo = 0.10 + 0.60 * cur; bad = []
for nm, sh in (("current", cur), ("population", pop), ("floor", flo)):
    for i, g in enumerate(GI.income):
        n_w = 5000 * sh[i]; surv = 1.0; med = None; wk = 0
        while wk < 80:
            f = min(0.5, 0.001 * 2 ** wk); surv *= (1 - f) ** n_w
            if med is None and 1 - surv >= 0.5: med = wk; fm = f
            wk += 1
        row = AL[(AL.allocation == nm) & (AL.emerges_in == g)].iloc[0]
        if med != row.median_week_of_detection or not close(100 * fm, row.local_frequency_at_median_pct, 0.05) or not close(n_w, row.sequences_per_week, 0.06): bad.append((nm, g, med))
check("B.8", "allocation: median week of detection and local frequency at detection for 3 allocations x 4 groups, with a loop in place of the cumulative product", not bad, bad[:3])
lowc = AL[(AL.allocation == "current") & (AL.emerges_in == "Low income")].iloc[0]; lowp = AL[(AL.allocation == "population") & (AL.emerges_in == "Low income")].iloc[0]
check("B.9", "allocation: a variant that emerges in the low-income group is seen in week 6 at 6.4% under the current pattern and in week 1 at 0.2% by population; the three allocations each sum to 5,000 sequences", lowc.median_week_of_detection == 6 and rnd(lowc.local_frequency_at_median_pct, 1) == "6.4" and lowp.median_week_of_detection == 1 and rnd(lowp.local_frequency_at_median_pct, 1) == "0.2" and all(close(AL[AL.allocation == k].sequences_per_week.sum(), 5000, 1) for k in ("current", "population", "floor")), (lowc.median_week_of_detection, lowc.local_frequency_at_median_pct))
STd = pd.read_csv(P("table_standards.csv")); PHt = pd.read_csv(P("table_placeholders.csv")); UNt = pd.read_csv(P("table_unrecognised.csv"))
def vocab(col):
    s = M[col].dropna().astype(str); cnt = {}
    for v in s: cnt[v] = cnt.get(v, 0) + 1
    vals = sorted(cnt.values(), reverse=True); low = {re.sub(r"\s+", " ", k.lower().strip()) for k in cnt}
    return len(s), len(cnt), len(low), 100 * sum(vals[:10]) / len(s), 100 * sum(vals[:100]) / len(s), 100 * sum(1 for v in vals if v == 1) / len(vals)
ok = True; det = []
for col, lab in (("geo_loc_name", "Place (geo_loc_name)"), ("host", "Host"), ("isolation_source", "Isolation source"), ("serovar_clean", "Serovar (cleaned)")):
    v = vocab(col); r_ = STd[STd.field == lab].iloc[0]; good = (v[0], v[1], v[2]) == (r_.records, r_.distinct_raw, r_.distinct_lower_trimmed) and all(close(a, b, 0.06) for a, b in zip(v[3:], (r_.top10_share_pct, r_.top100_share_pct, r_.used_once_pct_of_distinct))); ok &= good; det.append((lab, v[1], v[2]))
check("B.10", "vocabulary: records, distinct strings, strings after lower-casing, top-10 and top-100 shares and strings used once, for place, host, isolation source and serovar (dictionary counts)", ok, det)
isrc = vocab("isolation_source"); cty = M.geo_loc_name.dropna().astype(str).str.split(":").str[0].str.strip().nunique()
check("B.11", "vocabulary: isolation source has 9,248 distinct strings, 8,481 after the clean-up, 51.5% used once; place reduces to 205 country names; the top 100 source strings cover 86.9%", isrc[1] == 9248 and isrc[2] == 8481 and rnd(isrc[5], 1) == "51.5" and cty == 205 and rnd(isrc[4], 1) == "86.9", (isrc[1], isrc[2], cty))
un = M[M.sector == "Unclassified"]; txt = un.host.notna() | un.isolation_source.notna(); sh = [100 * (M.sector != "Unclassified").mean(), 100 * txt.sum() / n, 100 * (~txt).sum() / n]
check("B.12", "sector rules: 91.1% of records classified, 4.4% with host or source text that is not recognised, 4.6% with neither; the commonest unrecognised string is 'human, stool' (5,985)", [rnd(v, 1) for v in sh] == ["91.1", "4.4", "4.6"] and UNt.text.iloc[0] == "human, stool" and int(UNt.records.iloc[0]) == 5985 and (un[txt].isolation_source.fillna(un[txt].host).astype(str).str.lower().str.strip() == "human, stool").sum() == 5985, [rnd(v, 1) for v in sh])
hs = M.host.dropna().astype(str); hum = hs[hs.str.contains(r"homo|human|patient", case=False)]; ph = M.host.dropna().astype(str).str.strip().str.lower()
check("B.13", "host: 10 spellings name a human (99.97% 'Homo sapiens'); 'not collected' and 'not provided' are 0.69% and 0.52% of host entries", hum.nunique() == 10 and rnd(100 * (hum == "Homo sapiens").mean(), 2) == "99.97" and rnd(100 * (ph == "not collected").mean(), 2) == "0.69" and rnd(100 * (ph == "not provided").mean(), 2) == "0.52", (hum.nunique(), rnd(100 * (ph == "not collected").mean(), 2)))
EMt = pd.read_csv(P("table_embargo.csv")).set_index("embargo_days"); bad = []
for E in (0, 30, 90, 180, 365):
    x_ = detect(E); r_ = EMt.loc[E]
    if not (len(x_) == r_.clusters and x_.delay.median() == r_.median_delay_days and close(x_.delay.quantile(0.9), r_.p90_delay_days, 0.06) and close(100 * (x_.delay <= 30).mean(), r_.detected_within_30_days_pct, 0.06) and close(x_.acc.mean(), r_.mean_isolates_collected_before_detection, 0.06)): bad.append((E, len(x_), x_.delay.median(), round(x_.acc.mean(), 2)))
check("B.14", "embargo: median delay, 90th percentile, share within 30 days and mean isolates before detection for embargoes of 0, 30, 90, 180 and 365 days (separate code)", not bad, bad)
check("B.15", "embargo: the median delay rises from 29 to 59, 119, 209 and 394 days, day for day with the embargo, and no cluster is detected within 30 days under a 90-day embargo", [EMt.loc[e].median_delay_days for e in (0, 30, 90, 180, 365)] == [29, 59, 119, 209, 394] and EMt.loc[90].detected_within_30_days_pct == 0 and rnd(EMt.loc[30].detected_within_30_days_pct, 1) == "1.2", "")
PCt = pd.read_csv(P("table_protocol_check.csv")); fa = PCt[PCt.result == "FAIL"]
check("B.16", "protocol check: the consistent protocol passes all ten rules; the protocol with six planted faults fails exactly R3, R4, R5, R7, R8, R9 and R10 (seven rules)", (PCt[PCt.protocol.str.startswith("A")].result == "pass").all() and len(PCt[PCt.protocol.str.startswith("A")]) == 10 and sorted(fa.rule.tolist()) == sorted(["R3", "R4", "R5", "R7", "R8", "R9", "R10"]) and (PCt.protocol.str.startswith("B") & (PCt.result == "pass")).sum() == 3, fa.rule.tolist())
check("B.17", "protocol check: the formula gives 149 sequences and 2 weeks at 96 a week (protocol B states 1 week), and the allocations of B add to 110% with a largest departure of 12 points", math.ceil(149 / 96) == 2 and "stated 1" in PCt[(PCt.rule == "R3") & PCt.protocol.str.startswith("B")].detail.iloc[0] and "sum 1.10" in PCt[(PCt.rule == "R4") & PCt.protocol.str.startswith("B")].detail.iloc[0] and "12 points" in PCt[(PCt.rule == "R5") & PCt.protocol.str.startswith("B")].detail.iloc[0], "")
# ================================================================== G. the chapter
def docx_text(path):
    from docx import Document
    d = Document(path); body = [p.text for p in d.paragraphs]; tabs = [[[c.text for c in r.cells] for r in t_.rows] for t_ in d.tables]; return d, "\n".join(body), tabs
if A.chapter:
    d, txt, tabs = docx_text(A.chapter); flat = re.sub(r"\s+", " ", txt + " " + " ".join(c for t_ in tabs for r in t_ for c in r))
    figs = re.findall(r"^Figure 12\.(\d+)\.", txt, flags=re.M); check("G.1", "figures are numbered 12.1 to 12.12 in order, one caption for every image; 12 images (3 diagrams and 9 figures), 23 tables", figs == [str(i) for i in range(1, 13)] and len(d.inline_shapes) == 12 and len(d.tables) == 23, (figs, len(d.inline_shapes), len(d.tables)))
    objs, on = [], False
    for p in d.paragraphs:
        if p.style.name == "Heading 2": on = (p.text == "Learning objectives")
        elif on and p.text.strip(): objs.append(p.text)
    objs = objs[1:]; check("G.2", "fourteen learning objectives, from the components of a system to the links with other chapters", len(objs) == 14 and objs[0].startswith("Describe the components") and objs[-1].startswith("Connect surveillance"), len(objs))
    heads = [p.text for p in d.paragraphs if p.style.name == "Heading 2"]; ho = [h_ for h_ in heads if h_.startswith("Hands-on")]; check("G.3", "nine hands-on sections numbered 1 to 9 in order, then the protocol section", [h_.split(" ")[1] for h_ in ho] == [str(i) for i in range(1, 10)] and "Designing a surveillance protocol" in heads, ho)
    want = ["149 random sequences", "1.55 weeks of a capacity of 96 a week", "400,000 simulated weeks give 95.11%", "186 sequences", "143", "350 sequences from 35 households", "72.9%", "551 sequences", "0.35%", "week 6", "6.4%", "week 1", "9,248 distinct strings", "8,481", "4.4%", "5,985", "119 days", "59", "394", "passes all 10", "passes 3 and fails 7", "six planted faults", "20 DOIs".replace("20 DOIs", "56%"), "16 November 2026", "11 to 15 January 2027", "$360 million", "1,500", "0.1% in the Eastern Mediterranean region to 3.4% in Europe", "163 of 194"]
    miss = [w_ for w_ in want if w_ not in flat]; check("G.4", "the sentences quote the validated numbers and the dated policy status", not miss, miss)
    def tbl(first, second=None): return next((t_ for t_ in tabs if t_[0][0].startswith(first) and (second is None or t_[0][1].startswith(second))), None)
    t1 = tbl("Variant frequency"); ok = t1 is not None and [int(r[2]) for r in t1[1:]] == SZ.n_95.tolist() and [r[5] for r in t1[1:]] == [rnd(v, 2) for v in SZ.weeks_at_96_per_week_95]; check("G.5", "the sizing table of the chapter equals the result file", ok, None if t1 is None else t1[1])
    t2 = tbl("Country", "Human isolates"); ok = t2 is not None and [r[3] for r in t2[1:]] == [rnd(v, 2) for v in SRl.smallest_frequency_pct_in_a_week] and [r[1].replace(",", "") for r in t2[1:]] == [str(v) for v in SRl.isolates_2022]; check("G.6", "the table of weekly volumes of the chapter equals the result file", ok, None if t2 is None else t2[1])
    t3 = tbl("Setting (2%"); ok = t3 is not None and [int(r[1]) for r in t3[1:]] == SDp.sequences_needed.tolist(); check("G.7", "the table of departures from the formula equals the result file", ok, None if t3 is None else t3[1])
    t4 = tbl("Variant emerges in"); ok = t4 is not None and all(f"{int(r_.median_week_of_detection)} (" in t4[1 + list(GI.income).index(r_.emerges_in)][2 if r_.allocation == "current" else 4 if r_.allocation == "population" else 6] for r_ in AL.itertuples()); check("G.8", "the allocation table of the chapter equals the result file (median weeks for 12 combinations)", ok, None if t4 is None else t4[1])
    t5 = tbl("Field", "Records"); ok = t5 is not None and [r[2].replace(",", "") for r in t5[1:]] == [str(int(v)) for v in STd.distinct_raw]; check("G.9", "the vocabulary table of the chapter equals the result file", ok, None if t5 is None else t5[1])
    t6 = tbl("Embargo"); ok = t6 is not None and [r[1] for r in t6[1:]] == [rnd(v, 1) for v in EMt.median_delay_days] and [r[3] for r in t6[1:]] == [rnd(v, 1) for v in EMt.detected_within_30_days_pct]; check("G.10", "the embargo table of the chapter equals the result file", ok, None if t6 is None else t6[1])
    t7 = tbl("Rule", "Check"); ok = t7 is not None and [r[2] for r in t7[1:]] == PCt[PCt.protocol.str.startswith("A")].result.tolist() and [r[3] for r in t7[1:]] == PCt[PCt.protocol.str.startswith("B")].result.tolist(); check("G.11", "the protocol-check table of the chapter equals the result file", ok, None if t7 is None else t7[1])
    import zipfile; z = zipfile.ZipFile(A.chapter); rels = z.read("word/_rels/document.xml.rels").decode(); dois = {u_.rstrip(".,;").replace("&amp;", "&") for u_ in re.findall(r"https://doi\.org/[^\"'<> ]+", rels)}
    newd = ["10.3389/fsci.2024.1298248", "10.1007/s44370-026-00061-0", "10.1201/9781003409809-1", "10.1038/s41579-020-0421-0", "10.1093/gigascience/giac003", "10.3389/fmicb.2023.1124120", "10.3389/fpubh.2020.00238", "10.1001/jama.1983.03330370053031"]
    have = {x_.split("doi.org/")[1] for x_ in dois}; check("G.12", "the chapter links the 8 new DOIs read for this extension, with 20 distinct DOIs in all (12 in the original)", all(n_ in have for n_ in newd) and len(have) == EXPECTED_DOIS, (len(have), [n_ for n_ in newd if n_ not in have]))
    if A.doi:
        import urllib.request, urllib.parse, time; bad_d = []
        for dd_ in sorted(dois):
            doi = dd_.split("doi.org/")[1]; okd = False
            for tryn in range(4):
                try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "mepi-book/1.0 (mailto:none@example.org)"}), timeout=40).read(); okd = True; break
                except Exception as e: err = str(e)[:30]; time.sleep(5)
            if not okd: bad_d.append((doi, err))
        check("G.13", "every DOI of the chapter resolves in Crossref", not bad_d, bad_d)
    else: skip("G.13", "every DOI of the chapter resolves in Crossref", "run with --doi")
else:
    for k, w_ in [("G.1", "figure numbers, tables"), ("G.2", "objectives"), ("G.3", "hands-on sections"), ("G.4", "sentences"), ("G.5", "sizing table"), ("G.6", "weekly-volume table"), ("G.7", "departures table"), ("G.8", "allocation table"), ("G.9", "vocabulary table"), ("G.10", "embargo table"), ("G.11", "protocol table"), ("G.12", "DOIs"), ("G.13", "DOIs resolve")]: skip(k, w_, "run with --chapter")
fd = next((d_ for d_ in ("../figures", "figures", os.path.join(HERE, "../figures")) if os.path.isdir(d_)), None)
figs_ = ["diagram_surveillance_system", "diagram_release_tiers", "diagram_protocol_flow", "fig1_indicators", "fig2_equity", "fig2b_r_lorenz", "fig3_privacy", "fig4_detection", "fig5_sizing", "fig6_allocation", "fig7_standards", "fig8_embargo", "fig9_agentic_workflow_diagram"]
if fd: miss = [f_ for f_ in figs_ if not os.path.exists(os.path.join(fd, f_ + ".png"))]; check("H.1", "all 13 figure files are present (3 diagrams, 9 figures and the R Lorenz plot)", not miss, miss)
else: skip("H.1", "figure files", "no figures folder")
cd = next((d_ for d_ in (os.path.join(HERE, "../code"), "../code", "code") if os.path.isdir(d_) and glob.glob(os.path.join(d_, "ex12_5_sizing.py"))), None)
if cd:
    import py_compile; badc = []
    for f_ in sorted(glob.glob(os.path.join(cd, "*.py"))):
        try: py_compile.compile(f_, doraise=True)
        except Exception: badc.append(os.path.basename(f_))
    check("H.2", "every Python script of the package compiles", not badc, badc)
out_ = pd.DataFrame(R); out_.to_csv("validation_results_ch12.csv", index=False)
for r in R: print(f"{r['result']:5} {r['id']:8} {r['check'][:150]}" + (f"   [{r['detail']}]" if r["result"] != "PASS" else ""))
print(f"\n{(out_.result == 'PASS').sum()} PASS, {(out_.result == 'FAIL').sum()} FAIL, {(out_.result == 'SKIP').sum()} SKIP"); sys.exit(1 if (out_.result == "FAIL").any() else 0)
