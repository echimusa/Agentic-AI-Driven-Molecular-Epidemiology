"""
Chapter 13, Hands-on 4 - Variant replacement and a two-variant renewal forecast from public data.
Data: England daily cases by specimen date (UKHSA data dashboard API) and UK SARS-CoV-2 sequences by Nextstrain clade and sampling date (CoV-Spectrum open LAPIS API, INSDC data).
Models (both renewal equations, generation interval gamma mean 5 d, SD 2.5 d, R fixed over the 14-day horizon):
  A  cases only: one reproduction number estimated from the last 10 days of total cases.
  B  cases plus genomes: total cases are split into Delta-like and Omicron by the logistic fit to the sequence shares (using sequences up to 7 days before the origin),
     each variant gets its own reproduction number, and the forecasts are summed.
Scoring: weighted interval score (Bracher et al. 2021) and log error of the 7-day mean of cases 8-14 days ahead, for origins from 29 Nov 2021 to 5 Jan 2022.
Output: table_variant_forecasts.csv, table_variant_summary.csv, Figure 4.
"""
import json, os, urllib.request, urllib.parse
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from ebolacore import gamma_si
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(21); W = gamma_si(5.0, 2.5, 25); H, DRAWS, NOISE = 14, 2000, 0.10
if not os.path.exists("lapis_uk_clades.csv"):
    u = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated?" + urllib.parse.urlencode({"country": "United Kingdom", "dateFrom": "2021-10-01", "dateTo": "2022-02-28", "fields": "date,nextstrainClade"})
    pd.DataFrame(json.load(urllib.request.urlopen(u, timeout=120))["data"]).to_csv("lapis_uk_clades.csv", index=False)
if not os.path.exists("england_cases.csv"):
    u = "https://api.ukhsa-dashboard.data.gov.uk/themes/infectious_disease/sub_themes/respiratory/topics/COVID-19/geography_types/Nation/geographies/England/metrics/COVID-19_cases_casesByDay?page_size=365"; rows = []
    while u: d = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "ch13"}), timeout=60)); rows += d["results"]; u = d["next"]
    pd.DataFrame(rows)[["date", "metric_value"]].assign(date=lambda x: pd.to_datetime(x.date)).sort_values("date").to_csv("england_cases.csv", index=False)
cases = pd.read_csv("england_cases.csv", parse_dates=["date"]).set_index("date").metric_value.asfreq("D"); cases = cases["2021-10-01":"2022-02-15"]; sm = cases.rolling(7).mean()
g = pd.read_csv("lapis_uk_clades.csv", parse_dates=["date"]); g["om"] = g.nextstrainClade.isin(["21K", "21L", "21M"])
seq = g.groupby(["date", "om"])["count"].sum().unstack(fill_value=0).rename(columns={True: "omicron", False: "other"}).reindex(cases.index, fill_value=0)
print(f"{int(seq.sum().sum()):,} UK sequences, {int(cases.sum()):,} England cases; Omicron share of sequences: {seq.omicron.sum() / seq.sum().sum():.2f}")

def logistic_share(t0, window=28, lag=7):
    """Omicron share by day from a binomial logistic fit to the sequences dated in [t0-lag-window, t0-lag], extrapolated to t0 and beyond."""
    s = seq.loc[t0 - pd.Timedelta(days=lag + window): t0 - pd.Timedelta(days=lag)]; x = np.arange(len(s), dtype=float); n = (s.omicron + s.other).values + 1.0; k = s.omicron.values + 0.5
    b = np.array([-3.0, 0.1])
    for _ in range(60):
        eta = np.clip(b[0] + b[1] * x, -30, 30); p = 1 / (1 + np.exp(-eta)); Wt = n * p * (1 - p) + 1e-9; X = np.column_stack([np.ones_like(x), x]); z = eta + (k - n * p) / Wt
        nb = np.linalg.solve((X.T * Wt) @ X + 1e-4 * np.eye(2), (X.T * Wt) @ z)
        if np.max(np.abs(nb - b)) < 1e-8: b = nb; break
        b = nb
    return lambda days_after_start: 1 / (1 + np.exp(-np.clip(b[0] + b[1] * days_after_start, -30, 30))), len(s), b[1]

def cori_post(series, window=10, a=1.0, bb=5.0):
    T = len(series); lam = np.array([np.dot(W[: min(t, len(W))], series[t - 1::-1][: min(t, len(W))]) if t > 0 else 0.0 for t in range(T)]); I, L = series[-window:].sum(), lam[-window:].sum(); return a + I, 1 / bb + L
def project(hist, Rdraws):
    T = len(hist); out = np.zeros((len(Rdraws), T + H)); out[:, :T] = hist
    for h in range(H):
        t = T + h; s = min(t, len(W)); out[:, t] = Rdraws * (out[:, t - s:t][:, ::-1] * W[:s]).sum(axis=1)
    return out[:, T:]
def wis(y, q):
    qs = {a: np.quantile(q, [a / 2, 1 - a / 2]) for a in (0.05, 0.5)}; m = np.median(q); K = 2; s = 0.5 * abs(y - m)
    for a, (lo, hi) in qs.items(): s += (a / 2) * ((hi - lo) + (2 / a) * max(lo - y, 0) + (2 / a) * max(y - hi, 0))
    return s / (K + 0.5)

rows = []; origins = pd.date_range("2021-11-29", "2022-01-05", freq="2D")
for t0 in origins:
    hist = sm.loc[:t0].dropna().values; idx = sm.loc[:t0].dropna().index; shape, rate = cori_post(hist); RA = stats.gamma.rvs(a=shape, scale=1 / rate, size=DRAWS, random_state=rng)
    FA = project(hist, RA)
    f, nseq, slope = logistic_share(t0); start = t0 - pd.Timedelta(days=7 + 28); share = np.array([f((d - start).days) for d in idx]); histO, histD = hist * share, hist * (1 - share)
    sD, rD = cori_post(histD + 1e-9); sO, rO = cori_post(histO + 1e-9); RD = stats.gamma.rvs(a=sD, scale=1 / rD, size=DRAWS, random_state=rng); RO = stats.gamma.rvs(a=sO, scale=1 / rO, size=DRAWS, random_state=rng)
    FB = project(histD, RD) + project(histO, RO)
    y = cases.loc[t0 + pd.Timedelta(days=8): t0 + pd.Timedelta(days=14)].mean()
    rec = dict(origin=t0, observed=y, omicron_share_at_origin=float(share[-1]), growth_adv_per_day=slope, R_total=shape / rate, R_delta=sD / rD, R_omicron=sO / rO)
    for nm, F in [("A", FA), ("B", FB)]:
        pred = F[:, 7:14].mean(axis=1) * np.exp(rng.normal(0, NOISE, DRAWS)); q = np.quantile(pred, [0.025, 0.5, 0.975]); rec.update({f"{nm}_median": q[1], f"{nm}_lo": q[0], f"{nm}_hi": q[2], f"{nm}_wis": wis(y, pred), f"{nm}_logerr": np.log(q[1] / y), f"{nm}_in95": bool(q[0] <= y <= q[2])})
    rows.append(rec)
R = pd.DataFrame(rows); R.round(4).to_csv("table_variant_forecasts.csv", index=False)
S = pd.DataFrame({m: dict(mean_WIS=R[f"{m}_wis"].mean(), median_WIS=R[f"{m}_wis"].median(), mean_abs_log_error=R[f"{m}_logerr"].abs().mean(), mean_log_error=R[f"{m}_logerr"].mean(), coverage_95=R[f"{m}_in95"].mean() * 100) for m in "AB"}).T
S.index = ["A: cases only", "B: cases plus genomes"]; S.round(3).to_csv("table_variant_summary.csv"); print(S.round(3).to_string()); print(f"origins: {len(R)}; B better (lower WIS) at {(R.B_wis < R.A_wis).sum()} of {len(R)}")
ph = R[(R.origin >= "2021-12-05") & (R.origin <= "2021-12-27")]; print("during the replacement (5 to 27 Dec):", S.index[0], round(ph.A_wis.mean(), 0), round(ph.A_logerr.mean(), 3), "|", S.index[1], round(ph.B_wis.mean(), 0), round(ph.B_logerr.mean(), 3))
R["phase"] = pd.cut(R.omicron_share_at_origin, [-1, 0.3, 0.9, 1.01], labels=["Omicron share below 30%", "30% to 90%", "above 90%"])
P = R.groupby("phase", observed=True).agg(origins=("origin", "size"), A_mean_WIS=("A_wis", "mean"), B_mean_WIS=("B_wis", "mean"), A_mean_logerr=("A_logerr", "mean"), B_mean_logerr=("B_logerr", "mean")); P.round(3).to_csv("table_variant_phases.csv"); print(P.round(3).to_string())
print(R[["origin", "omicron_share_at_origin", "R_total", "R_delta", "R_omicron", "observed", "A_median", "B_median"]].round(2).iloc[::3].to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6), gridspec_kw={"width_ratios": [1.1, 1.1, 1]})
sh = seq.omicron.rolling(7, center=True).sum() / (seq.omicron + seq.other).rolling(7, center=True).sum(); a2 = ax[0].twinx(); ax[0].bar(cases.index, cases.values, width=1, color="#cbd5e0"); ax[0].plot(sm.index, sm.values, color="#2b6cb0", lw=2)
a2.plot(sh.index, sh * 100, color="#c53030", lw=2); a2.set_ylabel("Omicron share of sequences (%, red)"); a2.spines["top"].set_visible(False); ax[0].set(ylabel="England cases per day (blue: 7-day mean)", title="Cases and the genomic share"); ax[0].tick_params(axis="x", rotation=30)
ax[1].plot(cases.index, sm, color="black", lw=1.6, label="Observed (7-day mean)"); 
for nm, c in [("A", "#dd6b20"), ("B", "#2f855a")]:
    ax[1].plot(R.origin + pd.Timedelta(days=11), R[f"{nm}_median"], "o-", ms=3.5, color=c, label=("A: cases only" if nm == "A" else "B: cases + genomes")); ax[1].fill_between(R.origin + pd.Timedelta(days=11), R[f"{nm}_lo"], R[f"{nm}_hi"], color=c, alpha=0.12)
ax[1].set(xlim=(pd.Timestamp("2021-12-05"), pd.Timestamp("2022-01-25")), ylabel="Cases per day", title="Forecasts 8 to 14 days ahead, by target date"); ax[1].legend(frameon=False, fontsize=8.3, loc="upper left"); ax[1].tick_params(axis="x", rotation=30)
ax[2].axhline(0, color="black", lw=0.8); ax[2].plot(R.origin, R.A_logerr, "o-", ms=3.5, color="#dd6b20", label="A: cases only"); ax[2].plot(R.origin, R.B_logerr, "o-", ms=3.5, color="#2f855a", label="B: cases + genomes"); ax[2].set(ylabel="Log error of the forecast median", title="Error by forecast origin"); ax[2].legend(frameon=False, fontsize=8.5); ax[2].tick_params(axis="x", rotation=30)
for a_ in (ax[0], ax[1], ax[2]): a_.xaxis.set_major_locator(mdates.DayLocator(interval=14 if a_ is ax[0] else 10)); a_.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
fig.tight_layout(); fig.savefig("figs13/fig4_variant_forecast.png", dpi=200); plt.close(fig); print("figure written")
