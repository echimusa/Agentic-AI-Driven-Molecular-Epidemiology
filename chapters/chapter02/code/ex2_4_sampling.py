"""
Chapter 2, Hands-on 4 - Sampling for molecular studies: sample size for an association study, and how much of the world's pathogen sequencing is openly sampled.
Part A checks the sample-size formula of the draft chapter, N = (Z_a + Z_b)^2 / (2 p (1 - p) (ln OR)^2), against a simulation of a balanced case-control study of an additive genetic effect, and against the same formula
with the outcome variance pi (1 - pi) added. Part B uses public SARS-CoV-2 data (weekly cases from Our World in Data / WHO; sequences by sampling date from the open CoV-Spectrum API, INSDC data) for 36 countries in 2021-2022
to compute the openly shared sequencing fraction, the smallest variant prevalence a week's sequences could detect with 95% probability, and the share of weeks that could detect a 1% variant.
Output: table_sample_size.csv, table_sequencing_coverage.csv, Figure 4.
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(2)
# ---------------------------------------------------------------- Part A
def n_formula(p, OR, alpha, power=0.8, pi=None):
    z = stats.norm.isf(alpha / 2) + stats.norm.isf(1 - power); base = z ** 2 / (2 * p * (1 - p) * np.log(OR) ** 2); return base if pi is None else base / (pi * (1 - pi))
def sim_power(N, p, OR, alpha, reps=3000):
    f = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2]); fc = f * OR ** np.arange(3); fc = fc / fc.sum(); n1 = int(round(N / 2)); crit = stats.chi2.isf(alpha, 1); hits = 0
    ca = rng.multinomial(n1, fc, reps); co = rng.multinomial(n1, f, reps); s = np.arange(3); tot = ca + co; n = tot.sum(axis=1)
    # Cochran-Armitage trend test
    t = (ca * s).sum(axis=1) - n1 * (tot * s).sum(axis=1) / n                                   # observed minus expected score sum among cases
    var = (n1 * n1 / n) * (((tot * s * s).sum(axis=1) / n) - ((tot * s).sum(axis=1) / n) ** 2); stat = t ** 2 / var; return (stat > crit).mean()
rows = []
for p in (0.05, 0.1, 0.3, 0.5):
    for OR in (1.2, 1.5):
        for alpha in (0.05, 5e-8):
            n1 = n_formula(p, OR, alpha); n2 = n_formula(p, OR, alpha, pi=0.5); rows.append(dict(maf=p, OR=OR, alpha=alpha, N_draft_formula=n1, N_with_outcome_variance=n2, power_at_N_draft=sim_power(n1, p, OR, alpha), power_at_N_with_variance=sim_power(n2, p, OR, alpha)))
A = pd.DataFrame(rows); A.round(3).to_csv("table_sample_size.csv", index=False); print(A.round(3).to_string(index=False))
# ---------------------------------------------------------------- Part B
seq = pd.read_csv("lapis_daily_2021_2022.csv", parse_dates=["date"]); cases = pd.read_csv("owid_cases_2021_2022.csv", parse_dates=["date"])
weeks = pd.date_range("2021-01-03", "2022-12-25", freq="W-SUN"); out = []
for (lc, oc), g in seq.groupby(["lapis", "owid"]):
    sw = g.set_index("date").n.resample("W-SUN").sum().reindex(weeks, fill_value=0); cw = cases[cases.location == oc].set_index("date").new_cases.resample("W-SUN").sum().reindex(weeks, fill_value=0)
    for yr in (2021, 2022):
        m = weeks.year == yr; n = sw[m]; c = cw[m]; pmin = 1 - 0.05 ** (1 / n[n > 0]); out.append(dict(country=lc, year=yr, cases=c.sum(), sequences=n.sum(), fraction_pct=100 * n.sum() / max(c.sum(), 1), weeks_with_any=(n > 0).mean() * 100, median_weekly=n.median(), weeks_detect_1pct=(n >= 299).mean() * 100, median_min_detectable_pct=100 * pmin.median() if len(pmin) else np.nan))
C = pd.DataFrame(out); C.round(4).to_csv("table_sequencing_coverage.csv", index=False)
tot = C.groupby("year")[["cases", "sequences"]].sum(); print("all 36 countries:", {y: f"{r.sequences:,.0f} open sequences for {r.cases:,.0f} cases ({100 * r.sequences / r.cases:.2f}%)" for y, r in tot.iterrows()})
w = C.pivot(index="country", columns="year", values="fraction_pct"); print("countries above the 0.5% benchmark: 2021", int((w[2021] >= 0.5).sum()), "| 2022", int((w[2022] >= 0.5).sum()), "of", len(w)); print(C[C.year == 2022].sort_values("sequences", ascending=False).head(8)[["country", "cases", "sequences", "fraction_pct", "weeks_detect_1pct", "median_min_detectable_pct"]].round(3).to_string(index=False))
c22 = C[C.year == 2022]; print("2022: countries with a week able to detect a 1% variant at least once:", int((c22.weeks_detect_1pct > 0).sum()), "of", len(c22), "| countries where more than half of weeks could:", int((c22.weeks_detect_1pct > 50).sum()))
print("share of all open sequences from the four largest contributors (USA, UK, Germany, Denmark):", round(seq[seq.lapis.isin(["USA", "United Kingdom", "Germany", "Denmark"])].n.sum() / seq.n.sum(), 3))
fig, ax = plt.subplots(1, 3, figsize=(13.8, 5.0), gridspec_kw={"width_ratios": [1.15, 1, 1]})
o = C[C.year == 2022].sort_values("fraction_pct"); ax[0].barh(range(len(o)), np.maximum(o.fraction_pct, 1e-4), color=np.where(o.fraction_pct >= 0.5, "#2b6cb0", "#a0aec0")); ax[0].set_xscale("log"); ax[0].set_yticks(range(len(o))); ax[0].set_yticklabels(o.country, fontsize=7); ax[0].axvline(0.5, color="#c53030", ls="--", lw=1); ax[0].text(0.55, 1, "0.5% benchmark", color="#c53030", fontsize=8); ax[0].text(1.2e-4, 0.2, "no open sequences", fontsize=7, color="#4a5568"); ax[0].set(xlabel="Open sequences per 100 reported cases, 2022 (log scale)", title="Who is openly sampled")
d = C[C.year == 2022].sort_values("weeks_detect_1pct"); ax[1].barh(range(len(d)), d.weeks_detect_1pct, color="#2b6cb0"); ax[1].set_yticks(range(len(d))); ax[1].set_yticklabels(d.country, fontsize=7); ax[1].set(xlabel="Weeks in 2022 with at least 299 open sequences (%)", title="Weeks able to detect a 1% variant")
nn = np.logspace(0, 4.5, 200)
for pv, c in [(0.001, "#805ad5"), (0.005, "#2f855a"), (0.01, "#2b6cb0"), (0.05, "#dd6b20")]: ax[2].plot(nn, 1 - (1 - pv) ** nn, color=c, label=f"prevalence {pv * 100:g}%")
ax[2].axhline(0.95, color="black", ls=":", lw=1); ax[2].set_xscale("log"); ax[2].set(xlabel="Sequences sampled at random in a week", ylabel="Probability of detecting the variant", title="Detection needs sequences, not only cases"); ax[2].legend(frameon=False, fontsize=8.3, loc="lower right")
for lc in ["United Kingdom", "South Africa", "India"]:
    mw = C[(C.country == lc) & (C.year == 2022)].median_weekly.iloc[0]
    if mw > 0: ax[2].scatter([mw], [1 - 0.99 ** mw], color="black", s=18, zorder=3); ax[2].annotate(lc, (mw, 1 - 0.99 ** mw), xytext=(4, -12), textcoords="offset points", fontsize=7.5)
fig.tight_layout(); fig.savefig("figs2/fig4_sampling.png", dpi=200); plt.close(fig); print("figure written")
