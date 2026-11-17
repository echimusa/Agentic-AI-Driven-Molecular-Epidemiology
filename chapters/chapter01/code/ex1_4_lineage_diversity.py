"""
Chapter 1, Hands-on 4 - An application: lineage diversity and variant sweeps in SARS-CoV-2 from open sequence counts.
Data: monthly counts of sequences by Pango lineage and Nextstrain clade for the United Kingdom and the USA, January 2020 to December 2024, from the open CoV-Spectrum LAPIS API (INSDC data).
Simpson's index of diversity is computed on lineages each month, as observed and after rarefying to 500 sequences so that months with different sequencing effort can be compared.
Output: table_lineage_diversity.csv and Figure 7 of the chapter (figs1/fig4_lineage_diversity.png).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import hypergeom
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
d = pd.read_csv("lapis_lineage_monthly.csv", keep_default_na=False); d["lineage"] = d.lineage.replace("", "unassigned"); rng = np.random.default_rng(4); N0 = 500
def simpson(counts):
    counts = np.asarray(counts, float); N = counts.sum(); return 1 - (counts * (counts - 1)).sum() / (N * (N - 1)) if N > 1 else np.nan
rows = []
for (c, m), g in d.groupby(["country", "month"]):
    cnt = g.groupby("lineage")["count"].sum().values; N = int(cnt.sum()); rare = np.nan; rrich = np.nan                      # one type per lineage name: the same lineage under two clades is one type
    if N >= N0:
        subs = [rng.multivariate_hypergeometric(cnt.astype(int), N0) for _ in range(100)]; rare = np.mean([simpson(x) for x in subs]); rrich = float((1 - hypergeom(N, cnt.astype(int), N0).pmf(0)).sum())                               # exact expected richness of a draw of N0 sequences without replacement
    top = g.groupby("lineage")["count"].sum().sort_values(ascending=False); rows.append(dict(country=c, month=m, sequences=N, lineages=g.lineage.nunique(), top_lineage=top.index[0], top_share=top.iloc[0] / N, D_observed=simpson(cnt), D_rarefied=rare, lineages_rarefied=rrich))
T = pd.DataFrame(rows); T["month_dt"] = pd.PeriodIndex(T.month, freq="M").to_timestamp(); T.drop(columns="month_dt").round(4).to_csv("table_lineage_diversity.csv", index=False)
uk = T[T.country == "United Kingdom"].set_index("month_dt"); us = T[T.country == "USA"].set_index("month_dt")
print(f"UK sequences {int(uk.sequences.sum()):,}; USA {int(us.sequences.sum()):,}; months with at least {N0} sequences: UK {int((uk.sequences >= N0).sum())} of {len(uk)}, USA {int((us.sequences >= N0).sum())} of {len(us)}")
for lab, t in [("UK", uk), ("USA", us)]:
    ok = t[t.D_rarefied.notna()]; print(f"{lab}: largest gap between observed and rarefied D in any month: {(ok.D_observed - ok.D_rarefied).abs().max():.3f}; richness (lineages seen) vs rarefied richness: correlation with sequences {stats.spearmanr(t.sequences, t.lineages).correlation:.2f} against {stats.spearmanr(ok.sequences, ok.lineages_rarefied).correlation:.2f}")
    print(f"{lab}: rarefied D ranges {ok.D_rarefied.min():.3f} (month {ok.D_rarefied.idxmin():%Y-%m}, top lineage {ok.loc[ok.D_rarefied.idxmin(), 'top_lineage']}, {ok.loc[ok.D_rarefied.idxmin(), 'top_share']:.0%} of sequences) to {ok.D_rarefied.max():.3f} (month {ok.D_rarefied.idxmax():%Y-%m})")
    rho = stats.spearmanr(t.sequences, t.D_observed, nan_policy="omit"); print(f"{lab}: Spearman correlation of observed D with sequences per month: {rho.correlation:.2f}")
for m in ["2020-06", "2020-12", "2021-04", "2021-07", "2021-12", "2022-03", "2022-12", "2023-06", "2024-06"]:
    r = uk.loc[pd.Timestamp(m + "-01")]; print(f"UK {m}: {int(r.sequences):>6,} sequences; {int(r.lineages):>4} lineages; top {r.top_lineage} {r.top_share:.0%}; D observed {r.D_observed:.3f}; D rarefied {r.D_rarefied if r.D_rarefied == r.D_rarefied else float('nan'):.3f}")
def group(c):
    if c.startswith("19") or c in ("20A", "20B", "20C", "20D", "20E", "20F", "20G"): return "Early lineages"
    return {"20I": "Alpha", "20H": "Beta", "20J": "Gamma", "21A": "Delta", "21I": "Delta", "21J": "Delta", "21K": "Omicron BA.1", "21L": "Omicron BA.2"}.get(c, "Later Omicron and recombinants" if (c[:2] in ("22", "23", "24") or c == "recombinant") else "Other")
d["group"] = d.clade.map(group); sh = d[d.country == "United Kingdom"].groupby(["month", "group"])["count"].sum().unstack(fill_value=0); sh = sh.div(sh.sum(axis=1), axis=0); sh.index = pd.PeriodIndex(sh.index, freq="M").to_timestamp()
order = ["Early lineages", "Alpha", "Beta", "Gamma", "Delta", "Omicron BA.1", "Omicron BA.2", "Later Omicron and recombinants", "Other"]; sh = sh[[c for c in order if c in sh.columns]]; colors = ["#a0aec0", "#dd6b20", "#d69e2e", "#805ad5", "#2f855a", "#c53030", "#e53e3e", "#2b6cb0", "#cbd5e0"]
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [1.15, 1, 0.9]}); ax[0].stackplot(sh.index, [sh[c] * 100 for c in sh.columns], labels=sh.columns, colors=colors[: len(sh.columns)]); ax[0].set(ylabel="Share of UK sequences (%)", title="Variant sweeps, United Kingdom", ylim=(0, 100)); ax[0].legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3); ax[0].tick_params(axis="x", rotation=30)
for t, lab, c in [(uk, "UK", "#2b6cb0"), (us, "USA", "#dd6b20")]: ok = t[t.D_rarefied.notna()]; ax[1].plot(ok.index, ok.D_rarefied, "o-", ms=3, color=c, label=lab)
ax[1].set(ylabel="Simpson's index at lineage level (500 sequences)", title="Diversity falls during a sweep", ylim=(0, 1.02)); ax[1].legend(frameon=False, fontsize=8.5); ax[1].tick_params(axis="x", rotation=30)
ax[2].scatter(uk.sequences, uk.lineages, s=14, color="#c53030", label="Lineages seen", alpha=0.8); ok = uk[uk.D_rarefied.notna()]; ax[2].scatter(ok.sequences, ok.lineages_rarefied, s=14, color="#2b6cb0", label="Expected in 500 sequences", alpha=0.8); ax[2].set_xscale("log"); ax[2].axvline(N0, color="black", ls="--", lw=0.8)
ax[2].set(xlabel="Sequences in the month (log scale)", ylabel="Number of lineages, UK", title="Richness depends on sampling effort"); ax[2].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig("figs1/fig4_lineage_diversity.png", dpi=200); plt.close(fig); print("figure written")
