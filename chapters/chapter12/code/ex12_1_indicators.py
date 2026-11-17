"""
Chapter 12, Hands-on 1 - Surveillance indicators from public metadata: completeness and timeliness.
Input: sal_surv.pkl from ex12_0_prepare.py (NCBI Pathogen Detection Salmonella, 911,210 isolates).
Output: table_completeness.csv, table_timeliness.csv, Figure 1.
"""
import numpy as np, pandas as pd
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
m = load_surv(); order = ["Human", "Food", "Animal", "Environment"]; COL = {"Human": "#2b6cb0", "Food": "#dd6b20", "Animal": "#2f855a", "Environment": "#805ad5"}
print(f"{len(m):,} isolates")

# ---------------------------------------------------------------- completeness: share of records carrying each field
sub = m.geo_loc_name.fillna("").str.split(":", n=1).str[1].str.strip().replace("", np.nan); sub_ok = sub.notna() & (sub != m.country)
fields = {"Any collection date": m.date_precision != "missing", "Collection date to the day": m.date_precision == "day", "Country": m.country.notna(), "Place below country level": sub_ok,
          "Host": m.host.notna(), "Isolation source": m.isolation_source.notna(), "Sector can be assigned": m.sector != "Unclassified", "Coordinates": m.has_coords}
C = pd.DataFrame({k: v.groupby(m.sector).mean().reindex(order) for k, v in fields.items()}).T * 100; C["All isolates"] = [v.mean() * 100 for v in fields.values()]; C.round(1).to_csv("table_completeness.csv"); print(C.round(1).to_string())
Y = m[m.year.between(2012, 2024)].groupby("year").apply(lambda d: pd.Series({"any_date": (d.date_precision != "missing").mean() * 100, "day_level": (d.date_precision == "day").mean() * 100, "source": d.isolation_source.notna().mean() * 100}))

# ---------------------------------------------------------------- timeliness
d = m[m.date_precision.isin(["day", "month"]) & m.release.notna() & m.coll_month.notna()].copy()
d["months"] = (d.release.dt.year - d.coll_month.dt.year) * 12 + (d.release.dt.month - d.coll_month.dt.month)
print(f"records with a collection month and a release date: {len(d):,}; release before the collection month: {(d.months < 0).sum():,}"); d = d[d.months >= 0]
rows = []
for s in order + ["All isolates"]:
    x = d if s == "All isolates" else d[d.sector == s]; x = x[x.year.between(2015, 2023)]
    dd = x[x.date_precision == "day"]; dd = dd[dd.lag_days >= 0]
    rows.append(dict(sector=s, n_month_or_day=len(x), within_1_month=(x.months <= 1).mean() * 100, within_3_months=(x.months <= 3).mean() * 100, within_12_months=(x.months <= 12).mean() * 100, median_months=x.months.median(),
                     n_day_level=len(dd), median_days_day_level=dd.lag_days.median(), within_21_days_day_level=(dd.lag_days <= 21).mean() * 100 if len(dd) else np.nan))
T = pd.DataFrame(rows).set_index("sector"); T.round(1).to_csv("table_timeliness.csv"); print(T.round(1).to_string())
print("day-level share of isolates by sector:", (m.groupby("sector").apply(lambda g: (g.date_precision == "day").mean()).reindex(order) * 100).round(1).to_dict())

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.7), gridspec_kw={"width_ratios": [1.25, 1, 1]})
show = ["Any collection date", "Collection date to the day", "Place below country level", "Host", "Isolation source", "Coordinates"]; w = 0.2
for i, s in enumerate(order): ax[0].barh(np.arange(len(show)) + (i - 1.5) * w, C.loc[show, s].values, height=w, color=COL[s], label=s)
ax[0].set_yticks(range(len(show))); ax[0].set_yticklabels(show, fontsize=9); ax[0].invert_yaxis(); ax[0].set(xlabel="Records with the field (%)", title="Completeness by sector"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4)
for s in order:
    x = d[(d.sector == s) & d.year.between(2015, 2023)].months; xs = np.arange(0, 37); ax[1].plot(xs, [(x <= k).mean() * 100 for k in xs], color=COL[s], label=f"{s} (n={len(x):,})")
ax[1].axvline(1, color="#c53030", ls="--", lw=1); ax[1].text(1.6, 8, "1 month", color="#c53030", fontsize=8.5); ax[1].set(xlim=(0, 36), xlabel="Months from collection month to public release", ylabel="Cumulative % of isolates", title="Timeliness (collected 2015-2023)"); ax[1].legend(frameon=False, fontsize=8.2, loc="lower right")
for s in order:
    t = d[(d.sector == s)].groupby("year").months.median(); t = t[(t.index >= 2012) & (t.index <= 2023)]; n = d[d.sector == s].groupby("year").size().reindex(t.index); ax[2].plot(t.index[n >= 100], t[n >= 100], marker="o", ms=3.5, color=COL[s], label=s)
ax[2].set(xlabel="Collection year", ylabel="Median months from collection month to release", title="Lag by collection year", ylim=(0, None)); fig.tight_layout(); fig.savefig("figs12/fig1_indicators.png", dpi=200); plt.close(fig); print("figure written")
