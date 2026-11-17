"""
Chapter 12, Hands-on 3 - Privacy risk and utility: k-anonymity of place and date for human isolates.
A record is unique (k = 1) when no other record shares its quasi-identifiers. Quasi-identifiers are place, collection date and serovar, all public in the same table.
Output: table_privacy.csv and Figure 3.
"""
import numpy as np, pandas as pd
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
m = load_surv(); h = m[(m.sector == "Human") & m.country.notna() & m.serovar_clean.notna()].copy()
print(f"{len(h):,} human isolates with a country; collection date precision: {h.date_precision.value_counts().to_dict()}; coordinates present: {h.has_coords.mean():.1%}")
sub = h.geo_loc_name.fillna("").str.split(":", n=1).str[1].str.strip().replace("", np.nan); h["place_sub"] = np.where(sub.notna() & (sub != h.country), h.country + ": " + sub.fillna(""), h.country)
print(f"human isolates with place below country level: {(h.place_sub != h.country).mean():.1%}")
def kstats(d, qi):
    g = d.groupby(qi, observed=True).target_acc.transform("size"); return dict(n=len(d), unique=(g == 1).mean() * 100, below5=(g < 5).mean() * 100, median_group=float(g.median()))
rows = []
day = h[h.date_precision == "day"].copy(); day["week"] = day.coll_day.dt.to_period("W").astype(str); day["month"] = day.coll_day.dt.to_period("M").astype(str); day["year_"] = day.coll_day.dt.year.astype(str); day["day"] = day.coll_day.astype(str)
mon = h[h.date_precision == "month"].copy(); mon["month"] = mon.coll_month.dt.to_period("M").astype(str); mon["year_"] = mon.coll_month.dt.year.astype(str)
sets = {"country + date": ["country"], "place below country + date": ["place_sub"], "place + date + serovar": ["place_sub", "serovar_clean"]}
# temporal error introduced by each date resolution, on the day-level records
mid = {"day": day.coll_day, "week": day.coll_day.dt.to_period("W").dt.start_time + pd.Timedelta(days=3.5), "month": day.coll_day.dt.to_period("M").dt.start_time + pd.Timedelta(days=14.5), "year_": pd.to_datetime(day.coll_day.dt.year.astype(str) + "-07-01")}
err = {k: float((day.coll_day - v).abs().dt.total_seconds().mean() / 86400) for k, v in mid.items()}
for sname, qi in sets.items():
    for res in ["day", "week", "month", "year_"]:
        r = kstats(day, qi + [res]); rows.append(dict(records="day-level dates", qi_set=sname, resolution=res.rstrip("_"), temporal_error_days=err[res], **r))
    for res in ["month", "year_"]:
        r = kstats(mon, qi + [res]); rows.append(dict(records="month-level dates", qi_set=sname, resolution=res.rstrip("_"), temporal_error_days=np.nan, **r))
T = pd.DataFrame(rows); T.round(2).to_csv("table_privacy.csv", index=False); print(T.round(1).to_string(index=False))
print("temporal error (mean days) by resolution:", {k: round(v, 1) for k, v in err.items()})

# ---------------------------------------------------------------- Figure 3
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.5), gridspec_kw={"width_ratios": [1, 0.8, 1]}); cols = {"country + date": "#a0aec0", "place below country + date": "#2b6cb0", "place + date + serovar": "#c53030"}; res = ["day", "week", "month", "year"]
for sname, c in cols.items():
    t = T[(T.records == "day-level dates") & (T.qi_set == sname)].set_index("resolution").reindex(res); ax[0].plot(res, t.unique, marker="o", color=c, label=sname)
ax[0].set(ylabel="Records that are unique (%)", title=f"Day-level records (n={len(day):,})", ylim=(0, None)); ax[0].legend(frameon=False, fontsize=8.5)
for sname, c in cols.items():
    t = T[(T.records == "month-level dates") & (T.qi_set == sname)].set_index("resolution").reindex(["month", "year"]); ax[1].plot(["month", "year"], t.unique, marker="o", color=c)
ax[1].set(title=f"Month-level records (n={len(mon):,})", ylim=(0, None))
t = T[(T.records == "day-level dates") & (T.qi_set == "place + date + serovar")].set_index("resolution").reindex(res); ax[2].plot(t.temporal_error_days, t.unique, marker="o", color="#c53030")
for r_ in res: ax[2].annotate(r_, (t.temporal_error_days[r_], t.unique[r_]), xytext=(5, 6), textcoords="offset points", fontsize=9)
ax[2].set(xlabel="Mean temporal error from coarsening the date (days)", ylabel="Records that are unique (%)", title="Privacy against temporal detail", ylim=(0, None)); fig.tight_layout(); fig.savefig("figs12/fig3_privacy.png", dpi=200); plt.close(fig); print("figure written")
