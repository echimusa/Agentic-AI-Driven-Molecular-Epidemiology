"""
Chapter 7, Exercise 3 - Describe an outbreak in time, place and person, and size a tracing window (Hagelloch measles, 1861).
Same public data as Exercise 1 (CRAN package `outbreaks`; documented infector for 184 of 188 cases). Plain Python, no R.
1. Time: the epidemic curve by prodrome date, and the generation of each case (steps back along the documented infector links to a case without a recorded infector).
2. Person: age, sex and class of the cases, and the recorded fatality by age group (deaths with a date of death in the data; Wilson 95% intervals).
3. Place: cases per family and how often the documented infector shares the infectee's family or school class, and how far apart the houses are.
4. Case definition: how the choice of onset date (prodrome or rash) shifts the curve and changes the serial interval.
5. Tracing: for a look-back window of W days and a definition of contact, what share of documented infectors a tracer would find, and how many candidates per case that costs.
Writes table_generations.csv, table_person.csv, table_place.csv, table_tracing.csv, table_descriptive_summary.csv and figs3/fig11_descriptive.png.
"""
import os, math
import numpy as np, pandas as pd, pyreadr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs3", exist_ok=True)
d = pyreadr.read_r("measles_hagelloch_1861.RData")["measles_hagelloch_1861"].copy()
d["onset"] = pd.to_datetime(d.date_of_prodrome); d["rash"] = pd.to_datetime(d.date_of_rash); d["t"] = (d.onset - d.onset.min()).dt.days; d["infector"] = pd.to_numeric(d.infector); d = d.set_index("case_ID", drop=False)
print(f"{len(d)} cases, prodrome {d.onset.min().date()} to {d.onset.max().date()} ({d.t.max()} days), documented infector for {d.infector.notna().sum()}")
# ---------------------------------------------------------------- 1. time: epidemic curve and generations
gen = {}
def depth(c, seen=()):
    if c in gen: return gen[c]
    inf = d.infector.get(c)
    g = 0 if (pd.isna(inf) or int(inf) not in d.index or int(inf) in seen) else depth(int(inf), seen + (c,)) + 1
    gen[c] = g; return g
for c in d.index: depth(int(c))
d["generation"] = d.index.map(gen)
G = d.groupby("generation").agg(cases=("t", "size"), median_day=("t", "median"), first_day=("t", "min"), last_day=("t", "max")).reset_index(); G["gap_to_previous_median"] = G.median_day.diff(); G.round(2).to_csv("table_generations.csv", index=False)
pairs = d[d.infector.notna()].copy(); pairs["si"] = pairs.t.values - d.t.reindex(pairs.infector.astype(int)).values
print(f"time: peak day {int(d.t.value_counts().idxmax())} ({int(d.t.value_counts().max())} cases); {len(G)} generations; median onset day by generation {G.median_day.round(0).astype(int).tolist()}; mean gap between generation medians {G.gap_to_previous_median.mean():.1f} d")
# ---------------------------------------------------------------- 2. person
def wilson(k, n, z=1.96):
    p = k / n; den = 1 + z * z / n; c = (p + z * z / (2 * n)) / den; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den; return max(0.0, c - h), min(1.0, c + h)
d["dead"] = d.date_of_death.notna(); d["age_group"] = pd.cut(d.age, [-1, 4, 9, 100], labels=["0-4", "5-9", "10+"])
rows = [dict(group="all", n=len(d), deaths=int(d.dead.sum()))]
for g, s in d.groupby("age_group", observed=True): rows.append(dict(group=f"age {g}", n=len(s), deaths=int(s.dead.sum())))
for g, s in d.groupby("gender"): rows.append(dict(group=f"sex {g}", n=len(s), deaths=int(s.dead.sum())))
for g, s in d.groupby("class"): rows.append(dict(group=f"class {g}", n=len(s), deaths=int(s.dead.sum())))
P = pd.DataFrame(rows); P["cfr"] = P.deaths / P.n; P["cfr_lo"], P["cfr_hi"] = zip(*[wilson(k, n) for k, n in zip(P.deaths, P.n)]); P.round(4).to_csv("table_person.csv", index=False)
print(f"person: median age {d.age.median():.0f} y (IQR {d.age.quantile(.25):.0f}-{d.age.quantile(.75):.0f}); recorded deaths {int(d.dead.sum())} of {len(d)} ({100 * d.dead.mean():.1f}%); " + "; ".join(f"age {r.group[4:]} {100 * r.cfr:.1f}%" for r in P[P.group.str.startswith('age')].itertuples()))
# ---------------------------------------------------------------- 3. place
fam = d.groupby("family_ID").size(); multi = fam[fam >= 2]
fid, cls = d.family_ID, d["class"]; xy = d[["x_loc", "y_loc"]]
pairs["same_family"] = (fid.reindex(pairs.index).values == fid.reindex(pairs.infector.astype(int)).values); pairs["same_class"] = (cls.reindex(pairs.index).values == cls.reindex(pairs.infector.astype(int)).values)
pairs["dist"] = np.hypot(xy.reindex(pairs.index).x_loc.values - xy.reindex(pairs.infector.astype(int)).x_loc.values, xy.reindex(pairs.index).y_loc.values - xy.reindex(pairs.infector.astype(int)).y_loc.values)
allp = [(i, j) for j in d.index for i in d.index if d.t[i] < d.t[j]]; base_same_family = np.mean([fid[i] == fid[j] for i, j in allp]); base_same_class = np.mean([cls[i] == cls[j] for i, j in allp])
PL = pd.DataFrame([dict(measure="families with cases", value=len(fam)), dict(measure="families with 2 or more cases", value=len(multi)), dict(measure="share of cases in those families", value=multi.sum() / len(d)), dict(measure="documented pairs in the same family", value=pairs.same_family.mean()), dict(measure="all earlier-later pairs in the same family", value=base_same_family), dict(measure="documented pairs in the same class", value=pairs.same_class.mean()), dict(measure="all earlier-later pairs in the same class", value=base_same_class), dict(measure="median distance, documented pairs (map units)", value=pairs.dist.median()), dict(measure="median distance, all earlier-later pairs (map units)", value=float(np.median([np.hypot(xy.x_loc[i] - xy.x_loc[j], xy.y_loc[i] - xy.y_loc[j]) for i, j in allp])))]); PL.round(4).to_csv("table_place.csv", index=False)
print("place: " + "; ".join(f"{r.measure} {r.value:.3f}" for r in PL.itertuples()))
# ---------------------------------------------------------------- 4. case definition: prodrome or rash as onset
lag = (d.rash - d.onset).dt.days; sir = (d.rash.reindex(pairs.index).values - d.rash.reindex(pairs.infector.astype(int)).values) / np.timedelta64(1, "D")
print(f"onset: rash follows prodrome by a median of {lag.median():.0f} d (IQR {lag.quantile(.25):.0f}-{lag.quantile(.75):.0f}); serial interval {pairs.si.mean():.1f} d (SD {pairs.si.std():.1f}) by prodrome, {np.nanmean(sir):.1f} d (SD {np.nanstd(sir, ddof=1):.1f}) by rash; peak day by rash {int(((d.rash - d.onset.min()).dt.days).value_counts().idxmax())} against {int(d.t.value_counts().idxmax())} by prodrome")
# ---------------------------------------------------------------- 5. tracing windows and contact definitions
defs = {"anyone": lambda i, j: True, "same family": lambda i, j: fid[i] == fid[j], "same class": lambda i, j: cls[i] == cls[j], "family or class": lambda i, j: fid[i] == fid[j] or cls[i] == cls[j], "within 25 map units": lambda i, j: math.hypot(xy.x_loc[i] - xy.x_loc[j], xy.y_loc[i] - xy.y_loc[j]) <= 25, "within 100 map units": lambda i, j: math.hypot(xy.x_loc[i] - xy.x_loc[j], xy.y_loc[i] - xy.y_loc[j]) <= 100}
T_ = []
for W in (7, 10, 14, 21):
    for name, f in defs.items():
        hit, ncand = 0, []
        for j in d.index:
            cand = [i for i in d.index if 0 < d.t[j] - d.t[i] <= W and f(i, j)]; ncand.append(len(cand))
            if pd.notna(d.infector[j]) and int(d.infector[j]) in cand: hit += 1
        T_.append(dict(window_days=W, contact_definition=name, recall=hit / d.infector.notna().sum(), candidates_per_case=float(np.mean(ncand))))
T_ = pd.DataFrame(T_); T_["recall_per_candidate"] = T_.recall / T_.candidates_per_case; T_.round(4).to_csv("table_tracing.csv", index=False)
q = pairs.si.quantile([.5, .9, .95, .99]); print(f"tracing: serial interval median {q[.5]:.0f}, 90% {q[.9]:.0f}, 95% {q[.95]:.0f}, 99% {q[.99]:.0f} d; max {int(pairs.si.max())} d")
print(T_[T_.window_days == 14].round(3).to_string(index=False))
# ---------------------------------------------------------------- summary numbers quoted in the chapter
SM = dict(cases=len(d), onset_first=str(d.onset.min().date()), onset_last=str(d.onset.max().date()), span_days=int(d.t.max()), documented_infectors=int(d.infector.notna().sum()), peak_day=int(d.t.value_counts().idxmax()), peak_cases=int(d.t.value_counts().max()), generations=int(len(G)),
    median_age=float(d.age.median()), age_q1=float(d.age.quantile(.25)), age_q3=float(d.age.quantile(.75)), recorded_deaths=int(d.dead.sum()), rash_lag_median=float(lag.median()), rash_lag_q1=float(lag.quantile(.25)), rash_lag_q3=float(lag.quantile(.75)),
    si_mean_prodrome=float(pairs.si.mean()), si_sd_prodrome=float(pairs.si.std()), si_mean_rash=float(np.nanmean(sir)), si_sd_rash=float(np.nanstd(sir, ddof=1)), peak_day_rash=int(((d.rash - d.onset.min()).dt.days).value_counts().idxmax()),
    si_q50=float(q[.5]), si_q90=float(q[.9]), si_q95=float(q[.95]), si_q99=float(q[.99]), si_max=int(pairs.si.max()))
pd.DataFrame({"key": list(SM), "value": list(SM.values())}).to_csv("table_descriptive_summary.csv", index=False)
# ---------------------------------------------------------------- summary of the single numbers quoted in the chapter
sir_ = (d.rash.reindex(pairs.index).values - d.rash.reindex(pairs.infector.astype(int)).values) / np.timedelta64(1, "D")
summ = [("cases", len(d)), ("onset_first", str(d.onset.min().date())), ("onset_last", str(d.onset.max().date())), ("span_days", int(d.t.max())), ("documented_infectors", int(d.infector.notna().sum())), ("peak_day", int(d.t.value_counts().idxmax())), ("peak_cases", int(d.t.value_counts().max())), ("generations", len(G)), ("median_age", float(d.age.median())), ("age_q1", float(d.age.quantile(.25))), ("age_q3", float(d.age.quantile(.75))), ("recorded_deaths", int(d.dead.sum())), ("rash_lag_median", float(lag.median())), ("rash_lag_q1", float(lag.quantile(.25))), ("rash_lag_q3", float(lag.quantile(.75))), ("si_mean_prodrome", float(pairs.si.mean())), ("si_sd_prodrome", float(pairs.si.std())), ("si_mean_rash", float(np.nanmean(sir_))), ("si_sd_rash", float(np.nanstd(sir_, ddof=1))), ("peak_day_rash", int(((d.rash - d.onset.min()).dt.days).value_counts().idxmax())), ("si_q50", float(q[.5])), ("si_q90", float(q[.9])), ("si_q95", float(q[.95])), ("si_q99", float(q[.99])), ("si_max", int(pairs.si.max()))]
pd.DataFrame(summ, columns=["key", "value"]).to_csv("table_descriptive_summary.csv", index=False)
# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2)); cols = ["#0b3c49", "#1b998b", "#f3a712", "#e4572e", "#3d7ea6", "#2e294e", "#718096"]
bot = np.zeros(d.t.max() + 1)
for g in sorted(d.generation.unique()):
    c = np.bincount(d[d.generation == g].t, minlength=d.t.max() + 1); ax[0].bar(range(len(c)), c, bottom=bot, color=cols[g % len(cols)], width=1.0, label=f"generation {g}"); bot += c
ax[0].set(xlabel="Day since the first prodrome", ylabel="Cases per day", title="Epidemic curve by generation"); ax[0].legend(frameon=False, fontsize=8, ncol=2)
a_ = P[P.group.str.startswith("age")]; ax[1].bar(a_.group.str[4:], 100 * a_.cfr, color="#e4572e", yerr=[100 * (a_.cfr - a_.cfr_lo), 100 * (a_.cfr_hi - a_.cfr)], capsize=4); ax[1].set(xlabel="Age group (years)", ylabel="Recorded fatality (%)", title="Person: recorded fatality by age")
s14 = T_[T_.window_days == 14]; ax[2].scatter(s14.candidates_per_case, 100 * s14.recall, color="#1b998b", s=50)
for r in s14.itertuples(): ax[2].annotate(r.contact_definition, (r.candidates_per_case, 100 * r.recall), fontsize=8, xytext=(4, -9), textcoords="offset points")
ax[2].set(xscale="log", xlabel="Candidate infectors per case (log scale)", ylabel="Documented infectors found (%)", title="Tracing a 14-day look-back"); fig.tight_layout(); fig.savefig("figs3/fig11_descriptive.png", dpi=180)
