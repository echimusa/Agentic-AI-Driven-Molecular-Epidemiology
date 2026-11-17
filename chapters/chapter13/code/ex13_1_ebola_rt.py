"""
Chapter 13, Hands-on 1 - Reproduction number and growth rate from public case data: the 2014 Ebola outbreak in Sierra Leone.
Data: the 11,903-case line list of the CRAN package outbreaks (ebola_sierraleone_2014, exported to ebola_sl_linelist.csv). Serial interval: gamma, mean 15.3 d, SD 9.3 d.
Output: table_ebola_rt.csv, Figure 1.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ebolacore import *
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
d = pd.read_csv("ebola_sl_linelist.csv", parse_dates=["date_of_onset"]); day = d.groupby("date_of_onset").size().asfreq("D", fill_value=0)
print(f"{len(d):,} cases with onset {day.index.min():%d %b %Y} to {day.index.max():%d %b %Y}; peak day {day.idxmax():%d %b %Y} ({day.max()} cases); peak 7-day total {day.rolling(7).sum().max():.0f}")
inc = day.values.astype(float); w = gamma_si(); R = cori_rt(inc, w, window=14); idx = day.index[R[:, 0].astype(int)]
Rt = pd.DataFrame({"date": idx, "Rt": R[:, 1], "lo": R[:, 2], "hi": R[:, 3]}).set_index("date"); Rt.round(3).to_csv("table_ebola_rt.csv")
m = Rt.Rt.resample("MS").median(); print("median Rt by month:", {f"{k:%b %Y}": round(float(v), 2) for k, v in m["2014-06-01":"2015-03-01"].items()})
below = (Rt.Rt < 1).astype(int); run = below.rolling(28).sum(); stay = run[run == 28].index.min(); print(f"first date after which the Rt point estimate stays below 1 for 28 days: {stay:%d %b %Y}"); print(f"days with Rt interval entirely above 1: {(Rt.lo > 1).sum()}; entirely below 1: {(Rt.hi < 1).sum()}")
# exponential-growth window: the early epidemic, which is also when genomes were collected (June to mid-September 2014)
for lab, a, b in [("1 Jun to 15 Sep 2014", "2014-06-01", "2014-09-15"), ("1 Jun to 31 Aug 2014", "2014-06-01", "2014-08-31"), ("15 Jul to 15 Sep 2014", "2014-07-15", "2014-09-15")]:
    s = day[a:b]; r, se = poisson_growth(s.values); print(f"growth {lab}: r = {r:.4f}/day (SE {se:.4f}); doubling {np.log(2) / r:.1f} d (95% {np.log(2) / (r + 1.96 * se):.1f}-{np.log(2) / (r - 1.96 * se):.1f}); R = {r_to_R(r):.2f} ({r_to_R(r - 1.96 * se):.2f}-{r_to_R(r + 1.96 * se):.2f})")
r0, se0 = poisson_growth(day["2014-06-01":"2014-09-15"].values); np.save("r_cases.npy", np.array([r0, se0]))
print("sensitivity of R to the generation interval at this growth rate:", {f"mean {m} d": round(float(r_to_R(r0, m, m * 9.3 / 15.3)), 2) for m in (10, 12, 15.3, 18, 20)})

# rolling growth rate (21-day windows) for the right-hand panel
roll = []
for t in range(21, len(inc), 7):
    r, se = poisson_growth(inc[t - 21:t]); roll.append((day.index[t - 1], r, se))
roll = pd.DataFrame(roll, columns=["date", "r", "se"]).set_index("date")
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
ax[0].bar(day.index, day.values, width=1, color="#a0aec0"); ax[0].plot(day.index, day.rolling(7, center=True).mean(), color="#2b6cb0", lw=2); ax[0].axvspan(pd.Timestamp("2014-06-01"), pd.Timestamp("2014-09-15"), color="#fbd38d", alpha=0.4, lw=0)
ax[0].set(ylabel="Cases by onset date", title="Sierra Leone, 11,903 cases"); ax[0].tick_params(axis="x", rotation=30)
ax[1].fill_between(Rt.index, Rt.lo, Rt.hi, color="#2b6cb0", alpha=0.25); ax[1].plot(Rt.index, Rt.Rt, color="#2b6cb0", lw=1.6); ax[1].axhline(1, color="black", ls="--", lw=0.8); ax[1].set(ylabel="Instantaneous reproduction number", title="Rt, 14-day window", ylim=(0, 4)); ax[1].tick_params(axis="x", rotation=30)
ax[2].errorbar(roll.index, roll.r, yerr=1.96 * roll.se, fmt="o", ms=3, color="#2b6cb0", ecolor="#bee3f8"); ax[2].axhline(0, color="black", lw=0.8); ax[2].set(ylabel="Growth rate per day", title="Growth rate, 21-day windows"); ax[2].tick_params(axis="x", rotation=30)
fig.tight_layout(); fig.savefig("figs13/fig1_ebola_rt.png", dpi=200); plt.close(fig); print("figure written")
