"""
Chapter 13, Hands-on 6 - Different windows, different answers: the Session 13 think-critically question on the Ebola line list.
A. Cori reproduction number with trailing windows of 7, 14 and 28 days: the value on 1 October 2014, the width of its interval, and the first date after the incidence peak from which the point estimate stays below 1 for 14 days.
B. The serial interval: mean 10, 15.3 and 20 days with the coefficient of variation held fixed.
C. A falling reporting fraction: expected case counts thinned from 100% on 1 October 2014 to 50% on 31 October (a hypothetical, not a feature of the data), and what it does to the estimate.
Input: ebola_sl_linelist.csv. Output: table_windows.csv, table_si_sensitivity.csv, table_reporting_fall.csv, figs13/fig6_windows.png.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ebolacore import gamma_si, cori_rt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
d = pd.read_csv("ebola_sl_linelist.csv", parse_dates=["date_of_onset"]); day = d.groupby("date_of_onset").size().asfreq("D", fill_value=0); inc = day.values.astype(float); peak = day.idxmax()
def rt_frame(x, w, window):
    R = cori_rt(x, w, window=window); idx = day.index[R[:, 0].astype(int)]; return pd.DataFrame({"Rt": R[:, 1], "lo": R[:, 2], "hi": R[:, 3]}, index=idx)
def sustained_below(F, n=14):
    F = F[peak:]                                                               # only dates from the incidence peak onwards
    r = (F.Rt < 1).astype(int).rolling(n).sum(); ok = r[r == n].index; return (ok.min() - pd.Timedelta(days=n - 1)) if len(ok) else pd.NaT
w0 = gamma_si(); D0 = pd.Timestamp("2014-10-01"); rows = []; fr = {}
for win in (7, 14, 28):
    F = rt_frame(inc, w0, win); fr[win] = F; s = sustained_below(F)
    rows.append(dict(window_days=win, Rt_1_Oct=F.loc[D0, "Rt"], lo=F.loc[D0, "lo"], hi=F.loc[D0, "hi"], interval_width=F.loc[D0, "hi"] - F.loc[D0, "lo"], Rt_1_Sep=F.loc["2014-09-01", "Rt"], first_date_below_1_for_14_days=s.strftime("%Y-%m-%d"), days_after_peak=(s - peak).days))
A = pd.DataFrame(rows); A.to_csv("table_windows.csv", index=False); print(f"peak day {peak:%d %b %Y}"); print(A.round(3).to_string(index=False))
rows = []
for m in (10.0, 15.3, 20.0):
    F = rt_frame(inc, gamma_si(m, m * 9.3 / 15.3), 14); rows.append(dict(mean_interval_days=m, sd_days=round(m * 9.3 / 15.3, 2), Rt_1_Oct=F.loc[D0, "Rt"], lo=F.loc[D0, "lo"], hi=F.loc[D0, "hi"], Rt_1_Nov=F.loc["2014-11-01", "Rt"]))
B = pd.DataFrame(rows); B.to_csv("table_si_sensitivity.csv", index=False); print(B.round(3).to_string(index=False))
frac = np.ones(len(day)); t0 = (D0 - day.index[0]).days; frac[t0:t0 + 31] = np.linspace(1.0, 0.5, 31); frac[t0 + 31:] = 0.5
thin = inc * frac; Fo = fr[14]; Ft = rt_frame(thin, w0, 14); rows = []
for lab in ("2014-10-15", "2014-11-01", "2014-11-15", "2014-12-01"):
    rows.append(dict(date=lab, reporting_fraction=frac[(pd.Timestamp(lab) - day.index[0]).days], Rt_all_cases=Fo.loc[lab, "Rt"], Rt_thinned=Ft.loc[lab, "Rt"], ratio=Ft.loc[lab, "Rt"] / Fo.loc[lab, "Rt"]))
Cc = pd.DataFrame(rows); Cc.to_csv("table_reporting_fall.csv", index=False); print(Cc.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.2, 0.9, 1.0]})
for win, col in ((7, "#c53030"), (14, "#2b6cb0"), (28, "#2f855a")):
    F = fr[win]["2014-08-15":"2015-01-31"]; ax[0].plot(F.index, F.Rt, color=col, lw=1.6, label=f"{win}-day window")
ax[0].axhline(1, color="black", ls="--", lw=0.8); ax[0].axvline(peak, color="grey", ls=":", lw=1); ax[0].text(peak, 3.15, " incidence peak", fontsize=8.5, color="grey"); ax[0].set(ylim=(0, 3.4), ylabel="Instantaneous reproduction number", title="A. The window sets the lag"); ax[0].legend(frameon=False, fontsize=8.8); ax[0].tick_params(axis="x", rotation=30)
for k, (m, col) in enumerate(zip((10.0, 15.3, 20.0), ("#d69e2e", "#2b6cb0", "#2f855a"))): r_ = B.iloc[k]; ax[1].errorbar([k], [r_.Rt_1_Oct], yerr=[[r_.Rt_1_Oct - r_.lo], [r_.hi - r_.Rt_1_Oct]], fmt="o", color=col, ms=7, capsize=4)
ax[1].set(xticks=range(3), xticklabels=["10", "15.3", "20"], xlim=(-0.5, 2.5), xlabel="Mean serial interval (days)", ylabel="Rt on 1 October 2014 (14-day window)", title="B. The serial interval"); ax[1].axhline(1, color="black", ls="--", lw=0.8)
G = Fo["2014-09-01":"2014-12-31"]; H = Ft["2014-09-01":"2014-12-31"]; ax[2].plot(G.index, G.Rt, color="#2b6cb0", lw=1.8, label="All reported cases"); ax[2].plot(H.index, H.Rt, color="#c53030", lw=1.8, label="Reporting falls to 50% in October"); ax[2].axhline(1, color="black", ls="--", lw=0.8)
ax[2].set(ylim=(0, 2.5), ylabel="Instantaneous reproduction number", title="C. A falling reporting fraction"); ax[2].legend(frameon=False, fontsize=8.8); ax[2].tick_params(axis="x", rotation=30)
fig.tight_layout(); fig.savefig("figs13/fig6_windows.png", dpi=200); plt.close(fig); print("figure written")
