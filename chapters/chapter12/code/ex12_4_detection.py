"""
Chapter 12, Hands-on 4 - How reporting delay postpones cluster detection.
For US SNP clusters of 10 or more isolates, a cluster is "detected" when its k-th isolate is visible. The ideal date is the k-th collection date, as if every isolate were shared the
day it was collected; the observed date is the k-th public release date. Month-level collection dates are placed at day 1, 15 or 28 of the month to bound the error this causes.
Input: sal_surv.pkl and the Pathogen Detection SNP-cluster list. Output: table_detection.csv and Figure 4.
"""
import os
import numpy as np, pandas as pd
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
CL = "../ch9/salmonella_cluster_list.tsv" if os.path.exists("../ch9/salmonella_cluster_list.tsv") else "cluster_list_us.tsv"
m = load_surv(); cl = pd.read_csv(CL, sep="\t", usecols=["PDS_acc", "target_acc"])
d = m[(m.country == "USA") & m.date_precision.isin(["day", "month"]) & m.release.notna()].copy(); d = d.merge(cl, on="target_acc")
d["rel"] = d.release.values.astype("datetime64[D]").astype(np.int64); base = np.where(d.date_precision == "day", d.coll_day, d.coll_month); d["base"] = pd.to_datetime(base).values.astype("datetime64[D]").astype(np.int64)
d["is_month"] = d.date_precision == "month"; d = d[d.groupby("PDS_acc").target_acc.transform("size") >= 10]; print(f"{d.PDS_acc.nunique():,} US clusters with at least 10 dated isolates; {len(d):,} isolates; month-level dates: {d.is_month.mean():.1%}")
groups = {k: (g.base.values, g.is_month.values, g.rel.values) for k, g in d.groupby("PDS_acc")}
def run(K, offset):
    out = []
    for k_, (b, mo, r) in groups.items():
        c = b + np.where(mo, offset, 0); ideal = np.sort(c)[K - 1]; obs = np.sort(r)[K - 1]; out.append((ideal, max(obs - ideal, 0)))
    x = pd.DataFrame(out, columns=["ideal", "delay"]); x["date"] = pd.to_datetime(x.ideal, unit="D"); return x[(x.date >= "2019-01-01") & (x.date <= "2022-12-31")]
rows = []
for K in (2, 3, 5, 10):
    x = run(K, 14); rows.append(dict(setting=f"k = {K}, month dates at day 15", clusters=len(x), median_delay_days=x.delay.median(), p90_delay_days=x.delay.quantile(0.9), within_14_days=(x.delay <= 14).mean() * 100, within_30_days=(x.delay <= 30).mean() * 100))
for off, lab in [(0, "day 1"), (28, "day 29")]:
    x = run(3, off); rows.append(dict(setting=f"k = 3, month dates at {lab}", clusters=len(x), median_delay_days=x.delay.median(), p90_delay_days=x.delay.quantile(0.9), within_14_days=(x.delay <= 14).mean() * 100, within_30_days=(x.delay <= 30).mean() * 100))
T = pd.DataFrame(rows).set_index("setting"); T.round(1).to_csv("table_detection.csv"); print(T.round(1).to_string())
X3 = run(3, 14)
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.5), gridspec_kw={"width_ratios": [1.15, 0.9, 0.9]})
v = np.sort(X3.delay.values); ax[0].plot(v, np.arange(1, len(v) + 1) / len(v) * 100, color="#2b6cb0", lw=2)
for q, lab in [(0.5, "median"), (0.9, "90th percentile")]: ax[0].axvline(np.quantile(v, q), color="#c53030", ls="--", lw=1); ax[0].text(np.quantile(v, q) + 4, 12 if q == 0.5 else 30, f"{lab}: {np.quantile(v, q):.0f} days", fontsize=8.8, color="#c53030")
ax[0].set(xlim=(0, 365), ylim=(0, 100), xlabel="Detection delay of the third isolate (days)", ylabel="Cumulative % of clusters", title=f"{len(X3)} US clusters, ideal detection in 2019-2022")
ks = [2, 3, 5, 10]; ax[1].bar([str(k) for k in ks], [T.loc[f"k = {k}, month dates at day 15", "median_delay_days"] for k in ks], color="#2b6cb0"); ax[1].set(xlabel="Isolates needed to detect (k)", ylabel="Median delay (days)", title="By detection threshold", ylim=(0, None))
lab = ["day 1", "day 15", "day 29"]; ax[2].bar(lab, [T.loc["k = 3, month dates at day 1", "median_delay_days"], T.loc["k = 3, month dates at day 15", "median_delay_days"], T.loc["k = 3, month dates at day 29", "median_delay_days"]], color="#a0aec0"); ax[2].set(xlabel="Where month-level dates are placed", ylabel="Median delay (days)", title="Sensitivity to date precision", ylim=(0, None))
fig.tight_layout(); fig.savefig("figs12/fig4_detection.png", dpi=200); plt.close(fig); print("figure written")
