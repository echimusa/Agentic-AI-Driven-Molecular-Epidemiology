"""
Chapter 12, Hands-on 8 - What an embargo costs: the detection delay of US Salmonella clusters when public release is held back.
Same clusters as Hands-on 4 (US SNP clusters of 10 or more dated isolates, 2019-2022, month-level collection dates placed at day 15). A cluster is detected when its third isolate is public. An embargo of E days adds E days to every
release date. For each cluster the script also counts the isolates collected by the time of detection beyond the third, a floor on the cases that accrued while the cluster was unseen.
Input: sal_surv.csv.gz, cluster_list_us.tsv. Output: table_embargo.csv, figs12/fig8_embargo.png.
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
d["is_month"] = d.date_precision == "month"; d = d[d.groupby("PDS_acc").target_acc.transform("size") >= 10]
groups = {k: (g.base.values + np.where(g.is_month.values, 14, 0), g.rel.values) for k, g in d.groupby("PDS_acc")}
K = 3
def run(E):
    out = []
    for k_, (c, r) in groups.items():
        ideal = np.sort(c)[K - 1]; obs = np.sort(r)[K - 1] + E; out.append((ideal, max(obs - ideal, 0), max(int((c <= obs).sum()) - K, 0)))
    x = pd.DataFrame(out, columns=["ideal", "delay", "accrued"]); x["date"] = pd.to_datetime(x.ideal, unit="D"); return x[(x.date >= "2019-01-01") & (x.date <= "2022-12-31")]
rows = []
for E in (0, 30, 90, 180, 365):
    x = run(E); rows.append(dict(embargo_days=E, clusters=len(x), median_delay_days=x.delay.median(), p90_delay_days=x.delay.quantile(0.9), detected_within_30_days_pct=100 * (x.delay <= 30).mean(), median_isolates_collected_before_detection=x.accrued.median(), mean_isolates_collected_before_detection=x.accrued.mean()))
T = pd.DataFrame(rows).set_index("embargo_days"); T.round(1).to_csv("table_embargo.csv"); print(T.round(1).to_string())
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.1, 0.9, 0.9]})
for E, col in ((0, "#2b6cb0"), (30, "#d69e2e"), (90, "#dd6b20"), (180, "#c53030")):
    v = np.sort(run(E).delay.values); ax[0].plot(v, np.arange(1, len(v) + 1) / len(v) * 100, color=col, lw=2, label=f"embargo {E} days")
ax[0].set(xlim=(0, 400), ylim=(0, 100), xlabel="Detection delay of the third isolate (days)", ylabel="Cumulative % of clusters", title="A. Delay by embargo"); ax[0].legend(frameon=False, fontsize=8.8, loc="lower right")
ax[1].bar([str(e) for e in T.index], T.median_delay_days, color="#a0aec0"); ax[1].set(xlabel="Embargo (days)", ylabel="Median delay (days)", title="B. Median delay", ylim=(0, None))
ax[2].bar([str(e) for e in T.index], T.mean_isolates_collected_before_detection, color="#2b6cb0"); ax[2].set(xlabel="Embargo (days)", ylabel="Mean isolates collected before detection\n(beyond the third)", title="C. Cases that accrue unseen", ylim=(0, None))
fig.tight_layout(); fig.savefig("figs12/fig8_embargo.png", dpi=200); plt.close(fig); print("figure written")
