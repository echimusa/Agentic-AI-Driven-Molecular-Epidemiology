"""
Chapter 11, Hands-on 6 - vitamin D and hospitalised COVID-19: the "think critically" question of Session 11, with the same estimators as Hands-on 1 and 2.
Exposure: serum 25-hydroxyvitamin D (standardised, log scale) in 443,734 Europeans (Manousaki et al. 2020, GWAS Catalog GCST010144), the same study that Butler-Laporte et al. (2021) used. Instruments are the 138 conditionally independent
variants the study reports with P < 5e-8; their effects, standard errors and frequencies are read from the harmonised summary-statistics file on the EBI server with tabix, one position at a time.
Outcome: COVID-19 Host Genetics Initiative release 5, hospitalised vs population (GCST011081), read from hgi_r5_slim.npz of Chapter 10.
Offline from vitd_instruments.csv and vitd_hgi.csv. To rebuild them, put the catalogue associations table (assoc_full.tsv) and the npz in this folder, delete the two csv files and keep the network on (tabix and bcftools needed).
Output: Figure 5, table_vitd_methods.csv, table_vitd_instruments.csv (the five strongest), table_vitd_summary.csv.   Needs numpy, pandas, scipy, matplotlib and mr11.py.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010144/harmonised/GCST010144.h.tsv.gz"
if not os.path.exists("vitd_instruments.csv"):
    import subprocess
    d = pd.read_csv("assoc_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["PUBMEDID", "DISEASE/TRAIT", "SNPS", "CHR_ID", "CHR_POS", "P-VALUE"]); d = d[(d.PUBMEDID == "32059762") & (d["DISEASE/TRAIT"] == "Serum 25-Hydroxyvitamin D levels")].copy()
    d["px"] = d["P-VALUE"].astype(float); d = d[(d.px < 5e-8) & d.SNPS.str.fullmatch(r"rs\d+")].drop_duplicates("SNPS"); d[["CHR_ID", "CHR_POS", "CHR_POS"]].to_csv("vitd_regions.tsv", sep="\t", header=False, index=False)
    out = subprocess.run(["tabix", URL, "-R", "vitd_regions.tsv"], capture_output=True, text=True).stdout; want = set(d.SNPS); rows = []        # one call, so the 1.5 MB index is fetched once
    for ln in out.splitlines():
        f = ln.split("\t")
        if f[12].strip() in want: rows.append(dict(rsid=f[12].strip(), chr=f[0], pos=int(f[1]), ea=f[2].upper(), oa=f[3].upper(), bx=float(f[4]), sx=float(f[5]), eaf_x=float(f[6]), px=float(f[7]), n=int(f[9])))
    pd.DataFrame(rows).to_csv("vitd_instruments.csv", index=False)
if not os.path.exists("vitd_hgi.csv"):
    z = np.load("hgi_r5_slim.npz"); X0 = pd.read_csv("vitd_instruments.csv"); rs = z["rsid"]; m = np.isin(rs, X0.rsid.values)
    pd.DataFrame({"rsid": rs[m], "oa1": z["a1"][m], "oa2": z["a2"][m], "by": z["beta"][m], "sy": z["se"][m], "eaf_y": z["eaf"][m]}).to_csv("vitd_hgi.csv", index=False)
X = pd.read_csv("vitd_instruments.csv"); Y = pd.read_csv("vitd_hgi.csv").drop_duplicates("rsid"); print(f"{len(X)} instruments read from the exposure file; {len(X.merge(Y, on='rsid'))} found in the COVID-19 file")
D = X.merge(Y, on="rsid"); f = D.eaf_x.clip(0.001, 0.999)
keep, sign = mr.harmonise(D.ea.values, D.oa.values, D.eaf_x.values, D.oa1.values, D.oa2.values, D.eaf_y.values); D = D[keep].copy(); D["by"] = D["by"].values * sign[keep]
bx, sx, by, sy = D.bx.values, D.sx.values, D.by.values, D.sy.values; neg = bx < 0; bx, by = np.abs(bx), np.where(neg, -by, by)       # orient so that every effect on the exposure is positive
D["bx"], D["by"] = bx, by; D["F"] = mr.fstat(bx, sx); print(f"{len(D)} instruments after harmonisation; mean F {D.F.mean():.0f}, minimum F {D.F.min():.0f}; variance explained {np.sum(2 * f[keep].values * (1 - f[keep].values) * bx ** 2):.4f} (standardised trait)")
bw, sw = mr.wald(bx, by, sy); I = mr.ivw(bx, by, sy); E = mr.egger(bx, by, sy); W = mr.weighted_median(bw, sw, nboot=1000); M = mr.weighted_mode(bw, sw, nboot=500); loo = mr.leave_one_out(bx, by, sy)
o = lambda b, se: (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se))
rows = [dict(method="IVW, random effects", beta=I["beta"], se=I["se_random"], assumes="No directional pleiotropy"), dict(method="MR-Egger", beta=E["slope"], se=E["se_slope"], assumes="Pleiotropy independent of instrument strength"),
        dict(method="Weighted median", beta=W["beta"], se=W["se"], assumes="At least half of the weight from valid instruments"), dict(method="Weighted mode", beta=M["beta"], se=M["se"], assumes="The largest group of instruments is valid")]
R = pd.DataFrame(rows); R[["odds_ratio", "ci_low", "ci_high"]] = [o(b, s) for b, s in zip(R.beta, R.se)]; R["p"] = 2 * __import__("scipy").stats.norm.sf(np.abs(R.beta / R.se)); R.to_csv("table_vitd_methods.csv", index=False); print(R.round(3).to_string(index=False))
power_or = np.exp(2.8 * I["se_random"]); S = pd.DataFrame([dict(instruments=len(D), mean_F=D.F.mean(), min_F=D.F.min(), Q=I["Q"], Q_df=I["df"], Q_p=I["Q_p"], I2=I["I2"], egger_intercept=E["intercept"], egger_intercept_se=E["se_intercept"], egger_intercept_p=E["p_intercept"],
                                 loo_min=loo.min(), loo_max=loo.max(), smallest_or_80pct_power=power_or, smallest_protective_or_80pct_power=1 / power_or, share_same_sign=float(np.mean(np.sign(bw) == np.sign(I["beta"]))), n_p_lt_0_05=int((2 * __import__("scipy").stats.norm.sf(np.abs(by / sy)) < 0.05).sum()))])
S.to_csv("table_vitd_summary.csv", index=False); print(S.round(4).T.to_string(header=False))
top = D.sort_values("F", ascending=False).head(5)[["rsid", "chr", "pos", "ea", "eaf_x", "bx", "F", "by", "sy"]].copy(); top["wald_or"] = np.exp(top.by / top.bx); top.to_csv("table_vitd_instruments.csv", index=False); print(top.round(3).to_string(index=False))
# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2)); lim = max(bx) * 1.05
ax[0].errorbar(bx, by, yerr=1.96 * sy, fmt="o", ms=3.2, color="#3D7EA6", ecolor="#BBD3E3", elinewidth=0.8, alpha=0.9); xs = np.array([0, lim])
for b_, c_, l_ in ((I["beta"], "#E4572E", "IVW"), (E["slope"], "#1B998B", "MR-Egger"), (M["beta"], "#2E294E", "Weighted mode")): ax[0].plot(xs, b_ * xs + (E["intercept"] if l_ == "MR-Egger" else 0), color=c_, lw=1.8, label=l_)
ax[0].axhline(0, color="#a0aec0", lw=0.8); ax[0].set(xlabel="Effect on 25-hydroxyvitamin D (catalogue units)", ylabel="Effect on hospitalised COVID-19 (log odds ratio)", title="Each point is one of the instruments"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper left")
prec = np.abs(bx) / sy; ax[1].scatter(bw, prec, s=14, color="#3D7EA6", alpha=0.85); ax[1].axvline(I["beta"], color="#E4572E", lw=1.5); ax[1].axvline(0, color="#a0aec0", lw=0.8); ax[1].set(xlabel="Wald ratio (log odds ratio)", ylabel="Precision of the ratio", title="Funnel plot"); ax[1].set_xlim(-6, 6)
k = np.argsort(loo)[::-1] if False else np.arange(len(loo)); ax[2].plot(np.sort(loo), np.arange(len(loo)), "o", ms=3, color="#3D7EA6"); ax[2].axvline(I["beta"], color="#E4572E", lw=1.5); ax[2].axvline(0, color="#a0aec0", lw=0.8); ax[2].set(xlabel="IVW estimate leaving one instrument out", ylabel="Instruments (sorted)", title="Leave-one-out")
fig.tight_layout(); fig.savefig("figs11/fig5_vitd.png", dpi=200); print("figure written")
