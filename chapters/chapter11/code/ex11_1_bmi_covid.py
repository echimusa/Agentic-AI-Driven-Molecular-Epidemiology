"""
Chapter 11, Hands-on 1 and 2 - Two-sample Mendelian randomisation of body mass index on hospitalised COVID-19, from public summary statistics.
Exposure: BMI, GIANT + UK Biobank meta-analysis of Pulit et al. (GWAS Catalog GCST009004, European ancestry); instruments are the BMI variants reported by the authors with P < 5e-8.
Outcome: COVID-19 Host Genetics Initiative release 5, hospitalised vs population (GCST011081), read from the Chapter 10 compact file hgi_r5_slim.npz (ex10_2_qc_manhattan.py).
Output: figures 1 and 2, table_mr_bmi_instruments.csv, table_mr_bmi_methods.csv.   Needs numpy pandas scipy matplotlib and mr11.py.
"""
import os, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mr11 as mr
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
EXP, ASSOC, NPZ = "bmi_pulit.tbl.gz", "bmi_assoc.csv", "../ch10/hgi_r5_slim.npz"

# ---------------------------------------------------------------- 1. candidate instruments: variants reported by the BMI study with P < 5e-8
A = pd.read_csv(ASSOC); A = A[(A["P-VALUE"] < 5e-8) & A.SNPS.astype(str).str.fullmatch(r"rs\d+")]; cand = set(A.SNPS); print(f"{len(cand)} distinct reported BMI variants with P < 5e-8")

# ---------------------------------------------------------------- 2. their effects on the exposure, read from the full summary statistics
if not os.path.exists("bmi_exposure_snps.csv"):
    t0, parts = time.time(), []
    for ch in pd.read_csv(EXP, sep="\t", chunksize=2_000_000, usecols=["SNP", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"]):
        ch["rsid"] = ch.SNP.str.split(":").str[0]; parts.append(ch[ch.rsid.isin(cand)]); print(f"scanned a chunk ({time.time() - t0:.0f} s)", flush=True)
    pd.concat(parts).to_csv("bmi_exposure_snps.csv", index=False)
X = pd.read_csv("bmi_exposure_snps.csv"); X = X.rename(columns={"Allele1": "ea", "Allele2": "oa", "Freq1": "eaf_x", "Effect": "bx", "StdErr": "sx", "P-value": "px"}); X[["ea", "oa"]] = X[["ea", "oa"]].apply(lambda c: c.str.upper()); X = X.drop_duplicates("rsid")
print(f"{len(X)} candidate variants found in the exposure summary statistics")

# ---------------------------------------------------------------- 3. their effects on the outcome
if not os.path.exists("hgi_for_bmi.csv"):
    z = np.load(NPZ); rs = z["rsid"]; m = np.isin(rs, list(cand))
    pd.DataFrame({"rsid": rs[m], "oa1": z["a1"][m], "oa2": z["a2"][m], "by": z["beta"][m], "sy": z["se"][m], "eaf_y": z["eaf"][m], "ncoh": z["ncoh"][m]}).to_csv("hgi_for_bmi.csv", index=False)
Y = pd.read_csv("hgi_for_bmi.csv"); Y = Y.drop_duplicates("rsid", keep=False); print(f"{len(Y)} candidate variants found in the outcome summary statistics")

# ---------------------------------------------------------------- 4. harmonise alleles, drop ambiguous palindromes, check instrument strength
D = X.merge(Y, on="rsid"); keep, sign = mr.harmonise(D.ea.values, D.oa.values, D.eaf_x.values, D.oa1.values, D.oa2.values, D.eaf_y.values)
n0 = len(D); D = D[keep].copy(); D["by"] = D["by"].values * sign[keep]; D = D[np.isfinite(D.by) & (D.sy > 0) & (D.sx > 0)]
D["F"] = mr.fstat(D.bx.values, D.sx.values); print(f"{n0} matched; {len(D)} instruments after harmonisation and removal of ambiguous palindromes; F mean {D.F.mean():.0f}, min {D.F.min():.0f}; I2_GX = {mr.i2gx(D.bx.values, D.sx.values):.3f}")
D.to_csv("table_mr_bmi_instruments.csv", index=False)

# ---------------------------------------------------------------- 5. estimators
bx, by, sy = D.bx.values, D.by.values, D.sy.values; bw, sw = mr.wald(bx, by, sy)
I = mr.ivw(bx, by, sy); E = mr.egger(bx, by, sy); W = mr.weighted_median(bw, sw); n = len(D)
pf = lambda b, se, df=None: 2 * (mr.stats.norm.sf(abs(b / se)) if df is None else mr.stats.t.sf(abs(b / se), df))
M = pd.DataFrame([dict(method="IVW, fixed effects", beta=I["beta"], se=I["se_fixed"], p=pf(I["beta"], I["se_fixed"])), dict(method="IVW, random effects", beta=I["beta"], se=I["se_random"], p=pf(I["beta"], I["se_random"], n - 1)),
                  dict(method="MR-Egger", beta=E["slope"], se=E["se_slope"], p=E["p_slope"]), dict(method="Weighted median", beta=W["beta"], se=W["se"], p=pf(W["beta"], W["se"]))])
M["OR"] = np.exp(M.beta); M["OR_lo"] = np.exp(M.beta - 1.96 * M.se); M["OR_hi"] = np.exp(M.beta + 1.96 * M.se); M.to_csv("table_mr_bmi_methods.csv", index=False)
print(M.round(4).to_string(index=False)); print(f"Cochran Q = {I['Q']:.1f} on {I['df']} df (P = {I['Q_p']:.2e}); I2 = {I['I2']:.3f}; Egger intercept {E['intercept']:.4f} (SE {E['se_intercept']:.4f}, P = {E['p_intercept']:.3f})")
print(f"variants with a nominal outcome P < 0.05: {(2 * mr.stats.norm.sf(np.abs(by / sy)) < 0.05).sum()} of {n}; same direction as the IVW estimate: {(np.sign(bw) == np.sign(I['beta'])).mean():.2f}")
loo = mr.leave_one_out(bx, by, sy); infl = loo - I["beta"]; o = np.argsort(-np.abs(infl)); print("most influential variants (change in IVW beta when removed):", [(D.rsid.values[i], round(float(infl[i]), 4)) for i in o[:4]]); print(f"leave-one-out IVW range: {loo.min():.3f} to {loo.max():.3f}")
# the strongest outcome signals: are any BMI instruments genome-wide significant for COVID-19?
sig = D[2 * mr.stats.norm.sf(np.abs(D.by / D.sy)) < 5e-8]; print("instruments with outcome P < 5e-8:", sig[["rsid", "bx", "by", "sy"]].round(4).to_dict("records"))

# ---------------------------------------------------------------- Figure 1: scatter and funnel
s = np.sign(bx); xs, ys = bx * s, by * s; fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.7)); xx = np.linspace(0, xs.max() * 1.05, 50)
ax[0].errorbar(xs, ys, yerr=1.96 * sy, fmt="o", ms=3, color="#718096", ecolor="#e2e8f0", elinewidth=0.8, zorder=1)
ax[0].plot(xx, I["beta"] * xx, color="#2b6cb0", lw=2, label="IVW"); ax[0].plot(xx, E["intercept"] + E["slope"] * xx, color="#c53030", lw=1.6, ls="--", label="MR-Egger"); ax[0].plot(xx, W["beta"] * xx, color="#2f855a", lw=1.6, ls=":", label="Weighted median")
ax[0].axhline(0, color="black", lw=0.6); ax[0].set(xlabel="Effect on BMI (per SD, oriented to increase BMI)", ylabel="Effect on hospitalised COVID-19 (log odds ratio)", title=f"{n} BMI variants"); ax[0].legend(frameon=False, fontsize=9)
prec = np.abs(bx) / sy; ax[1].scatter(bw, prec, s=9, color="#718096"); [ax[1].axvline(v, color=c, lw=1.4, ls=l) for v, c, l in [(I["beta"], "#2b6cb0", "-"), (E["slope"], "#c53030", "--"), (W["beta"], "#2f855a", ":")]]
ax[1].axvline(0, color="black", lw=0.6); ax[1].set(xlabel="Wald ratio (log odds ratio per SD BMI)", ylabel="Instrument precision, |bx| / se(by)", title="Funnel plot", xlim=(-3, 3)); fig.tight_layout(); fig.savefig("figs11/fig1_bmi_scatter_funnel.png", dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 2: methods and leave-one-out
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), gridspec_kw={"width_ratios": [1, 1.15]}); y = np.arange(len(M))[::-1]
ax[0].errorbar(M.OR, y, xerr=[M.OR - M.OR_lo, M.OR_hi - M.OR], fmt="o", color="#2b6cb0", capsize=3); ax[0].axvline(1, color="black", ls="--", lw=0.8); ax[0].set_yticks(y); ax[0].set_yticklabels(M.method); ax[0].set_xscale("log")
ax[0].set_xticks([0.8, 1, 1.5, 2, 3]); ax[0].get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter()); ax[0].minorticks_off(); ax[0].set(xlabel="Odds ratio of hospitalised COVID-19 per SD of BMI (95% CI)", title="Four estimators")
top = o[:12][::-1]; ax[1].barh(range(len(top)), infl[top], color=np.where(infl[top] > 0, "#dd6b20", "#2b6cb0")); ax[1].set_yticks(range(len(top))); ax[1].set_yticklabels(D.rsid.values[top], fontsize=8.5); ax[1].axvline(0, color="black", lw=0.8)
ax[1].set(xlabel="Change in the IVW estimate when the variant is removed (log odds ratio)", title="Leave-one-out: the 12 most influential variants"); fig.tight_layout(); fig.savefig("figs11/fig2_bmi_methods_loo.png", dpi=200); plt.close(fig); print("figures written")
