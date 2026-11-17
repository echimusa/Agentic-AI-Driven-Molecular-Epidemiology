"""
Chapter 11, Hands-on 4 - Drug-target Mendelian randomisation at IL6R with correlated variants.
Exposure: C-reactive protein (Said et al., GWAS Catalog GCST90029070, harmonised GRCh38). Outcome: hospitalised COVID-19, HGI release 5 (GCST011081).
Instruments: CRP-associated variants (P < 5e-8) in a 270 kb window around IL6R. Variants are correlated, so the estimate uses a generalised-least-squares IVW with a
European LD matrix from the 1000 Genomes high-coverage data (GRCh38), read for the window directly from the EBI server with bcftools.
Output: table_cis_mr.csv, Figure 4.
"""
import os, subprocess
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
CHR, LO, HI = 1, 154_300_000, 154_570_000; NPZ = "../ch10/hgi_r5_slim.npz"
D = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"; VCF = f"{D}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"

# ---------------------------------------------------------------- exposure and outcome statistics for the window
if not os.path.exists("crp_il6r_region.csv"):
    parts = []
    for ch in pd.read_csv("crp_said.h.tsv.gz", sep="\t", chunksize=2_000_000, usecols=["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "p_value"], dtype={"hm_chrom": str}):
        parts.append(ch[(ch.hm_chrom == str(CHR)) & ch.hm_pos.between(LO, HI)])
    pd.concat(parts).to_csv("crp_il6r_region.csv", index=False)
X = pd.read_csv("crp_il6r_region.csv").rename(columns={"hm_rsid": "rsid", "hm_pos": "pos", "hm_effect_allele": "xa1", "hm_other_allele": "xa2", "hm_beta": "bx", "standard_error": "sx", "p_value": "px"}).drop_duplicates("pos", keep=False)
if not os.path.exists("hgi_il6r_region.csv"):
    z = np.load(NPZ); m = (z["chrom"] == CHR) & (z["pos"] >= LO) & (z["pos"] <= HI)
    pd.DataFrame({"pos": z["pos"][m], "ya1": z["a1"][m], "ya2": z["a2"][m], "by": z["beta"][m], "sy": z["se"][m], "py": z["p"][m], "ncoh": z["ncoh"][m]}).to_csv("hgi_il6r_region.csv", index=False)
Y = pd.read_csv("hgi_il6r_region.csv").drop_duplicates("pos", keep=False); print(f"window chr{CHR}:{LO:,}-{HI:,}: {len(X):,} variants with CRP statistics, {len(Y):,} with COVID-19 statistics")

# ---------------------------------------------------------------- European LD from 1000 Genomes (GRCh38)
if not os.path.exists("chr1_il6r_eur.gt"):
    subprocess.run(f"bcftools view -r chr{CHR}:{LO}-{HI} -S eur_samples.txt --force-samples -m2 -M2 -v snps {VCF} | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' > chr1_il6r_eur.gt", shell=True, check=True)
G = pd.read_csv("chr1_il6r_eur.gt", sep="\t", header=None, dtype=str); dos = G.iloc[:, 3:].apply(lambda c: c.str.count("1")).values.astype(np.float32); L = G.iloc[:, :3]; L.columns = ["pos", "ref", "alt"]; L["pos"] = L.pos.astype(int); L["row"] = range(len(L))
print(f"{dos.shape[1]} European samples, {len(L):,} biallelic SNPs in the window")

# ---------------------------------------------------------------- merge, align every dataset to the reference-panel ALT allele (all files are forward-strand GRCh38)
J = X.merge(Y, on="pos").merge(L.drop_duplicates("pos", keep=False), on="pos")
def sgn(a1, a2, ref, alt): return np.where((a1 == alt) & (a2 == ref), 1, np.where((a1 == ref) & (a2 == alt), -1, 0))
J["sx_sign"] = sgn(J.xa1, J.xa2, J.ref, J.alt); J["sy_sign"] = sgn(J.ya1, J.ya2, J.ref, J.alt); J = J[(J.sx_sign != 0) & (J.sy_sign != 0)].reset_index(drop=True)
J["bx"] *= J.sx_sign; J["by"] *= J.sy_sign; J["maf"] = dos[J.row.values].mean(axis=1) / 2; J = J[(J.maf > 0.01) & (J.maf < 0.99)].reset_index(drop=True); J["F"] = (J.bx / J.sx) ** 2
print(f"{len(J):,} variants present in CRP, COVID-19 and the LD panel with matching alleles; {(J.F > 29.72).sum()} with CRP P < 5e-8")
J["lp_x"] = -(stats.norm.logsf(np.abs(J.bx / J.sx)) + np.log(2)) / np.log(10)                      # -log10 P from the Z score, which does not underflow
J["lp_y"] = -np.log10(J.py.clip(lower=1e-300)); R = np.corrcoef(dos[J.row.values]); lead = int(J.F.idxmax())
print(f"lead CRP variant (largest Z): position {J.pos[lead]:,.0f} {J.rsid[lead]} (CRP -log10 P {J.lp_x[lead]:.0f}, COVID-19 P {J.py[lead]:.2f})")

# ---------------------------------------------------------------- estimators
def gls_ivw(bx, by, sy, rho):
    Om = np.outer(sy, sy) * rho; Oi = np.linalg.pinv(Om); a = bx @ Oi @ bx; b = (bx @ Oi @ by) / a; res = by - b * bx; Q = res @ Oi @ res; df = len(bx) - 1
    se_f = 1 / np.sqrt(a); return dict(beta=b, se_fixed=se_f, se_random=se_f * np.sqrt(max(1, Q / df)) if df > 0 else se_f, Q=Q, df=df)
def prune(r2max):
    idx = list(J[J.F > 29.72].sort_values("F", ascending=False).index); keep = []                  # F > 29.72 is the same as P < 5e-8; strongest first
    for i in idx:
        if all(R[i, j] ** 2 < r2max for j in keep): keep.append(i)
    return sorted(keep)
rows = []; wr = J.by[lead] / J.bx[lead]; rows.append(dict(analysis="Lead variant only (Wald ratio)", variants=1, beta=wr, se=abs(J.sy[lead] / J.bx[lead])))
if (J.rsid == "rs2228145").any():                                   # a well-studied IL6R missense variant, used on its own
    q = J[J.rsid == "rs2228145"].iloc[0]; rows.append(dict(analysis="rs2228145 only (Wald ratio)", variants=1, beta=q.by / q.bx, se=abs(q.sy / q.bx)))
for r2max in (0.2, 0.4, 0.6):
    k = prune(r2max); bx, by, sy = J.bx.values[k], J.by.values[k], J.sy.values[k]; rho = R[np.ix_(k, k)]
    nv = ivw_b = np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2); se_naive = 1 / np.sqrt(np.sum(bx ** 2 / sy ** 2))
    g = gls_ivw(bx, by, sy, rho)
    rows.append(dict(analysis=f"IVW ignoring LD, r2 < {r2max}", variants=len(k), beta=ivw_b, se=se_naive)); rows.append(dict(analysis=f"Correlated IVW (GLS), r2 < {r2max}", variants=len(k), beta=g["beta"], se=g["se_random"], Q=g["Q"], df=g["df"]))
T = pd.DataFrame(rows); T["OR"] = np.exp(T.beta); T["OR_lo"] = np.exp(T.beta - 1.96 * T.se); T["OR_hi"] = np.exp(T.beta + 1.96 * T.se); T["p"] = 2 * stats.norm.sf(np.abs(T.beta / T.se)); T.to_csv("table_cis_mr.csv", index=False); print(T.round(4).to_string(index=False))
k = prune(0.4); print("variants in the r2 < 0.4 instrument set:", [(J.rsid[i], int(J.pos[i]), round(float(J.bx[i]), 3), round(float(J.F[i]), 0)) for i in k][:12])

# ---------------------------------------------------------------- Figure 4
fig, ax = plt.subplots(2, 1, figsize=(10.2, 7.4), sharex=False, gridspec_kw={"height_ratios": [1.6, 1]}); px = J.pos / 1e6
bins = [0, 0.2, 0.4, 0.6, 0.8, 1.01]; cols = ["#2b6cb0", "#63b3ed", "#68d391", "#f6ad55", "#e53e3e"]; cat = pd.cut(R[lead] ** 2, bins, labels=False, include_lowest=True)
for a, col, ttl in [(ax[0], "px", "C-reactive protein"), ]: 
    for c_, cc in enumerate(cols): s = cat == c_; a.scatter(px[s], J.lp_x[s], s=9, color=cc, label=f"r2 {bins[c_]:.1f} to {min(bins[c_ + 1], 1):.1f} with the lead")
a2 = ax[0].twinx(); a2.scatter(px, -np.log10(J.py), s=7, marker="x", color="#4a5568", alpha=0.6, label="Hospitalised COVID-19"); a2.set_ylabel("Hospitalised COVID-19, -log10(P) (grey crosses)"); a2.spines["top"].set_visible(False)
ax[0].axhline(-np.log10(5e-8), color="#a0aec0", ls="--", lw=1); ax[0].set(xlabel=f"Position on chromosome {CHR} (Mb, GRCh38)", ylabel="C-reactive protein, -log10(P)", title="IL6R window: CRP signal (coloured by LD with the lead) and COVID-19 signal"); ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
kk = prune(0.4); bx, by, sy = J.bx.values[kk], J.by.values[kk], J.sy.values[kk]; g = gls_ivw(bx, by, sy, R[np.ix_(kk, kk)]); sg = np.sign(bx); bxo, byo = bx * sg, by * sg; xs = np.linspace(0, bxo.max() * 1.1, 30)   # oriented so that every variant raises CRP
ax[1].errorbar(bxo, byo, yerr=1.96 * sy, fmt="o", color="#4a5568", ecolor="#cbd5e0", capsize=2); ax[1].plot(xs, g["beta"] * xs, color="#2b6cb0", lw=2, label=f"Correlated IVW: OR {np.exp(g['beta']):.2f}"); ax[1].axhline(0, color="black", lw=0.6)
ax[1].set(xlabel="Effect of the CRP-raising allele on C-reactive protein", ylabel="Effect on hospitalised COVID-19 (log OR)", title=f"{len(kk)} variants after pruning to r2 < 0.4"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig("figs11/fig4_il6r_cis_mr.png", dpi=200); plt.close(fig); print("figure written")
