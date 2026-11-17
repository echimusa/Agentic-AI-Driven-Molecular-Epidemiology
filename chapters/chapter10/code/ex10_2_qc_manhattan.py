"""
Chapter 10, Hands-on 2 - Quality control and loci from public summary statistics.
Input: GWAS Catalog GCST011081 (COVID-19 Host Genetics Initiative release 5, hospitalised vs population, GRCh38, harmonised file).
Output: QC numbers, Figure 2 (QQ and Manhattan), table_loci.csv, sumstats_p1e-4.tsv.gz (for Hands-on 4), qc_leads_for_r.csv (for the R cross-check).
"""
import os, sys, json, time, urllib.request
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011081/harmonised/32404885-GCST011081-MONDO_0100096.h.tsv.gz")
F, NPZ = "hgi_r5_hosp.h.tsv.gz", "hgi_r5_slim.npz"; N_CASES, N_CONTROLS = 9986, 1877672             # sample sizes listed for GCST011081 in the GWAS Catalog
if not os.path.exists(F): urllib.request.urlretrieve(URL, F)

# ---------------------------------------------------------------- read once into compact arrays
if not os.path.exists(NPZ):
    cols = ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency", "n"]
    parts, t0 = [], time.time()
    for ch in pd.read_csv(F, sep="\t", usecols=cols, chunksize=2_000_000, dtype={"hm_chrom": str, "hm_other_allele": str, "hm_effect_allele": str}, na_values=["NA"]):
        parts.append(ch); print(f"read {sum(len(p) for p in parts):,} rows ({time.time() - t0:.0f} s)", flush=True)
    d = pd.concat(parts, ignore_index=True); del parts
    chrom = d.hm_chrom.map(lambda c: 0 if not isinstance(c, str) else 23 if c == "X" else 24 if c == "Y" else 25 if c in ("MT", "M") else int(c) if c.isdigit() else 0).astype(np.int8)
    np.savez(NPZ, rsid=np.asarray(d.hm_rsid.fillna("").to_numpy(), dtype="U14"), chrom=chrom.values, pos=d.hm_pos.fillna(0).to_numpy().astype(np.int64), a1=np.asarray(d.hm_effect_allele.fillna("").to_numpy(), dtype="U8"), a2=np.asarray(d.hm_other_allele.fillna("").to_numpy(), dtype="U8"),
             beta=d.hm_beta.values.astype(np.float64), se=d.standard_error.values.astype(np.float64), p=d.p_value.values.astype(np.float64), eaf=d.hm_effect_allele_frequency.values.astype(np.float64), ncoh=d.n.values.astype(np.float64)); del d
z = np.load(NPZ); rsid, chrom, pos, a1, a2, beta, se, p, eaf, ncoh = (z[k] for k in ["rsid", "chrom", "pos", "a1", "a2", "beta", "se", "p", "eaf", "ncoh"])
n = len(p); print(f"\n{n:,} variants")

# ---------------------------------------------------------------- QC checks that every public file should pass before use
ok = np.isfinite(beta) & np.isfinite(se) & np.isfinite(p) & (se > 0) & (p > 0) & (p <= 1)
print(f"valid beta, SE and P: {ok.mean():.4f}; P = 0 or missing: {(~np.isfinite(p) | (p <= 0)).sum():,}; effect-allele frequency missing: {np.isnan(eaf).mean():.3f}; chromosomes present: {sorted(set(chrom.tolist()))[:3]} to {max(chrom)}")
dup = pd.Series(chrom.astype(np.int64) * 10**10 + pos).duplicated().sum(); print(f"duplicated positions: {dup:,}")
zs = beta[ok] / se[ok]; p_wald = 2 * stats.norm.sf(np.abs(zs)); m = ok.copy(); m[ok] = p_wald > 1e-300
agree = np.abs(np.log10(p[m]) - np.log10(2 * stats.norm.sf(np.abs(beta[m] / se[m])))) < 0.05; print(f"reported P agrees with the Wald P from beta/SE (within 0.05 log10 units): {agree.mean():.4f}")
print(f"effect-allele frequency below 1% or above 99%: {((eaf < 0.01) | (eaf > 0.99))[~np.isnan(eaf)].mean():.3f}; contributing cohorts per variant: median {np.nanmedian(ncoh):.0f}, max {np.nanmax(ncoh):.0f}")

# ---------------------------------------------------------------- genomic inflation
lam = stats.chi2.isf(np.median(p[ok]), 1) / stats.chi2.ppf(0.5, 1)             # the median of chi-square is the chi-square of the median P
lam1000 = 1 + (lam - 1) * (1 / N_CASES + 1 / N_CONTROLS) / (1 / 1000 + 1 / 1000)
reg3 = (chrom == 3) & (pos > 45_000_000) & (pos < 47_000_000); lam_ex = stats.chi2.isf(np.median(p[ok & ~reg3]), 1) / stats.chi2.ppf(0.5, 1)
print(f"lambda_GC = {lam:.3f}; lambda_1000 = {lam1000:.3f}; lambda_GC excluding 3p21.31 = {lam_ex:.3f}")

# ---------------------------------------------------------------- independent loci: greedy distance clumping, 500 kb each side
sig = np.where(ok & (p < 5e-8))[0]; order = sig[np.argsort(p[sig])]; taken = np.zeros(len(order), bool); leads = []
for i, v in enumerate(order):
    if taken[i]: continue
    near = (chrom[order] == chrom[v]) & (np.abs(pos[order] - pos[v]) <= 500_000); taken |= near; leads.append((v, int(near.sum())))
L = pd.DataFrame({"chr": [int(chrom[v]) for v, _ in leads], "pos_GRCh38": [int(pos[v]) for v, _ in leads], "rsid": [rsid[v] for v, _ in leads], "effect_allele": [a1[v] for v, _ in leads], "other_allele": [a2[v] for v, _ in leads],
                  "OR": [round(float(np.exp(beta[v])), 3) for v, _ in leads], "eaf": [round(float(eaf[v]), 3) for v, _ in leads], "beta": [float(beta[v]) for v, _ in leads], "se": [float(se[v]) for v, _ in leads], "P": [float(p[v]) for v, _ in leads], "n_significant": [k for _, k in leads]})
# nearest protein-coding gene from Ensembl (GRCh38), 250 kb either side
def nearest(c, ps):
    ch = "X" if c == 23 else str(c); u = f"https://rest.ensembl.org/overlap/region/human/{ch}:{max(ps - 250000, 1)}-{ps + 250000}?feature=gene;content-type=application/json"
    for _ in range(3):
        try:
            g = [x for x in json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "ch10"}), timeout=30)) if x.get("biotype") == "protein_coding"]
            if not g: return ""
            g.sort(key=lambda x: 0 if x["start"] <= ps <= x["end"] else min(abs(x["start"] - ps), abs(x["end"] - ps))); return ", ".join(x.get("external_name", x["id"]) for x in g[:2])
        except Exception: time.sleep(2)
    return "lookup failed"
L["nearest_genes"] = pd.read_csv("table_loci.csv").nearest_genes.fillna("").tolist() if os.path.exists("table_loci.csv") and len(pd.read_csv("table_loci.csv")) == len(L) else [nearest(c, q) for c, q in zip(L.chr, L.pos_GRCh38)]; L.to_csv("table_loci.csv", index=False); print(f"\n{len(L)} independent loci with P < 5e-8:"); print(L[["chr", "pos_GRCh38", "rsid", "effect_allele", "OR", "eaf", "P", "n_significant", "nearest_genes"]].to_string(index=False))

# ---------------------------------------------------------------- files for later steps
region = ok & (p < 1e-4); pd.DataFrame({"rsid": rsid[region], "chr": chrom[region], "pos": pos[region], "a1": a1[region], "a2": a2[region], "beta": beta[region], "se": se[region], "p": p[region], "eaf": eaf[region]}).to_csv("sumstats_p1e-4.tsv.gz", sep="\t", index=False)
pd.DataFrame({"rsid": L.rsid, "beta": L.beta, "se": L.se, "p_reported": L.P}).to_csv("qc_leads_for_r.csv", index=False)

# ---------------------------------------------------------------- Figure 2: QQ and Manhattan
rng = np.random.default_rng(2); pp = p[ok]; obs = -np.log10(np.sort(pp)); exp_ = -np.log10((np.arange(1, len(pp) + 1) - 0.5) / len(pp))
keep = np.r_[np.arange(0, 20000), np.unique(np.geomspace(20000, len(pp) - 1, 6000).astype(int))]; ci = len(pp)
fig, ax = plt.subplots(1, 2, figsize=(13.2, 4.6), gridspec_kw={"width_ratios": [1, 2.6]})
ax[0].plot(exp_[keep], obs[keep], ".", ms=3, color="#2b6cb0"); ax[0].plot([0, exp_.max()], [0, exp_.max()], color="black", lw=0.8)
ax[0].set(xlabel="Expected -log10(P)", ylabel="Observed -log10(P)", title=f"QQ plot, lambda = {lam:.3f}"); ax[0].text(0.05, 0.9, f"lambda_1000 = {lam1000:.3f}", transform=ax[0].transAxes, fontsize=9)
sel = ok & ((p < 0.01) | (rng.random(n) < 0.02)); x = np.zeros(n); off = {}; cum = 0
for c in range(1, 24):
    if (chrom == c).any(): off[c] = cum; cum += pos[chrom == c].max() + 2_000_000
for c in off: x[chrom == c] = pos[chrom == c] + off[c]
y = -np.log10(p[sel]); cols = np.where(chrom[sel] % 2 == 0, "#2b6cb0", "#a0aec0")
ax[1].scatter(x[sel], y, s=2.5, c=cols, rasterized=True); ax[1].axhline(-np.log10(5e-8), color="#c53030", ls="--", lw=1)
ax[1].set_xticks([off[c] + pos[chrom == c].max() / 2 for c in off]); ax[1].set_xticklabels([("X" if c == 23 else str(c)) for c in off], fontsize=8, rotation=90); ax[1].set(xlabel="Chromosome", ylabel="-log10(P)", title=f"Hospitalised COVID-19 vs population: {len(L)} loci with P < 5e-8")
for r in L.itertuples():
    if r.P < 1e-9 or r.rsid in set(L.sort_values("P").rsid.head(12)): ax[1].annotate(r.nearest_genes.split(",")[0] if r.nearest_genes else r.rsid, (off[r.chr] + r.pos_GRCh38, -np.log10(r.P)), xytext=(3, 5), textcoords="offset points", fontsize=7.5)
fig.tight_layout(); fig.savefig("figs10/fig2_qc_manhattan.png", dpi=200); plt.close(fig); print("figure written")
