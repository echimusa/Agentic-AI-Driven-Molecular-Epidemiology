"""
Chapter 10, Hands-on 4 - Fine-map the 3p21.31 locus with approximate Bayes factors (Wakefield) and European LD.
Summary statistics: COVID-19 HGI release 5, hospitalised vs population (GCST011081, GRCh38), from hgi_r5_slim.npz written by ex10_2.
LD: 1000 Genomes high-coverage phased genotypes (GRCh38), European samples, read for the region directly from the EBI server with bcftools.
Output: credible-set table, table_credible_set.csv, Figure 5.  Single causal variant assumption: stated in the chapter.
"""
import os, json, subprocess, urllib.request, io
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
CHR, LEAD_POS, HALF = 3, 45_848_429, 500_000; LO, HI = LEAD_POS - HALF, LEAD_POS + HALF; LEAD = "rs35081325"
D = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
VCF = f"{D}/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.filtered.shapeit2-duohmm-phased.vcf.gz"

# ---------------------------------------------------------------- summary statistics for the region
if os.path.exists("hgi_r5_slim.npz"):                               # written by ex10_2; the small region file shipped with the chapter is used otherwise
    z = np.load("hgi_r5_slim.npz"); z = {k: z[k] for k in ["rsid", "chrom", "pos", "a1", "a2", "beta", "se", "p"]}
else:
    r_ = pd.read_csv("hgi_chr3_region.tsv.gz", sep="\t", keep_default_na=False, na_values=[""]); z = {k: r_[k].values for k in ["rsid", "chrom", "pos", "a1", "a2", "beta", "se", "p"]}
m = (z["chrom"] == CHR) & (z["pos"] >= LO) & (z["pos"] <= HI) & np.isfinite(z["beta"]) & (z["se"] > 0)
S = pd.DataFrame({"rsid": z["rsid"][m], "pos": z["pos"][m], "a1": z["a1"][m], "a2": z["a2"][m], "beta": z["beta"][m], "se": z["se"][m], "p": z["p"][m]}); S = S[(S.a1.str.len() == 1) & (S.a2.str.len() == 1)].drop_duplicates("pos", keep=False).reset_index(drop=True)
print(f"{len(S):,} single-nucleotide variants with statistics in chr{CHR}:{LO:,}-{HI:,}; best P {S.p.min():.1e} at {S.loc[S.p.idxmin(), 'rsid']}")

# ---------------------------------------------------------------- European LD from 1000 Genomes (GRCh38)
if not os.path.exists("eur_samples.txt"):
    ped = pd.read_csv(f"{D}/20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+"); ped[ped.Superpopulation == "EUR"].SampleID.to_csv("eur_samples.txt", index=False, header=False)
if not os.path.exists("chr3_region_eur.gt.gz"):
    cmd = f"bcftools view -r chr{CHR}:{LO}-{HI} -S eur_samples.txt --force-samples -m2 -M2 -v snps {VCF} | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' | gzip -9 > chr3_region_eur.gt.gz"; subprocess.run(cmd, shell=True, check=True)
G = pd.read_csv("chr3_region_eur.gt.gz", sep="\t", header=None, dtype=str); n_eur = G.shape[1] - 3
dos = G.iloc[:, 3:].apply(lambda c: c.str.count("1")).values.astype(np.float32); info = G.iloc[:, :3]; info.columns = ["pos", "ref", "alt"]; info["pos"] = info.pos.astype(int)
print(f"{n_eur} European samples, {len(info):,} biallelic SNPs in the region")

# match alleles; dosage of the effect allele a1
J = S.merge(info.assign(row=range(len(info))), on="pos"); J = J[((J.a1 == J.alt) & (J.a2 == J.ref)) | ((J.a1 == J.ref) & (J.a2 == J.alt))].reset_index(drop=True)
X = dos[J.row.values]; X = np.where((J.a1 == J.alt).values[:, None], X, 2 - X); freq = X.mean(axis=1) / 2; J = J[(freq > 0.005) & (freq < 0.995)].reset_index(drop=True); X = X[((freq > 0.005) & (freq < 0.995))]
lead = J.index[J.rsid == LEAD][0]; r = np.array([np.corrcoef(X[lead], X[i])[0, 1] for i in range(len(J))]); J["r2"] = r ** 2
print(f"{len(J):,} variants with both statistics and LD; lead {LEAD} at {J.loc[lead, 'pos']:,}")

# ---------------------------------------------------------------- approximate Bayes factors and credible set
def abf(beta, se, sd):
    V, W = se ** 2, sd ** 2; r_ = W / (V + W); return 0.5 * np.log(1 - r_) + 0.5 * r_ * (beta / se) ** 2                  # log BF for association (Wakefield 2009)
def cred(logbf, level=0.95):
    pp = np.exp(logbf - logbf.max()); pp /= pp.sum(); o = np.argsort(-pp); k = int(np.searchsorted(np.cumsum(pp[o]), level) + 1); return pp, o[:k]
out = {}
for sd in (0.1, 0.15, 0.2):
    pp, idx = cred(abf(J.beta.values, J.se.values, sd)); out[sd] = (pp, idx); print(f"prior sd {sd}: 95% credible set of {len(idx)} variants; top variant {J.rsid[idx[0]]} with posterior probability {pp[idx[0]]:.3f}")
pp, idx = out[0.15]; J["pp"] = pp; C = J.loc[idx, ["rsid", "pos", "a1", "a2", "beta", "se", "p", "r2", "pp"]].assign(OR=lambda d: np.exp(d.beta)).sort_values("pp", ascending=False); C.to_csv("table_credible_set.csv", index=False)
print(C.head(10).round(4).to_string(index=False)); print(f"lead variant {LEAD}: P {J.loc[lead, 'p']:.1e}, posterior probability {J.pp[lead]:.3f}; credible-set members with r2 > 0.8 to the lead: {(C.r2 > 0.8).sum()} of {len(C)}")
for rs in ("rs35044562",):
    if (J.rsid == rs).any(): q = J[J.rsid == rs].iloc[0]; print(f"{rs}: P {q.p:.1e}, r2 with lead {q.r2:.2f}, posterior probability {q.pp:.3f}")
print(f"variants with P < 5e-8 in the region: {(J.p < 5e-8).sum():,}; with r2 > 0.8 to the lead: {(J.r2 > 0.8).sum():,}; sum of posterior probability in the credible set {C.pp.sum():.3f}")

# ---------------------------------------------------------------- genes in the region (Ensembl, GRCh38)
if not os.path.exists("region_genes.json"):
    u = f"https://rest.ensembl.org/overlap/region/human/{CHR}:{LO}-{HI}?feature=gene;content-type=application/json"; open("region_genes.json", "w").write(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "ch10"}), timeout=60).read().decode())
genes = [g for g in json.load(open("region_genes.json")) if g.get("biotype") == "protein_coding"]

# ---------------------------------------------------------------- Figure 5
fig, ax = plt.subplots(3, 1, figsize=(10.4, 8.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.4, 1.1]}); px = J.pos / 1e6
bins = [0, 0.2, 0.4, 0.6, 0.8, 1.01]; cols = ["#2b6cb0", "#63b3ed", "#68d391", "#f6ad55", "#e53e3e"]; cat = pd.cut(J.r2, bins, labels=False, include_lowest=True)
for k, c in enumerate(cols): s = cat == k; ax[0].scatter(px[s], -np.log10(J.p[s]), s=9, color=c, label=f"r2 {bins[k]:.1f} to {min(bins[k+1], 1):.1f}")
ax[0].scatter(px[idx], -np.log10(J.p[idx]), s=34, facecolors="none", edgecolors="black", linewidths=0.8, label="95% credible set"); ax[0].scatter(px[lead], -np.log10(J.p[lead]), s=70, marker="D", color="#6b21a8", zorder=5, label=LEAD)
ax[0].axhline(-np.log10(5e-8), color="#a0aec0", ls="--", lw=1); ax[0].set(ylabel="-log10(P)", title="3p21.31: hospitalised COVID-19, with LD to the lead variant (European samples)"); ax[0].legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
ax[1].vlines(px, 0, J.pp, color="#4a5568", lw=0.8); ax[1].scatter(px[idx], J.pp[idx], s=14, color="#c53030", zorder=3); ax[1].set(ylabel="Posterior probability")
for k_, g in enumerate(sorted(genes, key=lambda g: g["start"])):
    row = k_ % 2; a, b = max(g["start"], LO) / 1e6, min(g["end"], HI) / 1e6; ax[2].plot([a, b], [row, row], lw=5, color="#2b6cb0", solid_capstyle="butt", alpha=0.7); ax[2].text((a + b) / 2, row + 0.28, g.get("external_name", ""), ha="center", fontsize=7)
ax[2].set(ylim=(-0.6, 1.9), yticks=[], xlabel=f"Position on chromosome {CHR} (Mb, GRCh38)"); ax[2].spines["left"].set_visible(False); ax[2].set_xlim(LO / 1e6, HI / 1e6)
fig.tight_layout(); fig.savefig("figs10/fig5_finemap.png", dpi=200); plt.close(fig); print("figure written")
