"""
Chapter 10, Hands-on 3 - Population stratification with a known truth.
Real genotypes: 1000 Genomes phase 3, chromosome 22, 2,504 people from five continental groups. A trait is simulated with 10 known causal SNPs plus an
environmental effect that follows ancestry. The association scan is run three ways in PLINK: no covariates, 2 principal components, 10 principal components.
Output: Figures 3 and 4, table_stratification.csv.  Needs bcftools and plink1.9 (apt-get install bcftools plink1.9), pandas, numpy, matplotlib.
"""
import os, subprocess, time, urllib.request
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
G = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"; VCF = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"; PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
sh = lambda c: subprocess.run(c, shell=True, check=True, capture_output=True, text=True).stdout
for f in (VCF, PANEL):
    if not os.path.exists(f): urllib.request.urlretrieve(f"{G}/{f}", f)
t0 = time.time()
if not os.path.exists("chr22.bed"):                                     # biallelic SNPs with minor-allele frequency of at least 5%
    if not os.path.exists("chr22_common.vcf.gz"): sh(f"bcftools view -m2 -M2 -v snps -q 0.05:minor {VCF} -Oz -o chr22_common.vcf.gz")     # about 90 s
    sh("plink1.9 --vcf chr22_common.vcf.gz --double-id --set-missing-var-ids @:#:\\$1:\\$2 --make-bed --out chr22 --silent")     # the VCF has no variant IDs, so build them from position and alleles
    print(f"filtered in {time.time() - t0:.0f} s")
if not os.path.exists("chr22.eigenvec"):
    sh("plink1.9 --bfile chr22 --indep-pairwise 50 5 0.2 --out prune --silent && plink1.9 --bfile chr22 --extract prune.prune.in --pca 10 --out chr22 --silent")
fam = pd.read_csv("chr22.fam", sep=" ", header=None, names=["fid", "iid", "p", "m", "sex", "ph"]); n_snp = sum(1 for _ in open("chr22.bim")); panel = pd.read_csv(PANEL, sep="\t").set_index("sample")
pcs = pd.read_csv("chr22.eigenvec", sep=" ", header=None, names=["fid", "iid"] + [f"PC{i}" for i in range(1, 11)]); pcs["pop"] = pcs.iid.map(panel.super_pop)
print(f"{len(fam):,} people, {n_snp:,} common SNPs, {sum(1 for _ in open('prune.prune.in')):,} after LD pruning; groups: {pcs['pop'].value_counts().to_dict()}")
eig = pd.read_csv("chr22.eigenval", header=None)[0]; print("variance share of PC1 to PC4 (of the 10 PCs):", (eig / eig.sum()).round(3).head(4).tolist())

# ---------------------------------------------------------------- simulate: 10 causal SNPs (2% of variance each) + an ancestry-driven environment (10% of variance) + noise
rng = np.random.default_rng(10); bim = pd.read_csv("chr22.bim", sep="\t", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
frq = pd.read_csv(sh("plink1.9 --bfile chr22 --freq --out frq --silent; cat frq.frq").splitlines() and "frq.frq", sep=r"\s+")
cand = frq[(frq.MAF > 0.15)].SNP.tolist(); causal = sorted(rng.choice(cand, 10, replace=False)); open("causal.txt", "w").write("\n".join(causal))
sh("plink1.9 --bfile chr22 --extract causal.txt --recode A --out causal --silent"); g = pd.read_csv("causal.raw", sep=r"\s+").set_index("IID").drop(columns=["FID", "PAT", "MAT", "SEX", "PHENOTYPE"])
gs = (g - g.mean()) / g.std(); genetic = gs.sum(axis=1) * np.sqrt(0.02) / np.sqrt(10) * np.sqrt(10) * 1.0               # 10 SNPs, each with a standardised effect of sqrt(0.02)
pc = pcs.set_index("iid")[["PC1", "PC2"]].reindex(g.index); env = (0.8 * pc.PC1 / pc.PC1.std() + 0.6 * pc.PC2 / pc.PC2.std()); env = env / env.std() * np.sqrt(0.10)
y = genetic + env + rng.normal(0, np.sqrt(1 - 0.20 - 0.10), len(g)); y = (y - y.mean()) / y.std()
print(f"simulated trait: genetic share {genetic.var() / y.var():.2f}, ancestry share {env.var() / y.var():.2f}")
pd.DataFrame({"FID": fam.fid, "IID": fam.iid, "y": y.reindex(fam.iid).values}).to_csv("pheno.txt", sep="\t", index=False)
pcs[["fid", "iid"] + [f"PC{i}" for i in range(1, 11)]].to_csv("covar10.txt", sep="\t", index=False, header=["FID", "IID"] + [f"PC{i}" for i in range(1, 11)])
pcs[["fid", "iid", "PC1", "PC2"]].to_csv("covar2.txt", sep="\t", index=False, header=["FID", "IID", "PC1", "PC2"])

# ---------------------------------------------------------------- three association scans
def scan(name, cov):
    c = f"--covar {cov} --covar-name " + ("PC1,PC2" if cov == "covar2.txt" else ",".join(f"PC{i}" for i in range(1, 11))) if cov else ""
    sh(f"plink1.9 --bfile chr22 --pheno pheno.txt --pheno-name y {c} --linear hide-covar --allow-no-sex --out scan_{name} --silent")
    d = pd.read_csv(f"scan_{name}.assoc.linear", sep=r"\s+").dropna(subset=["P"]); return d
res = {"No covariates": scan("none", None), "2 PCs": scan("pc2", "covar2.txt"), "10 PCs": scan("pc10", "covar10.txt")}
cpos = bim[bim.id.isin(causal)].pos.values; far = lambda d: ~np.any(np.abs(d.BP.values[:, None] - cpos[None, :]) < 1_000_000, axis=1)   # exclude 1 Mb around every causal SNP
rows = []
for k, d in res.items():
    lam = stats.chi2.isf(np.median(d.P), 1) / stats.chi2.ppf(0.5, 1); f = far(d); lam_far = stats.chi2.isf(np.median(d.P[f]), 1) / stats.chi2.ppf(0.5, 1)
    fp = int((d.P[f] < 5e-8).sum()); hit = int(sum((d[d.SNP == s].P < 5e-8).any() for s in causal)); rows.append(dict(model=k, lambda_all=round(lam, 2), lambda_away_from_causal=round(lam_far, 2), false_positives=fp, tested_away=int(f.sum()), causal_found=hit))
T = pd.DataFrame(rows); T.to_csv("table_stratification.csv", index=False); print(T.to_string(index=False)); print("causal SNPs:", [str(c) for c in causal])

# ---------------------------------------------------------------- Figures 3 and 4
fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4), gridspec_kw={"width_ratios": [1.5, 1]}); cm = {"AFR": "#dd6b20", "AMR": "#d69e2e", "EAS": "#2f855a", "EUR": "#2b6cb0", "SAS": "#805ad5"}
for s, d in pcs.groupby("pop"): ax[0].scatter(d.PC1, d.PC2, s=7, color=cm[s], label=s, alpha=0.8)
ax[0].set(xlabel="PC1", ylabel="PC2", title="Principal components of chromosome 22 genotypes"); ax[0].legend(frameon=False, markerscale=2, fontsize=8.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.15))
ax[1].bar(range(1, 11), 100 * eig / eig.sum(), color="#2b6cb0"); ax[1].set(xlabel="Principal component", ylabel="% of variance among the 10", title="Share of the ten components", xticks=range(1, 11))
fig.tight_layout(); fig.savefig("figs10/fig3_pca.png", dpi=200); plt.close(fig)
fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.3), sharey=True)
for a, (k, d), col in zip(ax, res.items(), ["#c53030", "#d69e2e", "#2b6cb0"]):
    p = np.sort(d.P.values); e = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p)); o = -np.log10(p); idx = np.unique(np.r_[np.arange(0, len(p), 20), np.arange(len(p) - 400, len(p))]); idx = idx[idx < len(p)]
    a.plot(e[idx], o[idx], ".", ms=3, color=col); a.plot([0, e.max()], [0, e.max()], color="black", lw=0.8); a.axhline(-np.log10(5e-8), color="#a0aec0", ls=":", lw=1)
    r = T[T.model == k].iloc[0]; a.set(title=f"{k}: lambda = {r.lambda_all:.2f}", xlabel="Expected -log10(P)"); a.text(0.04, 0.84, f"P < 5e-8 away from causal SNPs: {r.false_positives:,}\ncausal SNPs found: {r.causal_found} of 10", transform=a.transAxes, fontsize=9)
ax[0].set_ylabel("Observed -log10(P)"); fig.tight_layout(); fig.savefig("figs10/fig4_stratification_qq.png", dpi=200); plt.close(fig); print("figures written")
