"""
Chapter 10, Hands-on 5 - the quality-control command of Session 10 on real genotypes, Hardy-Weinberg equilibrium across ancestries, FST and linkage disequilibrium.
Real genotypes: 1000 Genomes phase 3, chromosome 22, 2,504 people from five continental groups (AFR, AMR, EAS, EUR, SAS). A random 100,000 of the biallelic SNPs (seed 10) are converted to PLINK format.
  Step 1  the QC command of the session:  plink --bfile data --geno 0.02 --mind 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out qc   (each filter alone and all together)
  Step 2  a known-fault test of --geno and --mind: missing genotypes are injected into a small copy, and the filters must find exactly those
  Step 3  Hardy-Weinberg exact tests in the pooled sample and inside each group (Wahlund effect)
  Step 4  Hudson FST between the five groups, from genotype counts, on all variants and on common variants; Weir-Cockerham FST from PLINK for one pair as a cross-check
  Step 5  linkage disequilibrium (mean r2) by distance in each group, on the same common variants
Needs bcftools, plink1.9, pandas, numpy, scipy.  Downloads the 197 MB chromosome 22 VCF and the sample panel when they are absent.
Output: table_qc_steps.csv, table_qc_injected.csv, table_hwe_pooled_within.csv, table_fst_pairs.csv, table_ld_by_distance.csv, Figure 6, chr22_group_geno_counts.tsv.gz (genotype counts per variant and group, used by the later exercises and the validator).
"""
import os, re, subprocess, time, urllib.request, itertools
import numpy as np, pandas as pd
G = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"; VCF = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"; PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
sh = lambda c: subprocess.run(c, shell=True, check=True, capture_output=True, text=True).stdout
for f in (VCF, PANEL):
    if not os.path.exists(f): urllib.request.urlretrieve(f"{G}/{f}", f)
GROUPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]; panel = pd.read_csv(PANEL, sep="\t").set_index("sample"); t0 = time.time()
if not os.path.exists("sub.bed"):
    pos = pd.read_csv(__import__("io").StringIO(sh(f"bcftools view -G -m2 -M2 -v snps {VCF} | bcftools query -f '%CHROM\\t%POS\\n'")), sep="\t", header=None, names=["chr", "pos"]).drop_duplicates()
    pick = np.sort(np.random.default_rng(10).choice(len(pos), 100_000, replace=False)); pos.iloc[pick].to_csv("sample_pos.tsv", sep="\t", header=False, index=False)
    sh(f"bcftools view -m2 -M2 -v snps -T sample_pos.tsv {VCF} -Oz -o sub100k.vcf.gz")
    sh("plink1.9 --vcf sub100k.vcf.gz --double-id --set-missing-var-ids @:#:\\$1:\\$2 --make-bed --out sub --silent")
    print(f"{len(pos):,} biallelic SNPs on chromosome 22; 100,000 sampled; converted in {time.time() - t0:.0f} s")
n_var = sum(1 for _ in open("sub.bim")); n_ind = sum(1 for _ in open("sub.fam")); print(f"{n_var:,} variants, {n_ind:,} people")
def run(args, out):
    log = sh(f"plink1.9 --bfile sub {args} --out {out} --silent; cat {out}.log"); return log
def removed(log, key):
    m = re.search(rf"(\d+) variants removed due to {key}", log); return int(m.group(1)) if m else 0
def people_removed(log):
    m = re.search(r"(\d+) people removed due to missing genotype data", log); return int(m.group(1)) if m else 0
# ---------------------------------------------------------------- step 1: the QC command of the session, one filter at a time and all together
steps = []
for name, args in [("--geno 0.02", "--geno 0.02 --make-bed"), ("--mind 0.02", "--mind 0.02 --make-bed"), ("--maf 0.01", "--maf 0.01 --make-bed"), ("--hwe 1e-6 (pooled, all 2,504 people)", "--hwe 1e-6 --make-bed"),
                   ("all four together", "--geno 0.02 --mind 0.02 --maf 0.01 --hwe 1e-6 --make-bed")]:
    log = run(args, "qc_tmp"); nv = sum(1 for _ in open("qc_tmp.bim")); ni = sum(1 for _ in open("qc_tmp.fam"))
    steps.append(dict(filter=name, variants_left=nv, variants_removed=n_var - nv, people_left=ni, people_removed=n_ind - ni))
S = pd.DataFrame(steps); S.to_csv("table_qc_steps.csv", index=False); print(S.to_string(index=False))
# ---------------------------------------------------------------- step 2: known-fault test of --geno and --mind
bim = pd.read_csv("sub.bim", sep="\t", header=None); keep_v = bim.iloc[:2000, 1].tolist(); fam = pd.read_csv("sub.fam", sep=" ", header=None)
rng = np.random.default_rng(11); keep_i = sorted(rng.choice(len(fam), 500, replace=False)); open("small_v.txt", "w").write("\n".join(keep_v)); open("small_i.txt", "w").write("\n".join(f"{fam.iloc[i, 0]} {fam.iloc[i, 1]}" for i in keep_i))
sh("plink1.9 --bfile sub --extract small_v.txt --keep small_i.txt --recode --out small --silent")
ped = pd.read_csv("small.ped", sep=" ", header=None, dtype=str); bad_i = sorted(rng.choice(500, 10, replace=False)); bad_v = sorted(rng.choice(2000, 40, replace=False))
G_ = ped.iloc[:, 6:].values.reshape(500, 2000, 2).copy()
for i in bad_i: G_[i, rng.random(2000) < 0.15, :] = "0"              # 10 people with about 15% missing genotypes
for v in bad_v: G_[rng.random(500) < 0.05, v, :] = "0"             # 40 variants with about 5% missing genotypes
ped.iloc[:, 6:] = G_.reshape(500, 4000); ped.to_csv("small_faulty.ped", sep=" ", header=False, index=False); os.replace("small.map", "small_faulty.map")
log = sh("plink1.9 --file small_faulty --geno 0.02 --mind 0.02 --make-bed --out small_qc --silent; cat small_qc.log")
inj = pd.DataFrame([dict(injected_people=len(bad_i), found_people=people_removed(log), injected_variants=len(bad_v), found_variants=removed(log, "missing genotype data"),
                         variants_left=sum(1 for _ in open("small_qc.bim")), people_left=sum(1 for _ in open("small_qc.fam")))])
inj.to_csv("table_qc_injected.csv", index=False); print(inj.to_string(index=False))
# ---------------------------------------------------------------- step 3: Hardy-Weinberg, pooled and inside each group
def hardy(keep=None, tag="all"):
    if keep is not None: open(f"keep_{tag}.txt", "w").write("\n".join(f"{i} {i}" for i in keep))
    sh(f"plink1.9 --bfile sub --keep-allele-order {'--keep keep_' + tag + '.txt' if keep is not None else ''} --hardy --out hw_{tag} --silent")      # without --keep-allele-order PLINK makes A1 the minor allele of each subset, so groups would count different alleles
    h = pd.read_csv(f"hw_{tag}.hwe", sep=r"\s+")[["SNP", "A1", "A2", "GENO", "P"]]; g = h.GENO.str.split("/", expand=True).astype(int); g.columns = ["n11", "n12", "n22"]; return pd.concat([h[["SNP", "A1", "A2", "P"]], g], axis=1)
hp = hardy(); grp = {g: hardy(panel.index[panel.super_pop == g].tolist(), g) for g in GROUPS}
pos = bim.set_index(1)[3]
cnt = hp[["SNP", "A1", "A2"]].copy(); cnt["pos"] = cnt.SNP.map(pos)
for g in GROUPS:
    for c in ("n11", "n12", "n22"): cnt[f"{g}_{c}"] = grp[g][c].values
cnt.to_csv("chr22_group_geno_counts.tsv.gz", sep="\t", index=False)
thr = 1e-6; maf_pool = (2 * hp.n22 + hp.n12) / (2 * (hp.n11 + hp.n12 + hp.n22)); maf_pool = np.minimum(maf_pool, 1 - maf_pool)
ok = maf_pool >= 0.01; fail_pool = (hp.P < thr) & ok; fail_in = {g: (grp[g].P < thr) for g in GROUPS}; any_in = np.logical_or.reduce([fail_in[g] for g in GROUPS])
H = pd.DataFrame([dict(test="pooled sample", variants_tested_maf_ge_0_01=int(ok.sum()), fail_p_1e6=int(fail_pool.sum())), *[dict(test=f"inside {g} only", variants_tested_maf_ge_0_01=int(ok.sum()), fail_p_1e6=int((fail_in[g] & ok).sum())) for g in GROUPS],
                  dict(test="fail pooled but pass in all five groups", variants_tested_maf_ge_0_01=int(ok.sum()), fail_p_1e6=int((fail_pool & ~any_in).sum())), dict(test="fail in at least one group but pass pooled", variants_tested_maf_ge_0_01=int(ok.sum()), fail_p_1e6=int((ok & ~fail_pool & any_in).sum()))])
H.to_csv("table_hwe_pooled_within.csv", index=False); print(H.to_string(index=False))
# ---------------------------------------------------------------- step 4: Hudson FST between groups from counts (ratio of averages), on all 100,000 variants and on the variants with MAF >= 0.05; Weir-Cockerham FST from PLINK for EUR vs EAS
sh("plink1.9 --bfile sub --freq --out frq_all --silent"); fq = pd.read_csv("frq_all.frq", sep=r"\s+"); common = fq[fq.MAF >= 0.05].SNP.tolist(); open("common.txt", "w").write("\n".join(common)); is_common = hp.SNP.isin(set(common)).values
def freq_n(g): n = grp[g].n11 + grp[g].n12 + grp[g].n22; return (2 * grp[g].n11 + grp[g].n12) / (2 * n), 2 * n
def hudson(a, b, mask):
    p1, n1 = freq_n(a); p2, n2 = freq_n(b); num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1); den = p1 * (1 - p2) + p2 * (1 - p1); good = (den > 0) & mask
    return float(num[good].sum() / den[good].sum())
rows = [dict(group_a=a, group_b=b, hudson_fst_all=hudson(a, b, np.ones(len(hp), bool)), hudson_fst_common=hudson(a, b, is_common)) for a, b in itertools.combinations(GROUPS, 2)]
F = pd.DataFrame(rows); sel = panel.super_pop.isin(["EUR", "EAS"]); within = pd.DataFrame({"fid": panel.index[sel], "iid": panel.index[sel], "g": panel.super_pop[sel].values})
within.to_csv("within_eur_eas.txt", sep=" ", index=False, header=False); out = sh("plink1.9 --bfile sub --extract common.txt --within within_eur_eas.txt --fst --out fst_ee --silent; cat fst_ee.log"); m = re.search(r"Mean Fst estimate: ([\d.eE+-]+)", out); wc = re.search(r"Weighted Fst estimate: ([\d.eE+-]+)", out)
F["plink_wc_mean_common"] = np.nan; F["plink_wc_weighted_common"] = np.nan; i = F[(F.group_a == "EAS") & (F.group_b == "EUR")].index[0]; F.loc[i, "plink_wc_mean_common"] = float(m.group(1)); F.loc[i, "plink_wc_weighted_common"] = float(wc.group(1))
F.to_csv("table_fst_pairs.csv", index=False); print(F.round(4).to_string(index=False))
# ---------------------------------------------------------------- step 5: mean r2 by distance in each group, on the variants with MAF >= 0.05 in the pooled sample
bins = [(0, 10_000), (10_000, 25_000), (25_000, 50_000), (50_000, 100_000), (100_000, 250_000)]; rows = []
for g in GROUPS:
    sh(f"plink1.9 --bfile sub --extract common.txt --keep keep_{g}.txt --r2 --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0 --out ld_{g} --silent")
    d = pd.read_csv(f"ld_{g}.ld", sep=r"\s+", usecols=["BP_A", "BP_B", "R2"]); d["dist"] = (d.BP_B - d.BP_A).abs()
    for lo, hi in bins:
        x = d[(d.dist > lo) & (d.dist <= hi)] if lo > 0 else d[d.dist <= hi]; rows.append(dict(group=g, bin_kb=f"{lo // 1000}-{hi // 1000}", pairs=len(x), mean_r2=float(x.R2.mean())))
L = pd.DataFrame(rows); L.to_csv("table_ld_by_distance.csv", index=False); print(len(common), "common variants;"); print(L.pivot(index="bin_kb", columns="group", values="mean_r2").round(3).to_string())
print(f"done in {time.time() - t0:.0f} s")

# ---------------------------------------------------------------- Figure 6
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
col = {"AFR": "#2f855a", "AMR": "#805ad5", "EAS": "#dd6b20", "EUR": "#2b6cb0", "SAS": "#c53030"}
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.2), gridspec_kw={"width_ratios": [1, 1.05, 1.1]})
hh = H.iloc[:6]; names = ["pooled"] + GROUPS; ax[0].bar(range(6), hh.fail_p_1e6, color=["#4a5568"] + [col[g] for g in GROUPS]); ax[0].set_xticks(range(6)); ax[0].set_xticklabels(names)
for i, v in enumerate(hh.fail_p_1e6): ax[0].text(i, v + 60, f"{int(v):,}", ha="center", fontsize=8.5)
ax[0].set(ylabel="Variants failing HWE at P < 1e-6 (of 16,312)", title="Pooling ancestries breaks HWE"); mat = np.zeros((5, 5)); 
for r in F.itertuples(): i, j = GROUPS.index(r.group_a), GROUPS.index(r.group_b); mat[i, j] = mat[j, i] = r.hudson_fst_common
im = ax[1].imshow(mat, cmap="Blues", vmin=0, vmax=0.2); ax[1].set_xticks(range(5)); ax[1].set_xticklabels(GROUPS); ax[1].set_yticks(range(5)); ax[1].set_yticklabels(GROUPS); [ax[1].spines[s_].set_visible(False) for s_ in ax[1].spines]
for i in range(5):
    for j in range(5): ax[1].text(j, i, "" if i == j else f"{mat[i, j]:.3f}", ha="center", va="center", fontsize=9, color="white" if mat[i, j] > 0.11 else "black")
ax[1].set_title("Hudson FST between groups (common variants)", fontsize=10.5)
for g in GROUPS:
    d = L[L.group == g]; ax[2].plot(range(len(d)), d.mean_r2, "o-", color=col[g], label=g)
ax[2].set_xticks(range(5)); ax[2].set_xticklabels(L[L.group == "AFR"].bin_kb.tolist()); ax[2].set(xlabel="Distance between variants (kb)", ylabel="Mean r²", title="Linkage disequilibrium is lowest in AFR"); ax[2].legend(frameon=False, ncol=2, fontsize=9)
fig.tight_layout(); fig.savefig("figs10/fig6_qc_structure.png", dpi=200); print("figure written")
