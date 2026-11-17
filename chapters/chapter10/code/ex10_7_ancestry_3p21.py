"""
Chapter 10, Hands-on 7 - why a signal may not replicate: the 3p21.31 variants in five ancestries.
Session 10 asks: a GWAS signal in one population does not replicate in another; is it different linkage disequilibrium, different allele frequencies, or a false positive?
For the lead variant rs35081325 and the other 10 variants of the 95% credible set of Hands-on 4 (GRCh38), the script reads the phased genotypes of the 2,504 people of the 1000 Genomes high-coverage release
for chr3:45,798,429-45,898,429 (50 kb either side of the lead), and reports for each continental group (i) the frequency of the risk allele and the share of carriers, (ii) r2 with the lead variant,
(iii) the number of variants in the window in strong LD with the lead, and (iv) the number of cases (with equal controls) needed to replicate the odds ratio of 1.63 with 80% power at P < 0.05 and at P < 5e-8.
Needs bcftools and an internet connection for the first run (the region file is then kept: chr3_region_all.gt.gz); table_credible_set.csv from Hands-on 4; the sample panel of Hands-on 5.
Output: table_3p21_by_group.csv, table_3p21_r2_credible.csv, table_3p21_replication.csv, Figure 8.
"""
import os, subprocess, io
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
D = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"; VCF = f"{D}/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.filtered.shapeit2-duohmm-phased.vcf.gz"
LEAD, LEAD_POS, HALF = "rs35081325", 45_848_429, 50_000; LO, HI = LEAD_POS - HALF, LEAD_POS + HALF; GROUPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
panel = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", sep="\t").set_index("sample"); ids = panel.index.tolist()
if not os.path.exists("chr3_region_all.gt.gz"):
    open("all_samples.txt", "w").write("\n".join(ids)); head = subprocess.run(f"bcftools view -h -r chr3:{LO}-{HI} -S all_samples.txt --force-samples {VCF} | tail -1", shell=True, check=True, capture_output=True, text=True).stdout.strip().split("\t")[9:]
    open("chr3_region_samples.txt", "w").write("\n".join(head))
    subprocess.run(f"bcftools view -r chr3:{LO}-{HI} -S all_samples.txt --force-samples -m2 -M2 -v snps {VCF} | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' | gzip > chr3_region_all.gt.gz", shell=True, check=True)
samples = open("chr3_region_samples.txt").read().split(); G = pd.read_csv("chr3_region_all.gt.gz", sep="\t", header=None, dtype=str); info = G.iloc[:, :3].copy(); info.columns = ["pos", "ref", "alt"]; info["pos"] = info.pos.astype(int)
alt = G.iloc[:, 3:].apply(lambda c: c.str.count("1")).values.astype(np.int8)            # copies of the ALT allele: the phased genotypes are 0|0, 0|1, 1|0 or 1|1
assert alt.shape[1] == len(samples); pop = panel.super_pop.reindex(samples).values; print(f"{alt.shape[1]:,} people, {len(info):,} biallelic SNPs in chr3:{LO:,}-{HI:,}")
cs = pd.read_csv("table_credible_set.csv"); cs = cs.merge(info.assign(row=range(len(info))), on="pos"); cs = cs[((cs.a1 == cs.alt) & (cs.a2 == cs.ref)) | ((cs.a1 == cs.ref) & (cs.a2 == cs.alt))].reset_index(drop=True)
X = alt[cs.row.values].astype(float); X = np.where((cs.a1 == cs.alt).values[:, None], X, 2 - X)       # dosage of the effect allele a1
lead = cs.index[cs.rsid == LEAD][0]; allX = alt.astype(float); lead_row = int(cs.row[lead]); flip_lead = cs.a1[lead] != cs.alt[lead]; lead_dos = (2 - alt[lead_row]) if flip_lead else alt[lead_row]
rows, r2rows, rep = [], [], []
def r2(a, b):
    if a.std() == 0 or b.std() == 0: return np.nan
    return float(np.corrcoef(a, b)[0, 1] ** 2)
ncrit = {a: stats.chi2.isf(a, 1) for a in (0.05, 5e-8)}
def n_needed(p0, orr, a, target=0.8):
    lo, hi = 10.0, 1e10
    crit = ncrit[a]
    for _ in range(100):
        mid = np.sqrt(lo * hi); p1 = orr * p0 / (1 - p0 + orr * p0); pb = (p0 + p1) / 2; nc = (p1 - p0) ** 2 / (pb * (1 - pb) * (1 / (2 * mid) + 1 / (2 * mid)))
        if stats.ncx2.sf(crit, 1, nc) < target: lo = mid
        else: hi = mid
    return hi
for g in GROUPS:
    m = pop == g; n = int(m.sum()); d = lead_dos[m]; f = float(d.mean() / 2); carriers = float((d >= 1).mean())
    rr = np.array([r2(lead_dos[m].astype(float), (2 - alt[i][m]) if False else allX[i][m]) for i in range(len(info))]); dosg = allX[:, m]
    rr = np.array([r2(lead_dos[m].astype(float), dosg[i]) for i in range(len(info))]); strong = int(np.nansum(rr >= 0.8)) - (1 if f > 0 else 0); mod = int(np.nansum(rr >= 0.5)) - (1 if f > 0 else 0)
    rows.append(dict(group=g, people=n, risk_allele_freq=f, carriers=carriers, window_variants_polymorphic=int(np.sum(dosg.std(axis=1) > 0)), others_r2_ge_0_8=max(strong, 0), others_r2_ge_0_5=max(mod, 0)))
    for i, rs in cs.iterrows(): r2rows.append(dict(group=g, rsid=rs.rsid, risk_allele_freq=float(X[i][m].mean() / 2), r2_with_lead=r2(lead_dos[m].astype(float), X[i][m])))
    if f > 0: rep.append(dict(group=g, risk_allele_freq=f, cases_p_0_05=n_needed(f, 1.63, 0.05), cases_p_5e8=n_needed(f, 1.63, 5e-8)))
    else: rep.append(dict(group=g, risk_allele_freq=f, cases_p_0_05=np.nan, cases_p_5e8=np.nan))
T = pd.DataFrame(rows); R = pd.DataFrame(r2rows); P = pd.DataFrame(rep); T.to_csv("table_3p21_by_group.csv", index=False); R.to_csv("table_3p21_r2_credible.csv", index=False); P.to_csv("table_3p21_replication.csv", index=False)
print(T.round(3).to_string(index=False)); print(R.pivot(index="rsid", columns="group", values="r2_with_lead").round(2).to_string()); print(P.round(3).to_string(index=False))
# ---------------------------------------------------------------- figure: frequency and r2 by group
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.2)); col = {"AFR": "#2f855a", "AMR": "#805ad5", "EAS": "#dd6b20", "EUR": "#2b6cb0", "SAS": "#c53030"}
ax[0].bar(range(5), 100 * T.carriers, color=[col[g] for g in T.group]); ax[0].set_xticks(range(5)); ax[0].set_xticklabels(T.group)
for i, (c, f) in enumerate(zip(T.carriers, T.risk_allele_freq)): ax[0].text(i, 100 * c + 1.2, f"{100 * c:.1f}%\nfreq {f:.3f}", ha="center", fontsize=8.5)
ax[0].set(ylabel="Carriers of the rs35081325 risk allele (%)", title=f"Carriage ranges from {100 * T.carriers.min():.1f}% to {100 * T.carriers.max():.1f}%", ylim=(0, max(100 * T.carriers) * 1.35))
piv = R.pivot(index="rsid", columns="group", values="r2_with_lead").loc[cs.rsid, GROUPS]; im = ax[1].imshow(piv.values, cmap="Blues", vmin=0, vmax=1, aspect="auto"); ax[1].set_xticks(range(5)); ax[1].set_xticklabels(GROUPS); ax[1].set_yticks(range(len(piv))); ax[1].set_yticklabels(piv.index, fontsize=8)
for i in range(piv.shape[0]):
    for j in range(5):
        v = piv.values[i, j]; ax[1].text(j, i, "–" if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=7.5, color="white" if (not np.isnan(v) and v > 0.55) else "black")
ax[1].set_title("r² of the credible-set variants with the lead", fontsize=10.5); [ax[1].spines[s].set_visible(False) for s in ax[1].spines]
fig.tight_layout(); fig.savefig("figs10/fig8_3p21_ancestry.png", dpi=200); print("figure written")
