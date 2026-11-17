"""
Chapter 3, Hands-on 1 - Assembly quality control on public Salmonella genomes.
Data: assembly statistics (length, number of contigs, N50) for 911,210 isolates of the NCBI Pathogen Detection Salmonella table (salmonella_assembly_qc.csv.gz, built from the files of Chapters 1 and 2).
Three questions: how variable is genome length (and which genomes fall outside a robust range), how do three illustrative rule sets sort the same genomes, and does the failure rate depend on the sequencing laboratory?
The rule sets are illustrations, not programme standards; real programmes publish their own. Output: table_assembly_qc.csv, table_assembly_qc_labs.csv, Figure 1.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
from asm3 import load_asm
d = load_asm(); d = d.dropna(subset=["asm_stats_length_bp", "asm_stats_n_contig", "asm_stats_contig_n50"]).copy()
d["len_mb"] = d.asm_stats_length_bp / 1e6; d["contigs"] = d.asm_stats_n_contig; d["n50_kb"] = d.asm_stats_contig_n50 / 1e3
print(f"{len(d):,} isolates with assembly statistics; length median {d.len_mb.median():.2f} Mb (5th-95th percentile {d.len_mb.quantile(.05):.2f}-{d.len_mb.quantile(.95):.2f}); contigs median {d.contigs.median():.0f} (95th percentile {d.contigs.quantile(.95):.0f}); N50 median {d.n50_kb.median():.0f} kb (5th percentile {d.n50_kb.quantile(.05):.0f})")
# 1. robust range for genome length: median +/- 3 scaled MADs on the log scale
lg = np.log10(d.len_mb); med = lg.median(); mad = 1.4826 * np.median(np.abs(lg - med)); lo, hi = 10 ** (med - 3 * mad), 10 ** (med + 3 * mad); d["len_flag"] = (d.len_mb < lo) | (d.len_mb > hi)
print(f"robust genome-length range {lo:.2f}-{hi:.2f} Mb; {d.len_flag.mean():.2%} of isolates fall outside it: {int((d.len_mb < lo).sum()):,} short, {int((d.len_mb > hi).sum()):,} long; extremes {d.len_mb.min():.2f} and {d.len_mb.max():.2f} Mb")
# 2. three illustrative rule sets
rules = {"Lenient": dict(c=500, n=20), "Moderate": dict(c=300, n=50), "Strict": dict(c=200, n=100)}
for k, r in rules.items(): d[f"fail_{k}"] = (d.contigs > r["c"]) | (d.n50_kb < r["n"]) | d.len_flag
T = pd.DataFrame([dict(rule=k, max_contigs=r["c"], min_n50_kb=r["n"], fail_pct=d[f"fail_{k}"].mean() * 100, fail_n=int(d[f"fail_{k}"].sum()), fails_on_contigs_pct=(d.contigs > r["c"]).mean() * 100, fails_on_n50_pct=(d.n50_kb < r["n"]).mean() * 100) for k, r in rules.items()]); T.round(3).to_csv("table_assembly_qc.csv", index=False); print(T.round(2).to_string(index=False))
only_strict = d[d["fail_Strict"] & ~d["fail_Moderate"]]; print(f"isolates failing Strict but passing Moderate: {len(only_strict):,} ({len(only_strict) / len(d):.1%})")
# 3. by year and by laboratory (Moderate rule)
yr = d[d.year.between(2012, 2024)].groupby("year").agg(isolates=("contigs", "size"), fail_pct=("fail_Moderate", lambda s: s.mean() * 100), median_contigs=("contigs", "median"), median_n50=("n50_kb", "median")); print(yr.round(1).T.to_string())
labs = d.groupby("lab").agg(isolates=("contigs", "size"), fail_moderate=("fail_Moderate", "mean"), fail_strict=("fail_Strict", "mean"), median_contigs=("contigs", "median"), median_n50=("n50_kb", "median"), median_len=("len_mb", "median")).query("isolates >= 2000").sort_values("isolates", ascending=False)
labs.round(3).to_csv("table_assembly_qc_labs.csv"); print(f"{len(labs)} laboratories with at least 2,000 isolates; failure rate (Moderate) ranges {labs.fail_moderate.min():.1%} to {labs.fail_moderate.max():.1%}; interquartile range {labs.fail_moderate.quantile(.25):.1%}-{labs.fail_moderate.quantile(.75):.1%}")
print(labs.sort_values("fail_moderate", ascending=False).head(4)[["isolates", "fail_moderate", "median_contigs", "median_n50"]].round(3).to_string())
by_serovar = d[d.serovar.notna()].groupby(d.serovar.str.lower()).agg(n=("contigs", "size"), med_len=("len_mb", "median")).query("n >= 3000").sort_values("med_len"); print("median genome length by serovar (serovars with 3,000+ isolates): shortest", by_serovar.index[0], round(by_serovar.med_len.iloc[0], 2), "Mb; longest", by_serovar.index[-1], round(by_serovar.med_len.iloc[-1], 2), "Mb")
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [1, 1.1, 1.1]})
ax[0].hist(d.len_mb.clip(3.5, 6.5), bins=120, color="#2b6cb0"); ax[0].axvspan(0, lo, color="#c53030", alpha=0.15); ax[0].axvspan(hi, 7, color="#c53030", alpha=0.15); ax[0].set_xlim(3.5, 6.5); ax[0].set(xlabel="Assembly length (Mb)", ylabel="Isolates", title="Genome length and its robust range"); ax[0].set_yscale("log")
h = ax[1].hist2d(np.log10(d.contigs.clip(1, 3000)), np.log10(d.n50_kb.clip(1, 6000)), bins=70, cmap="Blues", norm=matplotlib.colors.LogNorm()); fig.colorbar(h[3], ax=ax[1], pad=0.01).set_label("Isolates", fontsize=8.5)
for k, c in [("Lenient", "#2f855a"), ("Moderate", "#dd6b20"), ("Strict", "#c53030")]: r = rules[k]; ax[1].plot([np.log10(r["c"])] * 2, [np.log10(r["n"]), 3.8], color=c, lw=1.5, label=k); ax[1].plot([0, np.log10(r["c"])], [np.log10(r["n"])] * 2, color=c, lw=1.5)
ax[1].set(xlabel="Contigs (log10)", ylabel="N50 (log10 kb)", title="Contigs and N50 with three rule sets"); ax[1].legend(frameon=False, fontsize=8.3, title="Pass region is below and right", title_fontsize=7.5, loc="upper right")
tl = labs.head(14)[::-1]; ax[2].barh(range(len(tl)), tl.fail_moderate * 100, color="#dd6b20", label="Moderate"); ax[2].barh(range(len(tl)), tl.fail_strict * 100, height=0.4, color="#c53030", label="Strict"); short = {"Public Health Agency of Canada": "PHAC (Canada)", "United States Department of Agriculture|Food Safety and Inspection Service": "USDA FSIS", "New York State Department of Health": "NY State DoH"}
ax[2].set_yticks(range(len(tl))); ax[2].set_yticklabels([short.get(s, s)[:24] for s in tl.index], fontsize=7.5); ax[2].set(xlabel="Isolates failing the rule (%)", title="Failure rate by sequencing laboratory"); ax[2].legend(frameon=False, fontsize=8.3, loc="lower right")
fig.tight_layout(); fig.savefig("figs3/fig1_assembly_qc.png", dpi=200); plt.close(fig); print("figure written")
