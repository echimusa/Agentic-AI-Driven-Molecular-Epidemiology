"""
Chapter 3, Hands-on 3 - Depth, systematic error and mixed samples: how read data turn into SNP calls and SNP distances.
A real 18,961-position genome (an Ebola virus genome of Chapter 13, alignment gaps removed) is the reference A. Isolate B differs from A at 3 known positions. Reads from B are simulated as per-position
counts: depth is negative binomial (mean d, dispersion 4, as GC and mapping bias make real depth uneven); random errors occur at rate e; a fraction h of positions are error hotspots (homopolymers, damage,
mapping artefacts) where 25% of reads show a wrong base; and a fraction c of reads come from a contaminating genotype that differs from A at 50 other positions. A position is called if depth >= 10 and the
non-reference allele fraction is at least f. Observed SNP distance = alt calls; false SNPs = alt calls at positions where B equals A.
Output: table_snp_depth.csv, table_snp_calling.csv, table_snp_threshold.csv, Figure 3.
"""
import numpy as np, pandas as pd
from Bio import SeqIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(3); seq = [str(r.seq).upper() for r in SeqIO.parse("ebola_aln.fasta", "fasta")]
ref = min(seq, key=lambda s: sum(c not in "ACGT" for c in s)); ref = "".join(c for c in ref if c in "ACGT"); L = len(ref); print(f"reference genome: {L:,} unambiguous bases")
KTRUE, MCONT, DMIN, NB_R, HOT_RATE = 3, 50, 10, 4.0, 0.25

def replicate(d, e=0.005, h=0.0, c=0.0, f=0.8, ktrue=KTRUE):
    Y = rng.negative_binomial(NB_R, NB_R / (NB_R + d), L)                                         # reads covering each position
    true = rng.choice(L, ktrue, replace=False) if ktrue else np.array([], int); alt = np.zeros(L, int)
    Yc = rng.binomial(Y, c)                                                                         # reads from the contaminating genotype
    err = rng.binomial(Y, e / 3.0)                                                                  # reads showing one specific wrong base through random error
    alt += err
    if h > 0:
        hot = rng.choice(L, int(round(h * L)), replace=False); alt[hot] = rng.binomial(Y[hot], HOT_RATE)     # hotspots: 25% of reads carry a wrong base
    if c > 0:
        cont = rng.choice(L, MCONT, replace=False); alt[cont] = np.maximum(alt[cont], Yc[cont])             # contaminant alleles at 50 positions
    if ktrue: alt[true] = Y[true] - Yc[true] - rng.binomial(Y[true] - Yc[true], e)                   # B's allele; contaminant reads (matching A) dilute it
    ok = Y >= DMIN; called = ok & (alt >= f * np.maximum(Y, 1)) & (alt > 0)
    tp = int(called[true].sum()) if ktrue else 0; fp = int(called.sum() - tp)
    return dict(callable=ok.mean(), tp=tp, fp=fp, observed=int(called.sum()), true_callable=int(ok[true].sum()) if ktrue else 0)

def run(reps=400, **kw):
    r = pd.DataFrame([replicate(**kw) for _ in range(reps)]); return dict(callable_pct=100 * r.callable.mean(), sens_pct=100 * r.tp.mean() / max(kw.get("ktrue", KTRUE), 1), false_snps=r.fp.mean(), observed=r.observed.mean(), p_obs_le2=(r.observed <= 2).mean() * 100, p_exact=(r.observed == kw.get("ktrue", KTRUE)).mean() * 100)

# 1. depth
rows = [dict(depth=d, **run(d=d)) for d in (5, 10, 20, 30, 50, 100)]; A = pd.DataFrame(rows); A.round(2).to_csv("table_snp_depth.csv", index=False); print("depth (error 0.5%, no hotspots, no contamination, f = 0.8):"); print(A.round(2).to_string(index=False))
# 2. allele-fraction threshold with hotspots and contamination
rows = []
for h, c in [(0.0, 0.0), (0.005, 0.0), (0.0, 0.10), (0.0, 0.30), (0.005, 0.30)]:
    for f in (0.2, 0.4, 0.5, 0.7, 0.9):
        rows.append(dict(hotspots=h, contamination=c, f_min=f, **run(d=50, h=h, c=c, f=f)))
B = pd.DataFrame(rows); B.round(2).assign(hotspots=B.hotspots).to_csv("table_snp_calling.csv", index=False)       # keep the hotspot share exact: round(2) turned 0.005 into 0.0 and made the rows ambiguous; print("allele-fraction threshold (depth 50):"); print(B.round(2).to_string(index=False))
# 3. what depth and threshold do to a 2-SNP rule when the true distance is 3
rows = []
for d in (10, 20, 50):
    for c in (0.0, 0.1, 0.3):
        for f in (0.5, 0.9): rows.append(dict(depth=d, contamination=c, f_min=f, **run(d=d, c=c, f=f, h=0.005)))
C = pd.DataFrame(rows); C.round(2).to_csv("table_snp_threshold.csv", index=False); print("true distance 3, hotspots 0.5%: probability the observed distance is at most 2 (a '2 SNP' rule would call the pair linked)"); print(C[["depth", "contamination", "f_min", "observed", "false_snps", "p_obs_le2", "p_exact"]].round(1).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.6)); x = A.depth
ax[0].plot(x, A.callable_pct, "o-", color="#2b6cb0", label="Genome callable (depth >= 10)"); ax[0].plot(x, A.sens_pct, "o-", color="#c53030", label="True SNPs called"); ax[0].set_xscale("log"); ax[0].set_xticks(list(x)); ax[0].set_xticklabels([str(i) for i in x]); ax[0].set(xlabel="Mean depth", ylabel="Percent", title="Depth: what can be called", ylim=(0, 105)); ax[0].legend(frameon=False, fontsize=8.5, loc="lower right")
cols = {(0.0, 0.0): "#2b6cb0", (0.005, 0.0): "#dd6b20", (0.0, 0.30): "#2f855a", (0.005, 0.30): "#c53030"}
for (h, c), col in cols.items():
    s = B[(B.hotspots == h) & (B.contamination == c)]; ax[1].plot(s.f_min, s.false_snps, "o-", color=col, label=f"hotspots {h * 100:g}%, contamination {c * 100:g}%")
ax[1].set_yscale("symlog", linthresh=0.1); ax[1].set(xlabel="Minimum alternative allele fraction to call a SNP", ylabel="False SNPs per genome", title="Permissive thresholds create SNPs"); ax[1].legend(frameon=False, fontsize=7.6, loc="upper right")
for c, col in [(0.0, "#2b6cb0"), (0.10, "#2f855a"), (0.30, "#c53030")]:
    s = B[(B.hotspots == 0.0) & (B.contamination == c)] if c else B[(B.hotspots == 0.0) & (B.contamination == 0.0)]; ax[2].plot(s.f_min, s.sens_pct, "o-", color=col, label=f"contamination {c * 100:g}%")
ax[2].set(xlabel="Minimum alternative allele fraction to call a SNP", ylabel="True SNPs called (%)", title="Strict thresholds lose real SNPs in mixed samples", ylim=(0, 105)); ax[2].legend(frameon=False, fontsize=8.5, loc="lower left")
fig.tight_layout(); fig.savefig("figs3/fig3_snp_calling.png", dpi=200); plt.close(fig); print("figure written")
