import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, hypergeom
from Bio import SeqIO
# ---------------- H2: independent SNP distance computation (byte arrays, different code path from the chapter)
recs = list(SeqIO.parse("ebola_aln.fasta", "fasta")); ids = [r.id for r in recs]; X = np.array([np.frombuffer(str(r.seq).upper().encode(), dtype=np.uint8) for r in recs]); n, L = X.shape
dt = pd.read_csv("ebola_dates.csv").set_index("acc").date; dates = pd.to_datetime(dt.loc[ids]); print(f"H2 genomes {n}, positions {L:,}, dates {dates.min().date()} to {dates.max().date()}")
ok = np.isin(X, np.frombuffer(b"ACGT", dtype=np.uint8)); D = np.zeros((n, n), dtype=np.int32)
for i in range(n): D[i] = (((X[i] != X) & ok[i]) & ok).sum(1)
iu = np.triu_indices(n, 1); d = D[iu]; days = np.abs((dates.values[:, None] - dates.values[None, :]).astype("timedelta64[D]").astype(int))[iu]
print(f"pairs {len(d):,}; median {np.median(d):.0f}; range {d.min()}-{d.max()}; within 2 SNPs {100 * (d <= 2).mean():.1f}%; identical {100 * (d == 0).mean():.1f}%; Spearman(distance, days) {spearmanr(d, days)[0]:.3f}")
for lo, hi in [(0, 7), (8, 30), (31, 60), (61, 90), (91, 150)]:
    m = (days >= lo) & (days <= hi); print(f"  days {lo}-{hi}: pairs {m.sum():,}, median {np.median(d[m]):.0f}, within2 {100 * (d[m] <= 2).mean():.1f}%, within5 {100 * (d[m] <= 5).mean():.1f}%")
y = days <= 30; r = rankdata(-d); auc = (r[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum()); auc = 1 - auc if auc < 0.5 else auc
# AUC for 'smaller distance => closer in time': probability that a close pair has smaller distance than a distant pair (ties half)
from scipy.stats import mannwhitneyu
U = mannwhitneyu(d[y], d[~y], alternative="less").statistic; auc2 = 1 - U / (y.sum() * (~y).sum()); print(f"AUC (sampled within 30 days) {auc2:.3f}")
for t in (0, 1, 2, 3, 5, 10):
    p = d <= t; tp = (p & y).sum(); print(f"  threshold {t}: called {p.sum():,}, sensitivity {100 * tp / y.sum():.1f}%, specificity {100 * ((~p) & (~y)).sum() / (~y).sum():.1f}%, PPV {100 * tp / p.sum():.1f}%")
rate = 1.11e-3; per_gen = rate * L; per_day = per_gen / 365; print(f"clock: {per_gen:.1f} SNPs per genome-year, {per_day:.4f} per day; days to ancestor for 1,2,5,10 SNPs: {[round(k / (2 * per_day), 1) for k in (1, 2, 5, 10)]}")
# ---------------- H4: lineage diversity from the monthly counts
c = pd.read_csv("lapis_lineage_monthly.csv"); c["country"] = c.country.replace({"United Kingdom": "UK", "USA": "USA"}); T = pd.read_csv("table_lineage_diversity.csv")
def simp(x): x = np.asarray(x, float); N = x.sum(); return 1 - (x * (x - 1)).sum() / (N * (N - 1))
rows = []
for (ct, mo), g in c.groupby(["country", "month"]):
    cnt = g.groupby("lineage")["count"].sum(); cnt = cnt[cnt > 0].values; N = cnt.sum()
    exp_rich = (1 - hypergeom(N, cnt, 500).pmf(0)).sum() if N >= 500 else np.nan                                   # exact expected richness in a draw of 500 without replacement
    rows.append(dict(country=ct, month=mo, N=N, lineages=len(cnt), D=simp(cnt), top=cnt.max() / N, exp_rich=exp_rich))
V = pd.DataFrame(rows); print("H4 countries in raw file:", V.country.unique().tolist(), "| table countries:", T.country.unique().tolist())
T["country"] = T.country.replace({"United Kingdom": "UK"}); M = V.merge(T, on=["country", "month"], suffixes=("_v", "_t")); print(f"months matched {len(M)}; max |D_obs recomputed - table| = {(M.D - M.D_observed).abs().max():.4f}; max |sequences| diff {(M.N - M.sequences).abs().max()}; max |lineages| diff {(M.lineages - M.lineages_v if 'lineages_v' in M else M.lineages_x - M.lineages_y).abs().max() if False else 'n/a'}")
M["Dr_minus_Do"] = (M.D_rarefied - M.D_observed); print(f"rarefied minus observed D: UK max |diff| {M[M.country == 'UK'].Dr_minus_Do.abs().max():.4f}, USA {M[M.country == 'USA'].Dr_minus_Do.abs().max():.4f}")
M["rich_ratio_exact"] = M.lineages_t / M.exp_rich if 'lineages_t' in M else M.lineages / M.exp_rich
uk = M[(M.country == "UK") & (M.N > 10000)]; print(f"UK months with more than 10,000 sequences: {len(uk)}; median ratio observed/expected richness (exact hypergeometric) {uk.rich_ratio_exact.median():.2f}; max {uk.rich_ratio_exact.max():.2f} in {uk.loc[uk.rich_ratio_exact.idxmax(), 'month']}")
for ct in ("UK", "USA"):
    s = M[(M.country == ct) & M.D_rarefied.notna()]; lo = s.loc[s.D_rarefied.idxmin()]; hi = s.loc[s.D_rarefied.idxmax()]; print(f"{ct}: months with D_rarefied {len(s)}; min {lo.D_rarefied:.3f} {lo.month} (top share {lo.top_share:.2f}, {lo.top_lineage}); max {hi.D_rarefied:.3f} {hi.month}")
print("USA 2022-02:", T[(T.country == "USA") & (T.month == "2022-02")][["sequences", "top_lineage", "top_share", "D_observed", "D_rarefied"]].to_dict("records")); print("UK 2022-03:", T[(T.country == "UK") & (T.month == "2022-03")][["sequences", "top_lineage", "top_share", "D_observed", "D_rarefied"]].to_dict("records"))
print("totals:", c.groupby("country")["count"].sum().to_dict())
