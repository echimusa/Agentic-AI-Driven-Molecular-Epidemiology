"""
Chapter 8, Hands-on 4 - the population-genetic view of the same genomes: diversity, differentiation, isolation by distance and Tajima's D.

Reads the alignment of the 591 Zika genomes used in Hands-on 1 (zika_aln.fasta, minus the nine clock outliers) and the place of each genome (phylo8.load_states).
For each region and each country with enough genomes it computes: nucleotide diversity (mean pairwise p-distance with pairwise deletion of missing sites),
Tajima's D on the genomes with at least 99% of the well-covered columns called and the columns without missing calls among them, Dxy (mean pairwise distance between two groups) and Hudson's Fst = 1 - Hw/Hb with a permutation test,
and a Mantel test between Fst/(1 - Fst) and the great-circle distance between country centroids for the American countries.
Writes table_popgen_groups.csv, table_popgen_regions.csv, table_popgen_countries.csv, table_popgen_ibd.csv and figs8/fig8_popgen.png. About a minute.
"""
import os, json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from Bio import SeqIO
import phylo8
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs8", exist_ok=True)
MIN_COUNTRY = 12; PERMS = 999; MANTEL_PERMS = 9999; rng = np.random.default_rng(8)
# ---------------------------------------------------------------- data
st = phylo8.load_states(); out = set(pd.read_csv("zika_clock_outliers.csv").acc); recs = [r for r in SeqIO.parse("zika_aln.fasta", "fasta") if r.id in st.index and r.id not in out]
ids = [r.id for r in recs]; n = len(ids); L = len(recs[0].seq); lut = np.full(256, 255, np.uint8)
for k, ch in enumerate("ACGT"): lut[ord(ch)] = k
A = np.vstack([lut[np.frombuffer(str(r.seq).upper().encode(), np.uint8)] for r in recs]); meta = st.loc[ids]; region = np.array([phylo8.region_of(c) for c in meta.cty]); cty = meta.cty.to_numpy()
print(f"{n} genomes, {L} alignment columns, {(A == 255).mean() * 100:.1f}% of cells missing or ambiguous")
# ---------------------------------------------------------------- pairwise distances with pairwise deletion
valid = A != 255; diff = np.zeros((n, n), np.int32); nval = np.zeros((n, n), np.int32)
for i in range(n):
    V = valid & valid[i]; diff[i] = ((A != A[i]) & V).sum(axis=1); nval[i] = V.sum(axis=1)
P = diff / np.maximum(nval, 1); np.fill_diagonal(P, np.nan)
colfrac = valid.mean(axis=0); core = colfrac >= 0.90; good = valid[:, core].mean(axis=1) >= 0.99           # Tajima's D uses the genomes with at least 99% of the 'core' columns (at least 90% called overall) called
print(f"core columns {int(core.sum())}; genomes with at least 99% of the core called {int(good.sum())} of {n}")
def within(idx):
    m = P[np.ix_(idx, idx)]; return float(np.nanmean(m))
def between(i1, i2): return float(P[np.ix_(i1, i2)].mean())
def tajima_d(idx):
    sub = A[idx]; cols = (sub != 255).all(axis=0); sub = sub[:, cols]; m = len(idx)
    if m < 4 or sub.shape[1] == 0: return np.nan, np.nan, int(cols.sum())
    counts = np.stack([(sub == k).sum(axis=0) for k in range(4)]); seg = (counts > 0).sum(axis=0) > 1; S = int(seg.sum())
    pairs_diff = (m * m - (counts.astype(float) ** 2).sum(axis=0)) / 2; pi = pairs_diff.sum() / (m * (m - 1) / 2)
    if S == 0: return np.nan, 0, int(cols.sum())
    a1 = sum(1 / i for i in range(1, m)); a2 = sum(1 / i ** 2 for i in range(1, m)); b1 = (m + 1) / (3 * (m - 1)); b2 = 2 * (m * m + m + 3) / (9 * m * (m - 1))
    c1 = b1 - 1 / a1; c2 = b2 - (m + 2) / (a1 * m) + a2 / a1 ** 2; e1 = c1 / a1; e2 = c2 / (a1 ** 2 + a2)
    return float((pi - S / a1) / np.sqrt(e1 * S + e2 * S * (S - 1))), S, int(cols.sum())
def fst_hudson(i1, i2, labels=None):
    hw = (within(i1) + within(i2)) / 2; hb = between(i1, i2); return 1 - hw / hb, hb
def perm_p(i1, i2, obs):
    pool = np.concatenate([i1, i2]); cnt = 0
    for _ in range(PERMS):
        rng.shuffle(pool); f, _ = fst_hudson(pool[:len(i1)], pool[len(i1):]); cnt += f >= obs
    return (cnt + 1) / (PERMS + 1)
# ---------------------------------------------------------------- groups: regions and countries of the Americas with enough genomes
regs = ["Americas", "Asia-Pacific", "Europe/Africa"]; ridx = {r: np.where(region == r)[0] for r in regs}
am_c = [c for c in pd.Series(cty[region == "Americas"]).value_counts().index if (cty == c).sum() >= MIN_COUNTRY]; cidx = {c: np.where(cty == c)[0] for c in am_c}
rows = []
for name, idx, level in [(r, ridx[r], "region") for r in regs] + [(c, cidx[c], "country (Americas)") for c in am_c]:
    idx_d = idx[good[idx]]; d, S, sites = tajima_d(idx_d); rows.append(dict(group=name, level=level, n=len(idx), pi_per_kb=1000 * within(idx), n_tajima=len(idx_d), tajima_d=d, segregating_sites=S, sites_used=sites))
G = pd.DataFrame(rows); G.round(4).to_csv("table_popgen_groups.csv", index=False); print(G.round(3).to_string(index=False))
R = []
for i in range(3):
    for j in range(i + 1, 3):
        f, hb = fst_hudson(ridx[regs[i]], ridx[regs[j]]); R.append(dict(group1=regs[i], group2=regs[j], dxy_per_kb=1000 * hb, fst=f, p_perm=perm_p(ridx[regs[i]], ridx[regs[j]], f)))
pd.DataFrame(R).round(4).to_csv("table_popgen_regions.csv", index=False); print(pd.DataFrame(R).round(4).to_string(index=False))
# ---------------------------------------------------------------- countries and isolation by distance
geo = json.load(open("ne_50m_countries.geojson"))["features"]; from shapely.geometry import shape
NAMES = {"USA": "United States of America", "Dominican Republic": "Dominican Rep.", "French Polynesia": "Fr. Polynesia", "Viet Nam": "Vietnam", "Virgin Islands": "U.S. Virgin Is."}
polys = {f["properties"]["NAME"]: shape(f["geometry"]) for f in geo}
def centroid(c):
    g = polys[NAMES.get(c, c)]; parts = list(g.geoms) if g.geom_type == "MultiPolygon" else [g]; big = max(parts, key=lambda p: p.area); return big.centroid.y, big.centroid.x
def haversine(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (a[0], a[1], b[0], b[1])); h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2; return 6371 * 2 * np.arcsin(np.sqrt(h))
cen = {c: centroid(c) for c in am_c}; C = []
for i, c1 in enumerate(am_c):
    for c2 in am_c[i + 1:]:
        f, hb = fst_hudson(cidx[c1], cidx[c2]); C.append(dict(country1=c1, country2=c2, dxy_per_kb=1000 * hb, fst=f, km=haversine(cen[c1], cen[c2])))
Cd = pd.DataFrame(C); Cd["fst_lin"] = Cd.fst / (1 - Cd.fst); Cd.round(5).to_csv("table_popgen_countries.csv", index=False)
k = len(am_c); pos = {c: i for i, c in enumerate(am_c)}; Dg = np.zeros((k, k)); Dk = np.zeros((k, k))
for r in Cd.itertuples(): Dg[pos[r.country1], pos[r.country2]] = Dg[pos[r.country2], pos[r.country1]] = r.fst_lin; Dk[pos[r.country1], pos[r.country2]] = Dk[pos[r.country2], pos[r.country1]] = r.km
iu = np.triu_indices(k, 1); r_obs = float(np.corrcoef(Dg[iu], Dk[iu])[0, 1]); cnt = 0
for _ in range(MANTEL_PERMS):
    p = rng.permutation(k); cnt += np.corrcoef(Dg[np.ix_(p, p)][iu], Dk[iu])[0, 1] >= r_obs
from scipy import stats
rho = float(stats.spearmanr(Dg[iu], Dk[iu])[0]); IBD = pd.DataFrame([dict(countries=k, pairs=len(iu[0]), pearson_r=r_obs, p_mantel_one_sided=(cnt + 1) / (MANTEL_PERMS + 1), spearman_rho=rho, mean_fst=float(Cd.fst.mean()), min_fst=float(Cd.fst.min()), max_fst=float(Cd.fst.max()))])
IBD.round(4).to_csv("table_popgen_ibd.csv", index=False); print(IBD.round(4).to_string(index=False))
# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3)); g = G.copy(); cols = ["#1b998b" if l == "region" else "#3d7ea6" for l in g.level]
ax[0].barh(range(len(g)), g.pi_per_kb, color=cols); ax[0].set_yticks(range(len(g))); ax[0].set_yticklabels([f"{a} (n={b})" for a, b in zip(g.group, g.n)], fontsize=8.5); ax[0].invert_yaxis(); ax[0].set_xlabel("Nucleotide diversity (differences per kb)"); ax[0].set_title("Diversity within groups")
M = np.full((k, k), np.nan)
for r in Cd.itertuples(): M[pos[r.country1], pos[r.country2]] = M[pos[r.country2], pos[r.country1]] = r.fst
im = ax[1].imshow(M, cmap="Blues", vmin=0); ax[1].set_xticks(range(k)); ax[1].set_yticks(range(k)); ax[1].set_xticklabels(am_c, rotation=60, ha="right", fontsize=8); ax[1].set_yticklabels(am_c, fontsize=8); fig.colorbar(im, ax=ax[1], fraction=0.046, label="Hudson's Fst"); ax[1].set_title("Differentiation among American countries")
ax[2].scatter(Cd.km, Cd.fst_lin, s=22, color="#e4572e", alpha=0.8); b1, b0 = np.polyfit(Cd.km, Cd.fst_lin, 1); xs = np.linspace(Cd.km.min(), Cd.km.max(), 50); ax[2].plot(xs, b0 + b1 * xs, color="#2e294e", lw=1.2)
ax[2].set_xlabel("Great-circle distance between centroids (km)"); ax[2].set_ylabel("Fst / (1 - Fst)"); ax[2].set_title("Isolation by distance"); ax[2].text(0.97, 0.95, f"Mantel r = {r_obs:.2f}, p = {(cnt + 1) / (MANTEL_PERMS + 1):.2f}", transform=ax[2].transAxes, ha="right", va="top", fontsize=9)
fig.tight_layout(); fig.savefig("figs8/fig8_popgen.png", dpi=200); print("figure written")
