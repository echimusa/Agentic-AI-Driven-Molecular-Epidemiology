"""
Chapter 12, Hands-on 2 - Representativeness and equity: genomes per million people, a Lorenz curve and a Gini coefficient.
Isolates: NCBI Pathogen Detection Salmonella (sal_surv.pkl). Population and income group: World Bank open data (API, cached to JSON).
Output: table_country_equity.csv, table_income_groups.csv, Figure 2.
"""
import json, os, urllib.request, difflib, re
import numpy as np, pandas as pd
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
def wb(url, cache):
    if not os.path.exists(cache): open(cache, "w").write(urllib.request.urlopen(url, timeout=60).read().decode())
    return json.load(open(cache))[1]
C = pd.DataFrame([{"iso3": c["id"], "wb_name": c["name"], "region": c["region"]["value"], "income": c["incomeLevel"]["value"]} for c in wb("https://api.worldbank.org/v2/country?format=json&per_page=400", "wb_countries.json")])
C = C[C.region != "Aggregates"]
P = pd.DataFrame([{"iso3": r["countryiso3code"], "pop": r["value"]} for r in wb("https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL?format=json&per_page=400&mrv=1", "wb_pop.json") if r["value"]])
C = C.merge(P, on="iso3", how="left").dropna(subset=["pop"]); print(f"{len(C)} World Bank economies with a population; world population in the table: {C['pop'].sum() / 1e9:.2f} billion")

# ---------------------------------------------------------------- match Pathogen Detection country names to World Bank names
m = load_surv(); m = m[m.country.notna() & m.year.between(2015, 2024)]
ALIAS = {"USA": "United States", "Viet Nam": "Viet Nam", "Vietnam": "Viet Nam", "South Korea": "Korea, Rep.", "North Korea": "Korea, Dem. People's Rep.", "Russia": "Russian Federation", "Iran": "Iran, Islamic Rep.", "Egypt": "Egypt, Arab Rep.",
         "Turkey": "Turkiye", "Czech Republic": "Czechia", "Slovakia": "Slovak Republic", "Syria": "Syrian Arab Republic", "Venezuela": "Venezuela, RB", "Yemen": "Yemen, Rep.", "Gambia": "Gambia, The", "Bahamas": "Bahamas, The", "Laos": "Lao PDR",
         "Kyrgyzstan": "Kyrgyz Republic", "Hong Kong": "Hong Kong SAR, China", "Macau": "Macao SAR, China", "Democratic Republic of the Congo": "Congo, Dem. Rep.", "Republic of the Congo": "Congo, Rep.", "Ivory Coast": "Cote d'Ivoire", "Cote d'Ivoire": "Cote d'Ivoire",
         "Micronesia": "Micronesia, Fed. Sts.", "Saint Lucia": "St. Lucia", "Saint Kitts and Nevis": "St. Kitts and Nevis", "Saint Vincent and the Grenadines": "St. Vincent and the Grenadines", "Brunei": "Brunei Darussalam", "Cape Verde": "Cabo Verde", "Swaziland": "Eswatini", "Burma": "Myanmar", "Palestine": "West Bank and Gaza", "Curacao": "Curacao", "Macedonia": "North Macedonia"}
names = C.wb_name.tolist(); lower = {n.lower(): n for n in names}
def match(n):
    n = ALIAS.get(n, n)
    if n in names: return n
    if n.lower() in lower: return lower[n.lower()]
    g = difflib.get_close_matches(n, names, n=1, cutoff=0.88); return g[0] if g else None
cnt = m.groupby("country").size(); mp = {c: match(c) for c in cnt.index}; unmatched = cnt[[c for c in cnt.index if mp[c] is None]].sort_values(ascending=False)
print(f"isolates collected 2015-2024 with a country: {cnt.sum():,}; matched to a World Bank economy: {cnt[[c for c in cnt.index if mp[c]]].sum():,} ({cnt[[c for c in cnt.index if mp[c]]].sum() / cnt.sum():.1%}); unmatched labels: {len(unmatched)} (largest: {unmatched.head(5).to_dict()})")
cnt = cnt.rename(index=mp); cnt = cnt[cnt.index.notna()].groupby(level=0).sum(); C["isolates"] = C.wb_name.map(cnt).fillna(0).astype(int); C["per_million"] = C.isolates / (C["pop"] / 1e6)
C.sort_values("per_million", ascending=False).to_csv("table_country_equity.csv", index=False)

# ---------------------------------------------------------------- concentration and inequality
tot = C.isolates.sum(); top = C.sort_values("isolates", ascending=False); print(f"economies with at least one isolate: {(C.isolates > 0).sum()} of {len(C)}; share of the world population living in economies with none: {C[C.isolates == 0]['pop'].sum() / C['pop'].sum():.1%}")
print("top 3 economies hold", f"{top.isolates.head(3).sum() / tot:.1%}", "of isolates:", top.head(3)[["wb_name", "isolates"]].values.tolist()); print("median isolates per million:", round(C.per_million.median(), 2), "| mean:", round(C.isolates.sum() / (C['pop'].sum() / 1e6), 2))
def lorenz(counts, pops):
    o = np.argsort(counts / pops); c, p = counts[o], pops[o]; L = np.r_[0, np.cumsum(c) / c.sum()]; Pcum = np.r_[0, np.cumsum(p) / p.sum()]; gini = 1 - np.sum((Pcum[1:] - Pcum[:-1]) * (L[1:] + L[:-1]))
    return Pcum, L, gini
Pc, L, G = lorenz(C.isolates.values.astype(float), C["pop"].values.astype(float)); print(f"Gini coefficient of isolates over population: {G:.3f}")
rng = np.random.default_rng(12); boots = [lorenz(*(lambda i: (C.isolates.values[i].astype(float), C["pop"].values[i].astype(float)))(rng.integers(0, len(C), len(C))))[2] for _ in range(500)]; print(f"bootstrap 95% interval over economies: {np.percentile(boots, 2.5):.3f} to {np.percentile(boots, 97.5):.3f}")
I = C.groupby("income").agg(economies=("iso3", "size"), population=("pop", "sum"), isolates=("isolates", "sum")); I["pop_share"] = 100 * I.population / I.population.sum(); I["isolate_share"] = 100 * I.isolates / I.isolates.sum(); I["per_million"] = I.isolates / (I.population / 1e6)
I = I.reindex(["High income", "Upper middle income", "Lower middle income", "Low income"]); I.round(2).to_csv("table_income_groups.csv"); print(I.round(2).to_string())

# ---------------------------------------------------------------- Figure 2
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.7), gridspec_kw={"width_ratios": [1, 1.1, 1.15]})
ax[0].plot(Pc, L, color="#2b6cb0", lw=2.2); ax[0].plot([0, 1], [0, 1], color="black", lw=0.8, ls="--"); ax[0].fill_between(Pc, L, Pc, color="#2b6cb0", alpha=0.12); ax[0].text(0.05, 0.82, f"Gini = {G:.2f}", transform=ax[0].transAxes, fontsize=10.5)
ax[0].set(xlabel="Cumulative share of world population", ylabel="Cumulative share of isolates", title="Lorenz curve", xlim=(0, 1), ylim=(0, 1)); ax[0].text(0.05, 0.74, "economies ordered by isolates per person", transform=ax[0].transAxes, fontsize=8.5, color="#4a5568")
x = np.arange(len(I)); ax[1].bar(x - 0.2, I.pop_share, 0.4, color="#a0aec0", label="Share of population"); ax[1].bar(x + 0.2, I.isolate_share, 0.4, color="#2b6cb0", label="Share of isolates"); ax[1].set_xticks(x); ax[1].set_xticklabels(["High", "Upper\nmiddle", "Lower\nmiddle", "Low"]); ax[1].set(ylabel="%", title="By World Bank income group"); ax[1].legend(frameon=False, fontsize=9)
for i, r in enumerate(I.itertuples()): ax[1].text(i + 0.2, r.isolate_share + 1, f"{r.isolate_share:.1f}", ha="center", fontsize=8.5)
tp = C[C["pop"] > 5e6].sort_values("per_million", ascending=False); show = pd.concat([tp.head(8), tp.tail(5)]); cols = ["#2b6cb0"] * 8 + ["#dd6b20"] * 5; ax[2].barh(range(len(show)), np.where(show.per_million > 0, show.per_million, np.nan), color=cols); ax[2].set_xscale("log"); ax[2].set_xlim(0.01, 3000)
for i_, v_ in enumerate(show.per_million.values):
    if v_ == 0: ax[2].text(0.012, i_, "no isolates", va="center", fontsize=8.5, color="#c05621"); ax[2].set_yticks(range(len(show))); ax[2].set_yticklabels(show.wb_name, fontsize=8.5); ax[2].invert_yaxis()
ax[2].set(xlabel="Isolates per million people, 2015-2024 (log scale)", title="Highest and lowest (population above 5 million)"); fig.tight_layout(); fig.savefig("figs12/fig2_equity.png", dpi=200); plt.close(fig); print("figure written")
