"""
Chapter 6, Exercise 1 - Root-to-tip regression and a date-randomisation test on public SARS-CoV-2 genomes.

Data: near-complete SARS-CoV-2 genomes from NCBI GenBank, sampled at random within each GenBank release
month (Jan 2020 - Mar 2021); only genomes with a day-level collection date are kept. Root = the Wuhan-Hu-1 reference genome (NC_045512.2).
Everything runs in plain Python; genomes are aligned to the reference with minimap2 (mappy).
"""
import io, os, json, time, random, urllib.parse, urllib.request
import numpy as np, pandas as pd, mappy
from Bio import SeqIO
from scipy import stats
import matplotlib.pyplot as plt

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
REF_ACC, PER_MONTH, SEED = "NC_045512.2", 30, 11
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
C0, C1, C2 = "#2b6cb0", "#c05621", "#718096"

def get(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read().decode()
        except Exception: time.sleep(2 * (k + 1))
    raise RuntimeError(url)

# ------------------------------------------------------------------ 1. sample genomes
def sample_ids():
    """Sample genomes at random within each GenBank release month (Jan 2020 - Mar 2021)."""
    rnd = random.Random(SEED); ids = []
    for y, m in [(2020, k) for k in range(1, 13)] + [(2021, k) for k in range(1, 4)]:
        term = (f'"Severe acute respiratory syndrome coronavirus 2"[Organism] AND 29000:30000[SLEN] '
                f'AND {y}/{m:02d}/01:{y}/{m:02d}/28[PDAT]')
        r = json.loads(get(E + "esearch.fcgi?db=nuccore&retmode=json&retmax=2000&term=" + urllib.parse.quote(term)))
        pool = r["esearchresult"]["idlist"]; rnd.shuffle(pool); ids += pool[:PER_MONTH]; time.sleep(0.4)
    return ids

def download(cache="sars2.fasta", meta="sars2_meta.csv"):
    if os.path.exists(cache) and os.path.exists(meta): return
    ids, rows, fasta = sample_ids(), [], ""
    for i in range(0, len(ids), 40):
        ch = ids[i:i + 40]
        s = json.loads(get(E + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ch)))["result"]
        for u in s["uids"]:
            f = dict(zip(s[u]["subtype"].split("|"), s[u]["subname"].split("|")))
            rows.append(dict(acc=s[u]["accessionversion"], date=f.get("collection_date"), country=f.get("country")))
        fasta += get(E + "efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ch)); time.sleep(0.5)
    ref = get(E + f"efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={REF_ACC}")
    open(cache, "w").write(ref + fasta); pd.DataFrame(rows).to_csv(meta, index=False)

download()
seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse("sars2.fasta", "fasta")}
ref = seqs[REF_ACC]; L = len(ref)
meta = pd.read_csv("sars2_meta.csv")
day_level = meta.date.astype(str).str.match(r"^(\d{4}-\d{2}-\d{2}|\d{2}-[A-Za-z]{3}-\d{4})$")   # drop year-only and year-month dates
print(f"{int(day_level.sum())} of {len(meta)} sampled genomes have a day-level collection date")
meta["date"] = pd.to_datetime(meta.date.where(day_level), errors="coerce", format="mixed")
full = meta.dropna(subset=["date"]).copy()
full = full[(full.date >= "2019-12-01") & (full.date <= "2021-03-31")]
full = full[full.acc.isin(seqs) & (full.acc != REF_ACC)]
full["n_frac"] = [sum(c not in "ACGT" for c in seqs[a]) / len(seqs[a]) for a in full.acc]
full = full[full.n_frac < 0.02].drop_duplicates("acc").sort_values("date")
print(f"{len(full)} genomes with day-level dates and <2% ambiguous bases, "
      f"collected {full.date.min():%d %b} to {full.date.max():%d %b %Y}")

# ------------------------------------------------ 2. reference-anchored alignment and distances
aligner = mappy.Aligner(seq=ref, preset="asm5")
code = {"A": 0, "C": 1, "G": 2, "T": 3}
refcode = np.array([code.get(c, 4) for c in ref])
def to_row(s):
    row = np.full(L, 4, dtype=np.int8)
    for h in aligner.map(s):
        q, r = h.q_st, h.r_st
        for n, op in h.cigar:
            if op in (0, 7, 8): row[r:r + n] = [code.get(c, 4) for c in s[q:q + n]]; q += n; r += n
            elif op == 1: q += n
            elif op == 2: r += n
    return row
A = np.vstack([to_row(seqs[a]) for a in full.acc])
called = (A < 4) & (refcode < 4)
diffs = ((A != refcode) & called).sum(axis=1)
p = diffs / called.sum(axis=1)
dist = -0.75 * np.log(1 - 4.0 / 3.0 * p)                       # Jukes-Cantor distance per site
t = ((full.date - pd.Timestamp("2019-12-26")).dt.days / 365.25).to_numpy()  # years since the reference date (approx. 26 Dec 2019)

# ----------------------------------------------- 3. root-to-tip regression
fit = stats.linregress(t, dist)
print(f"root-to-tip: slope {fit.slope:.2e} substitutions/site/year (SE {fit.stderr:.1e}), "
      f"= {fit.slope * L:.1f} substitutions per genome per year, R2 = {fit.rvalue**2:.2f}, n = {len(t)}")
print(f"x-intercept (implied date of the root): {pd.Timestamp('2019-12-26') + pd.Timedelta(days=-fit.intercept / fit.slope * 365.25):%d %b %Y}")

# ----------------------------------------------- 4. date-randomisation test
rng = np.random.default_rng(1)
perm = np.array([stats.linregress(rng.permutation(t), dist).slope for _ in range(1000)])
pval = (np.sum(perm >= fit.slope) + 1) / (len(perm) + 1)
print(f"date randomisation (1000 permutations): observed slope exceeds {np.mean(perm < fit.slope):.1%} of permuted slopes "
      f"(permuted slope range {perm.min():.1e} to {perm.max():.1e}), one-sided p = {pval:.3f}")

# ------------------------------------------------------------------ figure
os.makedirs("figs6", exist_ok=True)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
ax[0].scatter(full.date, dist * L, s=16, color=C0, alpha=0.7)
xs = np.array([t.min(), t.max()])
ax[0].plot(pd.Timestamp("2019-12-26") + pd.to_timedelta(xs * 365.25, unit="D"), fit.intercept * L + fit.slope * L * xs, color=C1, lw=2)
ax[0].set(ylabel="Substitutions from the reference (per genome)", title=f"Root-to-tip: {fit.slope * L:.0f} per genome per year, R\u00b2 = {fit.rvalue**2:.2f}")
ax[0].tick_params(axis="x", rotation=30)
ax[1].hist(perm * L, bins=30, color=C2, alpha=0.8)
ax[1].axvline(fit.slope * L, color=C1, lw=2)
ax[1].set(xlabel="Slope (substitutions per genome per year)", ylabel="Permutations", title="Dates shuffled 1000 times")
ax[1].text(fit.slope * L, ax[1].get_ylim()[1] * 0.92, " observed", color=C1, ha="right")
fig.tight_layout(); fig.savefig("figs6/fig_ch6_root_to_tip.png", dpi=200); plt.close(fig)
full.assign(distance_per_site=dist, substitutions=diffs).to_csv("sars2_root_to_tip_table.csv", index=False)
print("figure written to figs6/")
