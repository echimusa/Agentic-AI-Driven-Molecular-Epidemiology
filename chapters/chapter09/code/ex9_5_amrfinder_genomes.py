"""
Chapter 9, Hands-on 5 - Call resistance determinants from genomes with AMRFinderPlus and compare with the public record.
Picks 8 Salmonella isolates (a multi-class and a determinant-free isolate from each sector, seeded), downloads their assemblies with the NCBI Datasets API,
runs AMRFinderPlus (-O Salmonella, so point mutations are called), and compares the elements with those NCBI Pathogen Detection reports for the same isolate.
Needs: the AMRFinderPlus 4.2 binaries and database in ./amrf (see README), pandas.  Results are cached per genome, so the script can be rerun to finish.
"""
import os, io, re, zipfile, subprocess, time, urllib.request
import numpy as np, pandas as pd
from amr9 import *

rng = np.random.default_rng(5); os.makedirs("genomes", exist_ok=True)
f = pd.read_csv("sal_us_frame.csv.gz"); f["genes"] = f.genes.fillna("")
f = f[f.asm_acc.notna() & f.asm_acc.astype(str).str.startswith("GCA_")]
picks = []
for sec in ["Human", "Food", "Animal", "Environment"]:
    s = f[f.sector == sec]
    for label, pool in [("multi-class", s[s.n_classes >= 4]), ("no determinant", s[s.n_classes == 0])]:
        picks.append(pool.iloc[[int(rng.integers(len(pool)))]].assign(profile=label))
P = pd.concat(picks).reset_index(drop=True); print(P[["asm_acc", "sector", "profile", "serovar_clean", "year", "n_classes"]].to_string(index=False))

def fetch(acc):
    fa = f"genomes/{acc}.fna"
    if os.path.exists(fa): return fa
    u = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
    z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(u, timeout=120).read()))
    name = [n for n in z.namelist() if n.endswith((".fna", ".fa")) and "/data/" in n][0]; open(fa, "wb").write(z.read(name)); return fa

def run(acc, fa):
    out = f"genomes/{acc}.amrfinder.tsv"
    if not os.path.exists(out):
        subprocess.run(["./amrf/amrfinder", "-n", fa, "-O", "Salmonella", "-d", "amrf/db", "--threads", "1", "-o", out], check=True, capture_output=True)
    return pd.read_csv(out, sep="\t")

def mine(df):
    d = df[(df["Type"] == "AMR") & (df["Scope"] == "core") & ~df["Method"].str.startswith("PARTIAL") & ~df["Method"].str.startswith("INTERNAL_STOP")]
    return sorted(set(d["Element symbol"]))

rows, t0 = [], time.time()
for r in P.itertuples():
    if time.time() - t0 > 170 and not os.path.exists(f"genomes/{r.asm_acc}.amrfinder.tsv"): print("time budget reached; rerun to finish"); continue
    fa = fetch(r.asm_acc); df = run(r.asm_acc, fa); a = set(mine(df)); b = set(g for g in r.genes.split(";") if g)
    rows.append(dict(asm_acc=r.asm_acc, sector=r.sector, profile=r.profile, serovar=r.serovar_clean, year=int(r.year), mb=round(sum(len(l) for l in open(fa) if not l.startswith(">")) / 1e6, 2), 
                     n_mine=len(a), n_ncbi=len(b), n_both=len(a & b), only_mine=";".join(sorted(a - b)), only_ncbi=";".join(sorted(b - a)), identical=(a == b)))
R = pd.DataFrame(rows); R.to_csv("table_amrfinder_check.csv", index=False)
print(R[["asm_acc", "sector", "profile", "serovar", "mb", "n_mine", "n_ncbi", "n_both", "only_mine", "only_ncbi"]].to_string(index=False))
print(f"\nidentical element sets in {int(R.identical.sum())} of {len(R)} genomes; elements called by both: {int(R.n_both.sum())} of {int(R.n_ncbi.sum())} reported by NCBI and {int(R.n_mine.sum())} found here")
