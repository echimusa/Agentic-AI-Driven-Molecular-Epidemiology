"""
Chapter 9, Hands-on 5 (extension) - two detection tools on the same eight genomes: AMRFinderPlus (NCBI) and RGI with the CARD database.
AMRFinderPlus: the output of ex9_5 (core AMR elements that are not partial hits or internal stops, as in ex9_5).
RGI 6.0.8, CARD 4.0.2:  rgi main -i genome.fna -o out -t contig -a DIAMOND --local --clean   (default: Perfect and Strict hits; Loose hits are not reported).
Step 1 counts what each tool reports.  Step 2 sets aside as 'background' every RGI hit seen in at least one of the four genomes in which AMRFinderPlus finds nothing.
Step 3 compares the remaining hits with the AMRFinderPlus elements after the two naming schemes are made comparable.
Step 4 looks the AMRFinderPlus-only elements up in a second RGI run that also reports Loose hits (--include_loose), for the four genomes that have elements.
Input: amrfinder_outputs/*.amrfinder.tsv, rgi_outputs/*.rgi.txt, rgi_outputs/*.loose.rgi.txt.   Output: table_two_tools.csv, table_two_tools_background.csv, table_two_tools_summary.csv, table_two_tools_loose.csv.
"""
import glob, os, re
import numpy as np, pandas as pd
accs = sorted(os.path.basename(p)[:15] for p in glob.glob("amrfinder_outputs/*.amrfinder.tsv")); assert len(accs) == 8
def norm(s):
    s = s.strip().lower(); s = re.sub(r"^bla", "", s); return s
def amr_elements(acc):
    d = pd.read_csv(f"amrfinder_outputs/{acc}.amrfinder.tsv", sep="\t")
    d = d[(d["Type"] == "AMR") & (d["Scope"] == "core") & ~d["Method"].str.startswith("PARTIAL") & ~d["Method"].str.startswith("INTERNAL_STOP")]
    return sorted(set(d["Element symbol"]))
def rgi_hits(acc):
    d = pd.read_csv(f"rgi_outputs/{acc}.rgi.txt", sep="\t"); out = []
    for r in d.itertuples():
        name, snp = r.Best_Hit_ARO, r.SNPs_in_Best_Hit_ARO if isinstance(r.SNPs_in_Best_Hit_ARO, str) else ""
        if "variant" in r.Model_type and snp:
            m = re.search(r"\b([a-z]{3}[A-Z]|[A-Z][a-z]{2,3}[A-Z]?|EF-Tu|PBP3|UhpT)\b(?= (?:conferring|mutants|with))", name); gene = m.group(1) if m else name.split()[-1]
            out.append(f"{gene}_{snp.split(',')[0].strip()}" if "," not in snp else f"{gene}_{snp.replace(', ', '+')}")
        else: out.append(name)
    return out
rows, rgi_all, amr_all = [], {}, {}
for acc in accs:
    a = amr_elements(acc); r = rgi_hits(acc); amr_all[acc] = a; rgi_all[acc] = sorted(set(r)); rows.append(dict(assembly=acc, amrfinder_elements=len(a), rgi_hits=len(r), rgi_distinct=len(set(r))))
S1 = pd.DataFrame(rows); print(S1.to_string(index=False))
empty = [acc for acc in accs if not amr_all[acc]]; background = set.union(*[set(rgi_all[acc]) for acc in empty]); print(f"\n{len(empty)} genomes with no AMRFinderPlus element; distinct RGI hits seen in at least one of them: {len(background)}")
every = set.intersection(*[set(v) for v in rgi_all.values()]); print(f"RGI hits present in all eight genomes: {len(every)}")
bgrows = [dict(background_hit=h, in_genomes_without_amrfinder_element=sum(h in rgi_all[a] for a in empty), in_all_eight=int(h in every)) for h in sorted(background)]; pd.DataFrame(bgrows).to_csv("table_two_tools_background.csv", index=False)
cmp = []
for acc in [a for a in accs if amr_all[a]]:
    A = {norm(x): x for x in amr_all[acc]}; R = {norm(x): x for x in rgi_all[acc] if x not in background}
    fam = lambda k: re.match(r"([a-z]+)", k).group(1)
    for k in sorted(set(A) | set(R)):
        if k in A and k in R: cmp.append(dict(assembly=acc, amrfinder=A[k], rgi=R[k], result="same name"))
        elif k in A:
            alt = [rk for rk in R if rk not in A and fam(rk) == fam(k) and k.startswith(fam(k)) and fam(k) in ("tem", "ctx", "cmy", "oxa", "shv")]
            cmp.append(dict(assembly=acc, amrfinder=A[k], rgi=R[alt[0]] if alt else "", result="same family, different allele" if alt else "AMRFinderPlus only"))
        elif not any(fam(k) == fam(ak) and fam(k) in ("tem", "ctx", "cmy", "oxa", "shv") and ak not in R for ak in A):
            cmp.append(dict(assembly=acc, amrfinder="", rgi=R[k], result="RGI only"))
C = pd.DataFrame(cmp); C.to_csv("table_two_tools.csv", index=False); print("\n" + C.result.value_counts().to_string()); print(C[C.result != "same name"].to_string(index=False))
tot = len(C); same = int((C.result == "same name").sum()); fam_ = int((C.result == "same family, different allele").sum()); a_only = int((C.result == "AMRFinderPlus only").sum()); r_only = int((C.result == "RGI only").sum())
n_amr = sum(len(v) for v in amr_all.values()); pd.DataFrame([dict(genomes=8, amrfinder_elements=n_amr, rgi_hits=int(S1.rgi_hits.sum()), rgi_background=len(background), compared=tot, same_name=same, same_family_other_allele=fam_, amrfinder_only=a_only, rgi_only=r_only)]).to_csv("table_two_tools_summary.csv", index=False)
print(f"\nAMRFinderPlus {n_amr} elements; RGI {int(S1.rgi_hits.sum())} hits ({len(background)} background hits); after removing the background: {same} same name, {fam_} same family with a different allele, {a_only} AMRFinderPlus only, {r_only} RGI only")

# ---------------------------------------------------------------- step 4: are the AMRFinderPlus-only elements Loose RGI hits?
lo = []
for r in C[C.result == "AMRFinderPlus only"].itertuples():
    d = pd.read_csv(f"rgi_outputs/{r.assembly}.loose.rgi.txt", sep="\t"); h = d[d.Best_Hit_ARO.map(norm) == norm(r.amrfinder)].sort_values("Best_Identities", ascending=False)
    lo.append(dict(assembly=r.assembly, element=r.amrfinder, loose_hit=h.Best_Hit_ARO.iloc[0] if len(h) else "", cut_off=h.Cut_Off.iloc[0] if len(h) else "", identity_pct=round(float(h.Best_Identities.iloc[0]), 1) if len(h) else np.nan,
                   length_pct=round(float(h["Percentage Length of Reference Sequence"].iloc[0]), 1) if len(h) else np.nan, loose_hits_in_genome=len(d), strict_or_perfect_in_genome=int((d.Cut_Off != "Loose").sum())))
L = pd.DataFrame(lo); L.to_csv("table_two_tools_loose.csv", index=False); print("\n" + L.to_string(index=False))
