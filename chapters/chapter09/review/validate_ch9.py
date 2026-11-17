"""
Chapter 9 validation harness. Recomputes the numbers quoted in the chapter with code that does not import the chapter's own scripts (amr9.py and ex9_*.py).
A  the original chapter: the frame, genotypic MDR, crude odds ratios, the mixed-cluster counts and a short permutation null, genotype--phenotype agreement, the eight AMRFinderPlus outputs
B  the two tools: RGI and AMRFinderPlus counts, the background hits, the comparison, the Loose hits
C  why genotype and phenotype disagree: the false-susceptible calls, the isolates with no determinant, azithromycin
D  the SNP thresholds: every pair distance recomputed with Bio.Phylo, the threshold table recounted
E  the metadata completeness table
F  the study audit
G  the chapter text (tables, figures, objectives, sentences, references)
USAGE   python validate_ch9.py [--chapter Chapter_9_....docx] [--doi] [--full]
        --chapter  also checks the chapter's tables, figures, objectives, sentences and references against the results (needs pandoc and python-docx)
        --doi      also resolves every DOI of the chapter against Crossref (needs the network)
        --full     also reruns every script that can run from the shipped data in a scratch folder and checks that its result tables are byte-identical to the shipped ones
                   (needs pandas, numpy, scipy, statsmodels, networkx, matplotlib, Biopython and R with ggplot2; about 2 minutes), and refits the R-companion model with statsmodels
Prints PASS/FAIL/SKIP per check, writes validation_results_ch9.csv and exits 1 if any check fails.
"""
import argparse, csv, glob, hashlib, math, os, re, shutil, subprocess, sys, tempfile, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--doi", action="store_true"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
def P(name):
    for d in (".", "../data", "data", os.path.join(HERE, "../data"), os.path.join(HERE, "data")):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []; EXPECTED_DOIS = 32
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if bool(ok) else "FAIL", detail=str(detail)[:170])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(round(float(v), 10))).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
SECT = ["Human", "Food", "Animal", "Environment"]
# ================================================================== I. full reruns (run first, while this process is small: ex9_2 peaks at about 2.5 GB)
tmp = None
if A.full:
    tmp = tempfile.mkdtemp(prefix="ch9_"); codes = os.path.dirname(os.path.abspath(P("sal_us_frame.csv.gz"))); src = os.path.join(HERE, "../code") if os.path.isdir(os.path.join(HERE, "../code")) else HERE
    for f_ in glob.glob(os.path.join(src, "*.py")) + glob.glob(os.path.join(src, "*.R")): shutil.copy(f_, tmp)
    for item in os.listdir(codes):
        s_, d_ = os.path.join(codes, item), os.path.join(tmp, item)
        if item.startswith("table_"): continue
        shutil.copytree(s_, d_) if os.path.isdir(s_) else shutil.copy(s_, d_)
    os.makedirs(os.path.join(tmp, "figs9"), exist_ok=True); steps = [("ex9_2", ["python3", "ex9_2_resistance_by_sector.py"]), ("ex9_2b", ["Rscript", "ex9_2b_glm.R"]), ("ex9_3", ["python3", "ex9_3_genotype_phenotype.py"]), ("ex9_4", ["python3", "ex9_4_clusters_cooccurrence.py"]), ("ex9_3b", ["python3", "ex9_3b_discordance.py"]),
        ("ex9_5b", ["python3", "ex9_5b_two_tools.py"]), ("ex9_5c", ["Rscript", "ex9_5c_summarise.R"]), ("ex9_6", ["python3", "ex9_6_snp_thresholds.py"]), ("ex9_7", ["python3", "ex9_7_metadata_completeness.py"]), ("ex9_8", ["python3", "ex9_8_study_audit.py"]), ("ex9_9", ["python3", "ex9_9_new_figures.py"])]
    for name, cmd in steps:
        r_ = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True); check(f"I.{name}", f"{name} reruns from the shipped data without error", r_.returncode == 0, (r_.stderr or "")[-120:])
    tabs_ = sorted(os.path.basename(f_) for f_ in glob.glob(os.path.join(codes, "table_*.csv"))); diffs = []
    for tb in tabs_:
        a_, b_ = os.path.join(tmp, tb), os.path.join(codes, tb)
        if tb == "table_sampling_frame.csv" or tb == "table_serovar_mix.csv" or tb == "table_amrfinder_check.csv": continue
        if not os.path.exists(a_) or hashlib.md5(open(a_, "rb").read()).hexdigest() != hashlib.md5(open(b_, "rb").read()).hexdigest(): diffs.append(tb)
    check("I.tables", f"the {len(tabs_) - 3} result tables of the rerun are byte-identical to the shipped ones", not diffs, diffs)

# ================================================================== A. the original chapter
fr = pd.read_csv(P("sal_us_frame.csv.gz"), dtype={"state": str}); cnt = fr.sector.value_counts()
check("A.1", "frame of 387,635 US isolates 2015-2024: 293,419 human, 60,909 food, 23,384 animal, 9,923 environmental; susceptibility results 21, 4,217, 568 and 0",
      len(fr) == 387635 and [int(cnt[s]) for s in SECT] == [293419, 60909, 23384, 9923] and [int(fr[fr.sector == s].has_ast.sum()) for s in SECT] == [21, 4217, 568, 0], (len(fr), cnt.to_dict()))
ncl = fr.classes.fillna("").map(lambda s: len([c for c in s.split(";") if c])); fr["mdr"] = (ncl >= 3).astype(int); crude = fr.groupby("sector").mdr.mean()
check("A.2", "crude genotypic MDR is 8.7% in human, 25.4% in food, 21.5% in animal and 7.2% in environmental isolates", [rnd(100 * crude[s], 1) for s in SECT] == ["8.7", "25.4", "21.5", "7.2"], crude.round(4).to_dict())
def odds(a, b): return (a.mdr.sum() / (len(a) - a.mdr.sum())) / (b.mdr.sum() / (len(b) - b.mdr.sum()))
hum = fr[fr.sector == "Human"]; ors = [rnd(odds(fr[fr.sector == s], hum), 2) for s in ("Food", "Animal", "Environment")]
check("A.3", "crude odds ratios of genotypic MDR against human isolates: food 3.59, animal 2.88, environment 0.81", ors == ["3.59", "2.88", "0.81"], ors)
cl = pd.read_csv(P("salmonella_cluster_list.tsv"), sep="\t", usecols=["PDS_acc", "target_acc"]); j = fr.merge(cl, on="target_acc", how="inner"); j = j[j.PDS_acc.map(j.PDS_acc.value_counts()) >= 2]
g = j.groupby("PDS_acc").sector.agg(lambda s: set(s)); anyn = sum(1 for s in g if "Human" in s and len(s - {"Human"}) > 0); hf = sum(1 for s in g if {"Human", "Food"} <= s); ha = sum(1 for s in g if {"Human", "Animal"} <= s)
check("A.4", "351,462 isolates sit in 15,166 clusters with two or more isolates in the frame; clusters with a human and any non-human isolate 2,707, human with food 1,806, human with animal 1,531",
      len(j) == 351462 and len(g) == 15166 and (anyn, hf, ha) == (2707, 1806, 1531), (len(j), len(g), anyn, hf, ha))
rng = np.random.default_rng(2026); j = j.reset_index(drop=True); cid = pd.factorize(j.PDS_acc)[0]; ncid = cid.max() + 1; strata = pd.factorize(j.serovar_clean.astype(str) + "|" + j.year.astype(int).astype(str))[0]
isH = (j.sector == "Human").values; groups = [np.where(strata == s)[0] for s in range(strata.max() + 1)]; vals = []
for _ in range(30):
    h = isH.copy()
    for gg in groups: h[gg] = isH[rng.permutation(gg)]
    nh = np.bincount(cid, weights=h.astype(float), minlength=ncid); n_ = np.bincount(cid, minlength=ncid); vals.append(int(((nh > 0) & (nh < n_)).sum()))
check("A.5", "an independent permutation null (30 shuffles of the human label within serovar and year) gives 6,058 to 6,213 clusters with a human and a non-human isolate, far above the 2,707 observed", min(vals) > 5950 and max(vals) < 6300 and 2707 < min(vals), (min(vals), float(np.mean(vals)), max(vals)))
# genotype-phenotype: own parser and own rules from the catalogue
cat = {}
with open(P("amrf_db_catalog/ReferenceGeneCatalog.txt"), encoding="utf8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["type"] == "AMR":
            for k in (r["gene_family"], r["allele"]):
                if k: cat[k] = (r["class"], r["subclass"])
with open(P("amrf_db_catalog/AMRProt-mutation.tsv"), encoding="utf8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["#taxgroup"] == "Salmonella": cat[r["standard_mutation_symbol"]] = (r["class"], r["subclass"])
def toks(s):
    if not isinstance(s, str): return []
    out = []
    for t_ in s.replace('"', "").split(","):
        n, _, suf = t_.partition("=")
        if n and suf in ("COMPLETE", "POINT", "ALLELE", "HMM", "EXACT", "BLAST", ""): out.append(n)
    return out
RULE = {"ampicillin": ("BETA-LACTAM", None, 0), "amoxicillin-clavulanic acid": ("BETA-LACTAM", "CEPHALOSPORIN", 0), "cefoxitin": ("BETA-LACTAM", "CEPHALOSPORIN", 0), "ceftiofur": ("BETA-LACTAM", "CEPHALOSPORIN", 0), "ceftriaxone": ("BETA-LACTAM", "CEPHALOSPORIN", 0),
        "streptomycin": ("AMINOGLYCOSIDE", "STREPTOMYCIN", 0), "gentamicin": ("AMINOGLYCOSIDE", "GENTAMICIN", 0), "kanamycin": ("AMINOGLYCOSIDE", "KANAMYCIN", 0), "amikacin": ("AMINOGLYCOSIDE", "AMIKACIN", 0), "tetracycline": ("TETRACYCLINE", None, 0),
        "chloramphenicol": ("PHENICOL", None, 0), "sulfamethoxazole": ("SULFONAMIDE", None, 0), "trimethoprim-sulfamethoxazole": ("TRIMETHOPRIM", None, 0), "nalidixic acid": ("QUINOLONE", None, 1), "ciprofloxacin": ("QUINOLONE", None, 0),
        "azithromycin": ("MACROLIDE", None, 0), "sulfisoxazole": ("SULFONAMIDE", None, 0), "meropenem": ("BETA-LACTAM", "CARBAPENEM", 0), "colistin": ("COLISTIN", None, 0)}
def hit(genes, drug):
    cls, sub, pt = RULE[drug]
    return any(g in cat and cat[g][0] == cls and (sub is None or sub in cat[g][1]) and (not pt or "_" in g) for g in genes)
ast = pd.read_csv(P("sal_ast_subset.csv.gz"), low_memory=False); ast["g"] = ast.AMR_genotypes_core.map(toks)
def ast_of(s):
    d = {}
    if isinstance(s, str):
        for t_ in s.replace('"', "").split(","):
            if "=" in t_: k, v = t_.rsplit("=", 1); d[k.strip()] = v.strip()
    return d
ast["ast"] = ast.AST_phenotypes.map(ast_of); agree = tot = 0; per = {}
for drug in RULE:
    o = ast.ast.map(lambda d: d.get(drug)); ok = o.isin(["S", "R"])
    if ok.sum() < 30: continue
    pr = ast.g[ok].map(lambda g_: hit(g_, drug)); y = o[ok] == "R"; tp, fn, fp, tn = int((pr & y).sum()), int((~pr & y).sum()), int((pr & ~y).sum()), int((~pr & ~y).sum()); per[drug] = (tp, fn, fp, tn); agree += tp + tn; tot += tp + fn + fp + tn
check("A.6", "genotype predicts the broth-microdilution result in 135,892 of 138,111 drug-isolate tests (98.4%) over 19 drugs and 9,977 isolates", len(ast) == 9977 and len(per) == 19 and (agree, tot) == (135892, 138111) and rnd(100 * agree / tot, 1) == "98.4", (len(per), agree, tot))
tp, fn, fp, tn = per["tetracycline"]; az = per["azithromycin"]; st = per["streptomycin"]
check("A.7", "tetracycline sensitivity 98.4% and specificity 99.3% over 9,816 tests; azithromycin 14 of 42 resistant isolates with a determinant (33%); streptomycin specificity 93.5% with 352 false-resistant calls",
      rnd(100 * tp / (tp + fn), 1) == "98.4" and rnd(100 * tn / (tn + fp), 1) == "99.3" and tp + fn + fp + tn == 9816 and (az[0], az[0] + az[1]) == (14, 42) and rnd(100 * az[0] / 42, 0) == "33" and rnd(100 * st[3] / (st[3] + st[2]), 1) == "93.5" and st[2] == 352, (per["tetracycline"], az, st))
shp = sorted(glob.glob(P("amrfinder_outputs") + "/*.amrfinder.tsv")); sets = {}
for f_ in shp:
    d = pd.read_csv(f_, sep="\t"); d = d[(d["Type"] == "AMR") & (d["Scope"] == "core") & ~d["Method"].str.startswith("PARTIAL") & ~d["Method"].str.startswith("INTERNAL_STOP")]; sets[os.path.basename(f_)[:15]] = sorted(set(d["Element symbol"]))
same = 0
for acc, s in sets.items():
    row = fr[fr.asm_acc == acc]
    if len(row) and set(s) == set(str(row.genes.iloc[0]).split(";") if isinstance(row.genes.iloc[0], str) else []): same += 1
check("A.8", "AMRFinderPlus on eight assemblies: 23 elements, and for every one of the eight the set equals the public record in the frame", len(sets) == 8 and sum(len(v) for v in sets.values()) == 23 and same == 8, (len(sets), sum(len(v) for v in sets.values()), same))
# ================================================================== B. the two tools
rg = {os.path.basename(f_)[:15]: list(csv.DictReader(open(f_, encoding="utf8"), delimiter="\t")) for f_ in sorted(glob.glob(P("rgi_outputs") + "/GCA_*.rgi.txt")) if "loose" not in f_}
nh = [len(v) for v in rg.values()]; cut = pd.Series([r["Cut_Off"] for v in rg.values() for r in v]).value_counts()
check("B.1", "RGI 6.0.8 with CARD 4.0.2 reports 261 hits on the eight genomes (37 Perfect and 224 Strict), 27 to 40 per genome; AMRFinderPlus reports 0, 5, 7, 0, 6, 5, 0 and 0 elements", sum(nh) == 261 and int(cut["Perfect"]) == 37 and int(cut["Strict"]) == 224 and min(nh) == 27 and max(nh) == 40 and [len(sets[a]) for a in sorted(sets)] == [0, 5, 7, 0, 6, 5, 0, 0], (nh, cut.to_dict()))
rows_n = 0; dup = 0
for f_ in shp:
    d = pd.read_csv(f_, sep="\t"); d = d[(d["Type"] == "AMR") & (d["Scope"] == "core") & ~d["Method"].str.startswith("PARTIAL") & ~d["Method"].str.startswith("INTERNAL_STOP")]; rows_n += len(d); dup += int((d.groupby("Element symbol").size() == 2).sum())
check("B.2", "the AMRFinderPlus output has 44 core rows for 23 distinct elements, and 21 elements are listed twice", rows_n == 44 and dup == 21, (rows_n, dup))
key = lambda r: (r["Best_Hit_ARO"], r["SNPs_in_Best_Hit_ARO"] if "variant" in r["Model_type"] else "")
empty = [a for a in sets if not sets[a]]; bg = set().union(*[{key(r) for r in rg[a]} for a in empty]); every = set.intersection(*[{key(r) for r in v} for v in rg.values()])
check("B.3", "30 distinct RGI hits are seen in the four genomes where AMRFinderPlus finds nothing, and 23 of them are in all eight genomes", len(empty) == 4 and len(bg) == 30 and len(bg & every) == 23 and len(every) == 23, (len(empty), len(bg), len(every)))
tt = pd.read_csv(P("table_two_tools.csv")); vc = tt.result.value_counts()
check("B.4", "in the four genomes with calls, 17 elements have the same name in both tools, 1 the same family with another allele, 5 are AMRFinderPlus only, 3 RGI only", [int(vc.get(k, 0)) for k in ("same name", "same family, different allele", "AMRFinderPlus only", "RGI only")] == [17, 1, 5, 3], vc.to_dict())
lo = {}
for acc in ("GCA_010670385.1", "GCA_011246915.1", "GCA_018413125.1", "GCA_019546675.1"):
    lo[acc] = list(csv.DictReader(open(P(f"rgi_outputs/{acc}.loose.rgi.txt"), encoding="utf8"), delimiter="\t"))
def best(acc, name): h = [r for r in lo[acc] if r["Best_Hit_ARO"].lower() == name.lower()]; r = max(h, key=lambda r: float(r["Best_Identities"])); return r["Cut_Off"], round(float(r["Best_Identities"]), 1), round(float(r["Percentage Length of Reference Sequence"]), 1)
got = [best("GCA_010670385.1", "floR"), best("GCA_010670385.1", "tet(G)"), best("GCA_011246915.1", "floR"), best("GCA_018413125.1", "floR"), best("GCA_019546675.1", "aph(3'')-Ib")]
check("B.5", "all five AMRFinderPlus-only elements are Loose RGI hits: floR 98.5, 99.3 and 99.3% identity over the full length, tet(G) 94.6%, aph(3'')-Ib 100% identity over 52.1% of the reference; 475 to 497 Loose-run hits per genome",
      got == [("Loose", 98.5, 100.0), ("Loose", 94.6, 100.0), ("Loose", 99.3, 100.0), ("Loose", 99.3, 100.0), ("Loose", 100.0, 52.1)] and sorted(len(v) for v in lo.values()) == [475, 482, 484, 497], (got, [len(v) for v in lo.values()]))
cc = pd.read_csv(P("table_class_counts_r.csv")); check("B.6", "the R summary counts 23 elements by class: aminoglycoside 8, beta-lactam 4, phenicol 3, sulfonamide 3, tetracycline 3, quinolone 1, trimethoprim 1", dict(zip(cc.Class, cc.elements)) == {"AMINOGLYCOSIDE": 8, "BETA-LACTAM": 4, "PHENICOL": 3, "SULFONAMIDE": 3, "TETRACYCLINE": 3, "QUINOLONE": 1, "TRIMETHOPRIM": 1}, dict(zip(cc.Class, cc.elements)))
# ================================================================== C. why genotype and phenotype disagree
fn_tot = sum(v[1] for v in per.values()); none = 0; per_none = {}
for drug in RULE:
    o = ast.ast.map(lambda d: d.get(drug)); ok = o.isin(["S", "R"])
    if ok.sum() < 30: continue
    n0 = int(((o[ok] == "R") & (ast.g[ok].map(len) == 0)).sum()); none += n0; per_none[drug] = n0
check("C.1", "823 resistant isolates are called susceptible by the rules; 208 (25.3%) of them carry no core determinant; tetracycline 37, streptomycin 59, ciprofloxacin 22, ampicillin 16, nalidixic acid 9, azithromycin 0",
      fn_tot == 823 and none == 208 and rnd(100 * none / fn_tot, 1) == "25.3" and [per_none[d] for d in ("tetracycline", "streptomycin", "ciprofloxacin", "ampicillin", "nalidixic acid", "azithromycin")] == [37, 59, 22, 16, 9, 0], (fn_tot, none))
o = ast.ast.map(lambda d: d.get("azithromycin")); azr = ast[o == "R"].copy(); azr["pred"] = azr.g.map(lambda g_: hit(g_, "azithromycin")); miss = azr[~azr.pred]; miss["year"] = pd.to_datetime(miss.collection_date, errors="coerce", format="mixed").dt.year
check("C.2", "28 azithromycin-resistant isolates carry no macrolide determinant, all carry other determinants (median 9), 18 are Infantis, 2013 to 2020, none carries acrB R717L or R717Q, and 23 were called with database 2025-12-03.1 and 5 with 2026-05-15.1",
      len(miss) == 28 and (miss.g.map(len) > 0).all() and int(miss.g.map(len).median()) == 9 and int((miss.serovar == "Infantis").sum()) == 18 and (int(miss.year.min()), int(miss.year.max())) == (2013, 2020) and not any(x.startswith("acrB_R717") for gl in miss.g for x in gl)
      and miss.refgene_db_version.value_counts().to_dict() == {"2025-12-03.1": 23, "2026-05-15.1": 5}, (len(miss), miss.refgene_db_version.value_counts().to_dict()))
azr["year"] = pd.to_datetime(azr.collection_date, errors="coerce", format="mixed").dt.year; per_ = [(int((azr[(azr.year >= a) & (azr.year <= b)].pred).sum()), len(azr[(azr.year >= a) & (azr.year <= b)])) for a, b in ((1900, 2015), (2016, 2019), (2020, 2100))]
check("C.3", "azithromycin sensitivity by period: 5 of 8 resistant isolates in 2000-2015, 8 of 31 in 2016-2019, 1 of 3 in 2020-2024; 18 susceptible isolates carry a macrolide determinant", per_ == [(5, 8), (8, 31), (1, 3)] and az[2] == 18, (per_, az[2]))
mut = [r for r in csv.DictReader(open(P("amrf_db_catalog/AMRProt-mutation.tsv"), encoding="utf8"), delimiter="\t") if r["#taxgroup"] == "Salmonella" and r["standard_mutation_symbol"].startswith("acrB_R717")]
check("C.4", "the catalogue lists two Salmonella acrB R717 variants (R717L and R717Q), both with class macrolide", sorted(r["standard_mutation_symbol"] for r in mut) == ["acrB_R717L", "acrB_R717Q"] and all(r["class"] == "MACROLIDE" for r in mut), [r["standard_mutation_symbol"] for r in mut])
# ================================================================== D. SNP thresholds
from Bio import Phylo
pr = pd.read_csv(P("table_snp_pairs.csv")); thr = pd.read_csv(P("table_snp_thresholds.csv")); ok_thr = True; det = []
for r in thr.itertuples():
    q = pr[pr.snp_distance <= r.threshold_snps]; got_ = (len(q), int((q.sector == "Food").sum()), int((q.sector == "Animal").sum()), int((q.sector == "Environment").sum()), q.cluster.nunique())
    if got_ != (r.pairs, r.pairs_human_food, r.pairs_human_animal, r.pairs_human_environment, r.clusters_with_a_pair): ok_thr = False; det.append((r.threshold_snps, got_))
check("D.1", "the threshold table is recounted from the pair file at 0, 1, 2, 5, 10, 21, 50 and 100 SNPs", ok_thr and len(thr) == 8, det)
exp = {5: (144, 31), 10: (590, 52), 21: (2248, 75), 50: (5171, 89), 100: (6307, 89)}; got2 = {t_: (int((pr.snp_distance <= t_).sum()), pr[pr.snp_distance <= t_].cluster.nunique()) for t_ in exp}
check("D.2", "144 pairs in 31 clusters within 5 SNPs, 590 in 52 within 10, 2,248 in 75 within 21, 5,171 in 89 within 50 and 6,307 in 89 within 100; 89 clusters in all", got2 == exp and pr.cluster.nunique() == 89, got2)
mins = pr.groupby("cluster").snp_distance.min(); check("D.3", "the closest cross-sector pair of a cluster is a median of 9 SNPs (interquartile range 4 to 17, largest 42)", (mins.median(), mins.quantile(.25), mins.quantile(.75), mins.max()) == (9.0, 4.0, 17.0, 42.0), (mins.median(), mins.quantile(.25), mins.quantile(.75), mins.max()))
bad = 0; checked = 0
for c_, g_ in pr.groupby("cluster"):
    tree = Phylo.read(P(f"snp_trees/{g_.tree.iloc[0]}.nwk"), "newick"); names = {(t_.name or "").strip("'"): t_ for t_ in tree.get_terminals()}; sub = g_.sample(min(len(g_), 25), random_state=1)
    for r in sub.itertuples():
        d_ = tree.distance(names[r.human], names[r.other]); checked += 1
        if abs(d_ - r.snp_distance) > 1e-6: bad += 1
check("D.4", "an independent computation with Bio.Phylo reproduces the distance of up to 25 pairs from each of the 89 clusters (1,191 pairs)", bad == 0 and checked == 1191, (checked, bad))
mm = pd.read_csv(P("table_snp_min_distance.csv")); check("D.5", "all 1,313 frame isolates of the 89 clusters were found in the trees", int(mm.frame_isolates.sum()) == 1313 and int(mm.frame_leaves.sum()) == 1313 and len(mm) == 89, (mm.frame_isolates.sum(), mm.frame_leaves.sum()))
# ================================================================== E. metadata completeness
cmp_ = {"Specimen type (isolation source)": fr.isolation_source.notna(), "Origin of the specimen (US state)": fr.state.notna(), "Phenotypic data (susceptibility result)": fr.has_ast == 1}; mc = pd.read_csv(P("table_metadata_completeness.csv")).set_index("field"); ok_m = True; det = []
for name, mask in cmp_.items():
    for s in SECT:
        v = rnd(100 * mask[fr.sector == s].mean(), 1)
        if v != rnd(mc.loc[name, s], 1): ok_m = False; det.append((name, s, v))
check("E.1", "completeness of specimen type, US state and susceptibility result by sector is recomputed from the frame", ok_m, det)
check("E.2", "human isolates: specimen type 19.3%, US state 1.0%, susceptibility result 0.0%; food 99.6%, 92.8%, 6.9%; animal 99.2%, 85.0%, 2.4%; environment 98.9%, 59.2%, 0.0%",
      [rnd(mc.loc[n, s], 1) for s in SECT for n in cmp_] == ["19.3", "1.0", "0.0", "99.6", "92.8", "6.9", "99.2", "85.0", "2.4", "98.9", "59.2", "0.0"], mc.round(1).to_dict())
# ================================================================== F. the study audit
au = pd.read_csv(P("table_study_audit.csv")); smry = pd.read_csv(P("table_study_audit_summary.csv"))
check("F.1", "35 figures recomputed (Klaas 20, Kitamura 3, McDermott slides 11, Muloi 1): 34 agree and the one that does not is 46 of 207 pairs, printed as approximately 25% and 22.2% from the counts",
      len(au) == 35 and int(au.agrees.sum()) == 34 and dict(zip(smry.study, smry.figures)) == {"Klaas 2026": 20, "Kitamura 2026": 3, "McDermott slides": 11, "Muloi 2023": 1} and au[~au.agrees].recomputed.round(1).tolist() == [22.2], (len(au), int(au.agrees.sum())))
check("F.2", "421 isolates give 88,410 pairs, so 207 close pairs are 0.234% (below the printed 0.25%); 114 + 46 + 47 = 207; the eight sources add up to 1,016 recovered and 421 sequenced isolates",
      math.comb(421, 2) == 88410 and rnd(100 * 207 / 88410, 3) == "0.234" and 114 + 46 + 47 == 207 and sum((109, 219, 207, 152, 76, 179, 30, 44)) == 1016 and sum((60, 60, 60, 60, 61, 60, 30, 30)) == 421, math.comb(421, 2))
check("F.3", "NARMS genomes with phenotype data: the rows add up to 25,703 + 15,382 + 7,698 + 2,097 + 855 = 51,735 and to 10,287 + 21,628 + 17,296 + 2,524 = 51,735; the Japanese isolates add up to 1,613, 682 and 72",
      25703 + 15382 + 7698 + 2097 + 855 == 51735 == 10287 + 21628 + 17296 + 2524 and 623 + 563 + 427 == 1613 and 261 + 57 + 364 == 682 and 38 + 34 == 72)
# ================================================================== G. the chapter text
def docx_text(path):
    from docx import Document
    d = Document(path); body = [p.text for p in d.paragraphs]; tabs = [[[c.text for c in r.cells] for r in t_.rows] for t_ in d.tables]; return d, "\n".join(body), tabs
if A.chapter:
    d, txt, tabs = docx_text(A.chapter); flat = re.sub(r"\s+", " ", txt + " " + " ".join(c for t_ in tabs for r in t_ for c in r))
    figs = re.findall(r"^Figure 9\.(\d+)\.", txt, flags=re.M); check("G.1", "figures are numbered 9.1 to 9.16 in order, one caption for every image, 16 images in all (4 diagrams and 12 figures), 21 tables", figs == [str(i) for i in range(1, 17)] and len(d.inline_shapes) == 16 and len(d.tables) == 21, (figs, len(d.inline_shapes), len(d.tables)))
    objs, on = [], False
    for p in d.paragraphs:
        if p.style.name == "Heading 2": on = (p.text == "Learning objectives")
        elif on and p.text.strip(): objs.append(p.text)
    objs = objs[1:]; check("G.2", "twelve learning objectives, from the structure of acquired genes to the connection with Chapters 4, 7, 8 and 12", len(objs) == 12 and objs[0].startswith("Explain how acquired genes") and objs[-1].startswith("Connect resistance genomics"), len(objs))
    heads = [p.text for p in d.paragraphs if p.style.name == "Heading 2"]; ho = [h for h in heads if h.startswith("Hands-on")]; check("G.3", "eight hands-on sections numbered 1 to 8 in order", [h.split(" ")[1].strip(":") for h in ho] == [str(i) for i in range(1, 9)], ho)
    want = ["261 hits", "23 elements", "30 distinct hits", "17", "475 to 497", "823 false-susceptible", "208 (25.3%)", "28 false-susceptible", "31 of 89", "75 within 21", "35 figures", "34 agree", "22.2%", "19.3%", "207 pairs", "1,313 isolates", "5 of 8 resistant isolates", "44 core rows but 23 distinct elements", "2,824 isolates", "51,735"]
    miss = [w for w in want if w not in flat]; check("G.4", "the sentences quote the validated numbers", not miss, miss)
    check("G.5", "no stray backslash, image path, header artefact or placeholder in the text", not re.search(r"\\\\|media/|@Emile|TODO|XXX|Nov 9, 2026", flat), re.findall(r"\\\\|media/|@Emile|TODO|XXX|Nov 9, 2026", flat)[:3])
    def tbl(prefix): return next((t_ for t_ in tabs if t_[0][0].startswith(prefix)), None)
    t_thr = tbl("SNP threshold"); ok = t_thr is not None and all(row[1].replace(",", "") == str(int(thr.pairs[i])) and row[5] == str(int(thr.clusters_with_a_pair[i])) for i, row in enumerate(t_thr[1:]))
    check("G.6", "the SNP-threshold table of the chapter equals the result file row by row", ok and len(t_thr) == 9, None if t_thr is None else t_thr[1][:3])
    t_comp = tbl("Field of the guidance"); ok = t_comp is not None and all(abs(float(row[j + 1]) - float(mc.loc[mc.index[[k for k in range(len(mc)) if mc.index[k] == row[0]][0]], s])) < 0.05 for row in t_comp[1:] for j, s in enumerate(SECT)); check("G.7", "the completeness table of the chapter equals the result file", ok)
    t_dis = tbl("Drug"); dd = pd.read_csv(P("table_discordance.csv")).set_index("drug"); ok = t_dis is not None and all(int(row[1]) == dd.loc[row[0], "fn"] and int(row[2]) == dd.loc[row[0], "fn_no_core_determinant"] and int(row[5]) == dd.loc[row[0], "fp"] for row in t_dis[1:] if row[0] in dd.index)
    check("G.8", "the table of false-susceptible and false-resistant calls equals the result file", ok)
    t_au = tbl("Source"); ok = t_au is not None and t_au[-1][1:] == ["35", "34", "1"]; check("G.9", "the audit table of the chapter ends with 35 figures, 34 agreeing, 1 not", ok, None if t_au is None else t_au[-1])
    import zipfile; z = zipfile.ZipFile(A.chapter); rels = z.read("word/_rels/document.xml.rels").decode(); dois = {u.rstrip(".,;") for u in re.findall(r"https://doi\.org/[^\"'<> ]+", rels)}
    dois = {x.replace("&amp;", "&") for x in dois}; check("G.10", f"the chapter links {EXPECTED_DOIS} distinct DOIs", len(dois) == EXPECTED_DOIS, len(dois))
    if A.doi:
        import urllib.request, urllib.parse, json; bad_d = []
        for dd_ in sorted(dois):
            doi = dd_.split("doi.org/")[1]
            try: urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "mepi-book/1.0 (mailto:none@example.org)"}), timeout=40).read()
            except Exception as e: bad_d.append((doi, str(e)[:30]))
        check("G.11", "every DOI of the chapter resolves in Crossref", not bad_d, bad_d)
    else: skip("G.11", "every DOI of the chapter resolves in Crossref", "run with --doi")
else:
    for k, w in [("G.1", "figure numbers, tables"), ("G.2", "objectives"), ("G.3", "hands-on sections"), ("G.4", "sentences"), ("G.5", "stray characters"), ("G.6", "SNP table"), ("G.7", "completeness table"), ("G.8", "discordance table"), ("G.9", "audit table"), ("G.10", "DOI count"), ("G.11", "DOI resolution")]: skip(k, w, "run with --chapter")
# ================================================================== H. the figures and scripts
figs_ = ["diagram_pipeline_at_a_glance", "diagram_resistance_mechanisms", "diagram_one_health_sectors", "diagram_integrated_surveillance", "fig1_sampling_frame", "fig2_serovar_mix", "fig3_metadata_completeness", "fig4_class_prevalence", "fig5_mdr_forest_trends", "fig5b_r_forest", "fig6_genotype_phenotype", "fig7_discordance", "fig8_mixed_clusters", "fig9_cooccurrence", "fig10_two_tools", "fig11_snp_thresholds", "fig12_agentic_workflow_diagram"]
fd = next((d_ for d_ in ("../figures", "figures", os.path.join(HERE, "../figures")) if os.path.isdir(d_)), None)
if fd: miss = [f_ for f_ in figs_ if not os.path.exists(os.path.join(fd, f_ + ".png"))]; check("H.1", "all 17 figure files are present (4 diagrams, 12 figures and the R forest plot)", not miss, miss)
else: skip("H.1", "figure files", "no figures folder")
cd = next((d_ for d_ in (".", HERE, os.path.join(HERE, "../code")) if os.path.exists(os.path.join(d_, "ex9_8_study_audit.py"))), None)
if cd:
    import py_compile; badc = []
    for f_ in sorted(glob.glob(os.path.join(cd, "*.py"))):
        try: py_compile.compile(f_, doraise=True)
        except Exception as e: badc.append(os.path.basename(f_))
    check("H.2", "every Python script of the package compiles", not badc, badc)
if A.full:
    import statsmodels.api as sm, statsmodels.formula.api as smf
    fd_ = pd.read_csv(P("frame_for_r.csv.gz")); fd_["sector"] = pd.Categorical(fd_.sector, ["Human", "Food", "Animal", "Environment"]); m_ = smf.glm("mdr_geno ~ sector + C(serovar_grp) + C(year)", data=fd_, family=sm.families.Binomial()).fit(); orr = np.exp(m_.params)
    check("I.glm", "a fresh logistic regression with statsmodels gives adjusted odds ratios of 3.220 (food), 2.577 (animal) and 0.762 (environment)", [rnd(orr["sector[T.Food]"], 3), rnd(orr["sector[T.Animal]"], 3), rnd(orr["sector[T.Environment]"], 3)] == ["3.220", "2.577", "0.762"], orr.round(3).to_dict())
    shutil.rmtree(tmp, ignore_errors=True)
else: skip("I.*", "reruns of the scripts and the model", "run with --full")
out = pd.DataFrame(R); out.to_csv("validation_results_ch9.csv", index=False)
for r in R: print(f"{r['result']:5} {r['id']:8} {r['check'][:150]}" + (f"   [{r['detail']}]" if r["result"] != "PASS" else ""))
print(f"\n{(out.result == 'PASS').sum()} PASS, {(out.result == 'FAIL').sum()} FAIL, {(out.result == 'SKIP').sum()} SKIP"); sys.exit(1 if (out.result == "FAIL").any() else 0)
