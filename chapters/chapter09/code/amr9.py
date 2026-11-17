"""
Chapter 9 helpers: One Health sectors, resistance genotypes, drug classes, genotype-to-phenotype rules, confidence intervals.
Input: the NCBI Pathogen Detection Salmonella metadata table, the AMRFinderPlus reference gene catalogue and point-mutation table.
"""
import re, csv, os
import numpy as np, pandas as pd

KEEP = ["target_acc", "asm_acc", "biosample_acc", "bioproject_acc", "bioproject_center", "collected_by", "collection_date", "epi_type", "geo_loc_name", "host",
        "isolation_source", "IFSAC_category", "source_type", "scientific_name", "serovar", "computed_types", "AST_phenotypes", "AMR_genotypes", "AMR_genotypes_core",
        "number_drugs_resistant", "number_drugs_tested", "number_amr_genes", "minsame", "mindiff", "amrfinder_version", "refgene_db_version"]

def load_metadata(path="salmonella_metadata.tsv", cache="sal_meta_raw.pkl", ast_subset="sal_ast_subset.csv.gz"):
    if os.path.exists(cache): m = pd.read_pickle(cache)
    elif not os.path.exists(path) and os.path.exists(ast_subset):          # the small file shipped with the chapter: only isolates with susceptibility results
        m = pd.read_csv(ast_subset, low_memory=False)
    else:
        m = pd.read_csv(path, sep="\t", na_values=["NULL"], low_memory=False, quoting=csv.QUOTE_NONE, usecols=lambda c: c.lstrip("#") in KEEP); m.to_pickle(cache)
    m.columns = [c.lstrip("#") for c in m.columns]
    return m

# ---------------------------------------------------------------- One Health sector of each isolate
ANIMAL_HOST = re.compile(r"bos taurus|bovine|cattle|\bcow\b|calf|gallus|chicken|sus scrofa|swine|porcine|\bpig\b|turkey|meleagris|avian|equus|horse|ovis|sheep|capra|goat|canis|\bdog\b|felis|\bcat\b|duck|anas |bird|reptile|snake|turtle|lizard", re.I)
FOOD_SRC = re.compile(r"chicken|beef|pork|turkey|meat|poultry|sausage|egg|vegetable|lettuce|spinach|tomato|cucumber|cantaloupe|melon|sprout|nut|almond|pepper|herb|spice|flour|food|fish|shrimp|seafood|cheese|milk|salad|fruit|ground", re.I)
ENV_SRC = re.compile(r"water|soil|sewage|wastewater|swab|surface|drain|environment|sediment|manure|litter|dust|pond|river|stream|creek|farm", re.I)

def sector_of(row):
    ifs = str(row.IFSAC_category) if pd.notna(row.IFSAC_category) else ""
    st = str(row.source_type).lower() if pd.notna(row.source_type) else ""
    host = str(row.host) if pd.notna(row.host) else ""
    src = str(row.isolation_source) if pd.notna(row.isolation_source) else ""
    parts = [p.strip() for p in ifs.split("|")] if ifs else []
    if ifs.startswith("clinical/research| human") or st == "human" or re.search(r"homo sapiens|^human$", host, re.I): return "Human"
    if any(p.startswith("veterinary") for p in parts) or st == "animal" or "animal feed" in ifs: return "Animal"
    if st == "food" or any(p and not p.startswith(("environmental", "clinical", "veterinary", "animal feed", "other")) for p in parts): return "Food"
    if st == "environmental" or any(p.startswith("environmental") for p in parts): return "Environment"
    if ANIMAL_HOST.search(host) and row.epi_type != "clinical": return "Animal"
    if row.epi_type == "environmental/other" and FOOD_SRC.search(src) and not ANIMAL_HOST.search(host): return "Food"
    if row.epi_type == "environmental/other" and ENV_SRC.search(src): return "Environment"
    if ANIMAL_HOST.search(host): return "Animal"
    return "Unclassified"

def norm_serovar(s, ct):
    """'enteritidis', 'Enteritidis' and 'Salmonella enterica subsp. enterica serovar Enteritidis' are the same serovar."""
    if pd.isna(s) and pd.notna(ct):
        mt = re.search(r'serotype=([^"]+)', str(ct)); s = mt.group(1) if mt else np.nan
    if pd.isna(s): return "Unknown"
    s = re.sub(r".*serovar\s+", "", str(s)).strip().replace(";", ",")
    return s if s.startswith(("I ", "II ", "IIIa", "IIIb", "IV ")) or not s else s[0].upper() + s[1:]

def add_sector(m):
    m = m.copy(); m["sector"] = m.apply(sector_of, axis=1)
    m["country"] = m.geo_loc_name.fillna("").str.split(":").str[0].str.strip().replace("", np.nan)
    m["serovar_clean"] = [norm_serovar(a, b) for a, b in zip(m.serovar, m.computed_types)]
    m["year"] = pd.to_numeric(m.collection_date.astype(str).str.extract(r"^(\d{4})")[0], errors="coerce")
    return m

# ---------------------------------------------------------------- genotypes and drug classes
GOOD_SUFFIX = {"", "COMPLETE", "POINT", "ALLELE", "HMM", "EXACT", "BLAST", "INTERNAL_STOP"}
def parse_genotype(s, keep_partial=False):
    """'aph(6)-Id,"gyrA_D87Y=POINT","fosA7=PARTIAL"' -> list of element symbols called complete or as point mutations."""
    if pd.isna(s): return []
    out = []
    for tok in str(s).replace('"', "").split(","):
        tok = tok.strip()
        if not tok: continue
        name, _, suf = tok.partition("=")
        if suf in ("COMPLETE", "POINT", "ALLELE", "HMM", "EXACT", "BLAST", "") or keep_partial: out.append(name)
    return out

def load_catalog(db=None):
    db = db or next(d for d in ("amrf/db", "amrf_db_catalog") if os.path.exists(d))
    c = pd.read_csv(f"{db}/ReferenceGeneCatalog.txt", sep="\t", dtype=str)
    c = c[c.type == "AMR"].rename(columns={"class": "cls"}); cat = {}
    for col in ("gene_family", "allele"):
        for r in c.dropna(subset=[col]).itertuples(): cat[getattr(r, col)] = (r.cls, r.subclass if isinstance(r.subclass, str) else "")
    mu = pd.read_csv(f"{db}/AMRProt-mutation.tsv", sep="\t", dtype=str).rename(columns={"class": "cls", "#taxgroup": "taxgroup"}); mu = mu[mu.taxgroup == "Salmonella"]
    for r in mu.itertuples(): cat[r.standard_mutation_symbol] = (r.cls, r.subclass if isinstance(r.subclass, str) else "")
    return cat

def classes_of(genes, cat):
    return sorted({cat[g][0] for g in genes if g in cat})

# genotype-to-phenotype rules: which determinants predict resistance to which drug (Salmonella, NARMS panel)
def has(genes, cat, cls=None, sub=None, point_only=False, any_of=None):
    for g in genes:
        if g not in cat: continue
        c, s = cat[g]
        if point_only and "_" not in g: continue
        if cls and c != cls: continue
        if sub and sub not in s: continue
        return True
    return False

RULES = {
 "ampicillin":                      lambda g, c: has(g, c, "BETA-LACTAM"),
 "amoxicillin-clavulanic acid":     lambda g, c: has(g, c, "BETA-LACTAM", "CEPHALOSPORIN"),
 "cefoxitin":                       lambda g, c: has(g, c, "BETA-LACTAM", "CEPHALOSPORIN"),
 "ceftiofur":                       lambda g, c: has(g, c, "BETA-LACTAM", "CEPHALOSPORIN"),
 "ceftriaxone":                     lambda g, c: has(g, c, "BETA-LACTAM", "CEPHALOSPORIN"),
 "streptomycin":                    lambda g, c: has(g, c, "AMINOGLYCOSIDE", "STREPTOMYCIN"),
 "gentamicin":                      lambda g, c: has(g, c, "AMINOGLYCOSIDE", "GENTAMICIN"),
 "kanamycin":                       lambda g, c: has(g, c, "AMINOGLYCOSIDE", "KANAMYCIN"),
 "amikacin":                        lambda g, c: has(g, c, "AMINOGLYCOSIDE", "AMIKACIN"),
 "tetracycline":                    lambda g, c: has(g, c, "TETRACYCLINE"),
 "chloramphenicol":                 lambda g, c: has(g, c, "PHENICOL"),
 "sulfamethoxazole":                lambda g, c: has(g, c, "SULFONAMIDE"),
 "trimethoprim-sulfamethoxazole":   lambda g, c: has(g, c, "TRIMETHOPRIM"),
 "nalidixic acid":                  lambda g, c: has(g, c, "QUINOLONE", point_only=True),
 "ciprofloxacin":                   lambda g, c: has(g, c, "QUINOLONE"),
 "azithromycin":                    lambda g, c: has(g, c, "MACROLIDE"),
 "sulfisoxazole":                   lambda g, c: has(g, c, "SULFONAMIDE"),
 "meropenem":                       lambda g, c: has(g, c, "BETA-LACTAM", "CARBAPENEM"),
 "colistin":                        lambda g, c: has(g, c, "COLISTIN"),
}

def parse_ast(s):
    if pd.isna(s): return {}
    d = {}
    for tok in str(s).replace('"', "").split(","):
        if "=" in tok: k, v = tok.rsplit("=", 1); d[k.strip()] = v.strip()
    return d

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)
