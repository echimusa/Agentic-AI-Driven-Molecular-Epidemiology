"""
Chapter 10, Hands-on 1 - What have GWAS of infection found, and in whom?
Uses the public NHGRI-EBI GWAS Catalog release tables (studies and ancestries). Infection traits are found by keyword on the disease/trait and mapped trait.
Output: table_catalog_groups.csv, table_catalog_ancestry.csv, Figure 1.
"""
import os, re, csv, urllib.request
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs10", exist_ok=True)
B = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest"
for f, u in {"gwas_studies_latest.tsv": "gwas-catalog-download-studies-v1.0.3.1.tsv", "gwas_ancestry_latest.tsv": "gwas-catalog-download-ancestries-v1.0.3.1.tsv"}.items():
    if not os.path.exists(f): urllib.request.urlretrieve(f"{B}/{u}", f)
rd = lambda f: pd.read_csv(f, sep="\t", low_memory=False, quoting=csv.QUOTE_NONE, on_bad_lines="skip")
s, a = rd("gwas_studies_latest.tsv"), rd("gwas_ancestry_latest.tsv")
print(f"{len(s):,} studies in the release (publication dates to {s['DATE'].max()}, added to {s['DATE ADDED TO CATALOG'].max()}); {len(a):,} ancestry rows")

# ---------------------------------------------------------------- infection traits by keyword (first match wins)
PAT = {"COVID-19": r"covid|sars-cov|coronavirus", "Tuberculosis": r"tubercul", "Malaria": r"malaria|plasmodium", "HIV": r"\bhiv\b|human immunodeficiency",
       "Viral hepatitis": r"hepatitis [abce]\b|viral hepatitis|\bhbv\b|\bhcv\b|hepatitis b|hepatitis c", "Influenza": r"influenza",
       "Other infection": r"sepsis|bacter[a]?emia|pneumonia|staphylococc|streptococc|meningitis|urinary tract infection|clostridi|salmonell|herpes|zoster|varicella|epstein|cytomegalovirus|papillomavirus|chlamydia|gonorrh|syphilis|leishman|dengue|ebola|schistosom|helminth|toxoplasm|polio|infection"}
txt = (s["DISEASE/TRAIT"].fillna("") + " | " + s["MAPPED_TRAIT"].fillna("")).str.lower(); grp = pd.Series("", index=s.index)
for g, p in PAT.items(): grp[(grp == "") & txt.str.contains(p, regex=True)] = g
# the publication-date field has no 2024 entries, so the year each study was added to the catalog is used
s["group"] = grp; s["year"] = pd.to_numeric(s["DATE ADDED TO CATALOG"].astype(str).str[:4], errors="coerce"); s["inf"] = s.group != ""
s["full_ss"] = s["FULL SUMMARY STATISTICS"].astype(str).str.lower().eq("yes")
t = s[s.inf].groupby("group").agg(studies=("STUDY ACCESSION", "size"), associations=("ASSOCIATION COUNT", "sum"), full_summary_stats=("full_ss", "sum")).sort_values("studies", ascending=False)
t["pct_with_summary_stats"] = (100 * t.full_summary_stats / t.studies).round(1); t.loc["All infection"] = [t.studies.sum(), t.associations.sum(), t.full_summary_stats.sum(), round(100 * t.full_summary_stats.sum() / t.studies.sum(), 1)]
t.loc["All traits"] = [len(s), s["ASSOCIATION COUNT"].sum(), s.full_ss.sum(), round(100 * s.full_ss.mean(), 1)]; t.to_csv("table_catalog_groups.csv"); print(t.astype({"studies": int, "associations": int, "full_summary_stats": int}).to_string())

# ---------------------------------------------------------------- ancestry of the discovery samples (initial stage), individuals split equally across listed categories
RULES = [("European", r"European"), ("African", r"African American|Afro-Caribbean|Sub-Saharan|African unspecified"), ("East Asian", r"East Asian"), ("South Asian", r"South Asian"),
         ("Hispanic or Latin American", r"Hispanic or Latin American")]
def cats(x):
    x = str(x); found = [g for g, p in RULES if re.search(p, x)]
    return found if found else ["Other or not reported"]
a = a[a["STAGE"].eq("initial")].copy(); a["n"] = pd.to_numeric(a["NUMBER OF INDIVIDUALS"], errors="coerce").fillna(0); a["cats"] = a["BROAD ANCESTRAL CATEGORY"].map(cats)
a = a.explode("cats").assign(n=lambda d: d.n / d.groupby(level=0).n.transform("size").clip(lower=1) if False else d.n); a["k"] = a.groupby(level=0).cats.transform("size"); a["n"] = a.n / a.k
a = a.merge(s[["STUDY ACCESSION", "group", "inf"]], on="STUDY ACCESSION", how="left"); order = ["European", "African", "East Asian", "South Asian", "Hispanic or Latin American", "Other or not reported"]
def share(d): v = d.groupby("cats").n.sum().reindex(order).fillna(0); return 100 * v / v.sum()
rows = {"All traits": share(a), "All infection": share(a[a.inf == True])}
for g in ["COVID-19", "HIV", "Tuberculosis", "Malaria"]: rows[g] = share(a[a.group == g])
anc = pd.DataFrame(rows).T[order].round(1); anc["individuals"] = [int(a[a.inf == True].n.sum()) if k == "All infection" else int(a[a.group == k].n.sum()) if k in PAT else int(a.n.sum()) for k in anc.index]   # people counted once per study
af = a.assign(is_af=a.cats.eq("African")).groupby("STUDY ACCESSION").is_af.any().rename("has_african").reset_index().merge(s[["STUDY ACCESSION", "group", "inf"]], on="STUDY ACCESSION")
anc["pct_studies_with_african_sample"] = [round(100 * (af.has_african.mean() if k == "All traits" else af[af.inf == True].has_african.mean() if k == "All infection" else af[af.group == k].has_african.mean()), 1) for k in anc.index]
anc = anc.rename(columns={"individuals": "person_study_entries"}); anc.to_csv("table_catalog_ancestry.csv"); print(anc.to_string())

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.5), gridspec_kw={"width_ratios": [1.2, 1.15]})
y = s[s.inf & s.year.between(2008, 2026)].groupby(["year", "group"]).size().unstack(fill_value=0).reindex(columns=list(PAT)).fillna(0); cm = plt.get_cmap("tab10"); bottom = np.zeros(len(y))
for i, g in enumerate(PAT): ax[0].bar(y.index, y[g], bottom=bottom, color=cm(i), label=g, width=0.8); bottom += y[g].values
ax[0].set(title="GWAS of infection traits added to the catalog, by year", ylabel="Studies", xlabel="Year added", xticks=range(2008, 2027, 2)); ax[0].legend(frameon=False, fontsize=8, ncol=1, loc="upper left")
cols = ["#2b6cb0", "#dd6b20", "#2f855a", "#805ad5", "#d69e2e", "#a0aec0"]; left = np.zeros(len(anc))
for c, col in zip(order, cols): ax[1].barh(anc.index[::-1], anc[c][::-1], left=left[::-1] if False else None, color=col) if False else None
left = np.zeros(len(anc)); idx = list(range(len(anc)))[::-1]
for c, col in zip(order, cols): ax[1].barh(idx, anc[c].values, left=left, color=col, label=c, edgecolor="white", linewidth=0.6); left += anc[c].values
ax[1].set_yticks(idx); ax[1].set_yticklabels(anc.index); ax[1].set(xlim=(0, 100), xlabel="Share of discovery-sample individuals (%)", title="Ancestry of the samples")
ax[1].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3); fig.tight_layout(); fig.savefig("figs10/fig1_catalog.png", dpi=200); plt.close(fig); print("figure written")
