"""
Chapter 9, Hands-on 1 (extension) - how complete are the public metadata against the core data fields of the Quadripartite guidance?
The Quadripartite guidance (WHO, FAO, UNEP and WOAH, 2025; Table 4) lists the core data an integrated surveillance system should report for every isolate: an isolate identifier, a specimen identifier, the
specimen type, the collection date, the origin of the specimen, the identification of the isolate (species, serotype), the phenotypic result, and, where genomes are used, the resistance gene profile and a BioProject.
This script counts, for each One Health sector of the frame of ex9_1 (United States, 2015-2024), the share of isolates that carry each field in the NCBI Pathogen Detection record.
The collection year is complete by construction, because the frame keeps only isolates with a year; the resistance gene profile is not counted, because an isolate with no determinant has an empty genotype and not a missing one.
Input: sal_us_frame.csv.gz.  Output: table_metadata_completeness.csv.
"""
import numpy as np, pandas as pd
f = pd.read_csv("sal_us_frame.csv.gz", dtype={"state": str}); order = ["Human", "Food", "Animal", "Environment"]
fields = [("Isolate identifier (target accession)", f.target_acc.notna()), ("Specimen identifier (BioSample)", f.biosample_acc.notna()), ("Specimen type (isolation source)", f.isolation_source.notna()),
          ("Collection date (year)", f.year.notna()), ("Origin of the specimen (US state)", f.state.notna()), ("Identification (serovar)", f.serovar_clean.notna() & (f.serovar_clean.astype(str) != "")),
          ("Phenotypic data (susceptibility result)", f.has_ast == 1), ("Project of origin (BioProject centre)", f.bioproject_center.notna())]
rows = []
for name, mask in fields:
    r = {"field": name}; r.update({s: round(100 * float(mask[f.sector == s].mean()), 1) for s in order}); r["all"] = round(100 * float(mask.mean()), 1); rows.append(r)
T = pd.DataFrame(rows); T.to_csv("table_metadata_completeness.csv", index=False); print(T.to_string(index=False)); print({s: int((f.sector == s).sum()) for s in order}, len(f))
