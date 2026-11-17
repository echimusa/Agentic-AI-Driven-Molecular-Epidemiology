"""
Chapter 1 data preparation: builds typing_all.csv.gz (serotype, PFGE patterns, SNP cluster and sector for every isolate) from the NCBI Pathogen Detection Salmonella metadata table,
its SNP cluster list and the parsed surveillance table of Chapter 12 (sal_surv.csv.gz). Needs the 883 MB metadata table. The shipped typing_pfge_subset.csv.gz is the part of this table
that carries a PFGE pattern, which is all that ex1_1_discrimination.py uses.
"""
import csv
import pandas as pd
cols = ["target_acc", "computed_types", "PFGE_PrimaryEnzyme_pattern", "PFGE_SecondaryEnzyme_pattern", "serovar", "minsame", "mindiff"]
m = pd.read_csv("salmonella_metadata.tsv", sep="\t", na_values=["NULL"], quoting=csv.QUOTE_NONE, usecols=cols, dtype=str, low_memory=False)
m["serotype"] = m.computed_types.str.extract(r'serotype=([^",]+)')[0]
s = pd.read_csv("sal_surv.csv.gz", usecols=["target_acc", "sector", "country", "year", "date_precision"], dtype={"country": str})
cl = pd.read_csv("salmonella_cluster_list.tsv", sep="\t", usecols=["PDS_acc", "target_acc"])
d = m.drop(columns="computed_types").merge(s, on="target_acc", how="left").merge(cl, on="target_acc", how="left").rename(columns={"PFGE_PrimaryEnzyme_pattern": "pfge_xbai", "PFGE_SecondaryEnzyme_pattern": "pfge_blni", "PDS_acc": "snp_cluster"})
d.to_csv("typing_all.csv.gz", index=False); d[d.pfge_xbai.notna()].to_csv("typing_pfge_subset.csv.gz", index=False); print(f"{len(d):,} isolates; {d.pfge_xbai.notna().sum():,} with a PFGE XbaI pattern")
