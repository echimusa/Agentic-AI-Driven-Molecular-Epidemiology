"""
Chapter 2 data preparation: builds asm_batch_enteritidis.csv.gz, the assembly metrics, sequencing laboratory, platform, serovar and sector of the Salmonella Enteritidis isolates in the NCBI Pathogen Detection table.
Needs the 883 MB metadata table (salmonella_metadata.tsv) and the parsed surveillance table of Chapter 12 (sal_surv.csv.gz). The shipped file is the output of this script.
"""
import csv
import numpy as np, pandas as pd
cols = ["target_acc", "Platform", "asm_stats_n_contig", "asm_stats_contig_n50", "asm_stats_length_bp", "number_amr_genes", "sequenced_by", "serovar", "bioproject_center"]
m = pd.read_csv("salmonella_metadata.tsv", sep="\t", na_values=["NULL"], quoting=csv.QUOTE_NONE, usecols=cols, dtype=str, low_memory=False)
s = pd.read_csv("sal_surv.csv.gz", usecols=["target_acc", "sector", "country", "year"], dtype={"country": str})
d = m.merge(s, on="target_acc", how="left")
for c in ["asm_stats_n_contig", "asm_stats_contig_n50", "asm_stats_length_bp", "number_amr_genes"]: d[c] = pd.to_numeric(d[c], errors="coerce")
d["plat"] = d.Platform.str.upper().str.replace("_", " "); print("platforms:", d.plat.value_counts().head(4).to_dict())
d[d.serovar.str.lower() == "enteritidis"].to_csv("asm_batch_enteritidis.csv.gz", index=False); print("Enteritidis isolates written:", int((d.serovar.str.lower() == "enteritidis").sum()))
