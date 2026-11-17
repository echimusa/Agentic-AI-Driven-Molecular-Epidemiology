"""Loader for the Salmonella assembly-quality table used in Hands-on 1 and 4: reads salmonella_assembly_qc.csv.gz if present, otherwise builds it from asm_batch.pkl (Chapter 2) and typing_all.csv.gz (Chapter 1)."""
import os
import pandas as pd
def load_asm():
    if os.path.exists("salmonella_assembly_qc.csv.gz"): return pd.read_csv("salmonella_assembly_qc.csv.gz", dtype={"country": str, "lab": str, "serovar": str})
    a = pd.read_pickle("asm_batch.pkl"); t = pd.read_csv("typing_all.csv.gz", usecols=["target_acc", "snp_cluster", "year"], low_memory=False)
    d = a.drop(columns=["year"], errors="ignore").merge(t, on="target_acc"); d["lab"] = d.sequenced_by.fillna(d.bioproject_center); d["in_cluster"] = d.snp_cluster.notna().astype(int)
    return d[["target_acc", "plat", "asm_stats_length_bp", "asm_stats_n_contig", "asm_stats_contig_n50", "number_amr_genes", "lab", "serovar", "year", "sector", "country", "in_cluster"]]
if __name__ == "__main__":
    d = load_asm(); d.to_csv("salmonella_assembly_qc.csv.gz", index=False); print(len(d), "rows written")
