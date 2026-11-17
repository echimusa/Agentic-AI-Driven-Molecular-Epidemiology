"""
Chapter 12 data preparation: reads the surveillance-relevant fields of the NCBI Pathogen Detection Salmonella table (the file downloaded for Chapter 9),
assigns One Health sectors with the Chapter 9 rules, parses collection and release dates, and writes sal_surv.pkl.
"""
import csv, re, os
import numpy as np, pandas as pd
from amr9 import add_sector
SRC = "../ch9/salmonella_metadata.tsv" if os.path.exists("../ch9/salmonella_metadata.tsv") else "salmonella_metadata.tsv"
cols = ["target_acc", "collection_date", "sra_release_date", "target_creation_date", "geo_loc_name", "lat_lon", "outbreak", "epi_type", "host", "isolation_source", "IFSAC_category", "source_type", "serovar", "computed_types", "bioproject_center", "sequenced_by", "collected_by"]
m = pd.read_csv(SRC, sep="\t", na_values=["NULL"], quoting=csv.QUOTE_NONE, usecols=cols, dtype=str, low_memory=False)
m["epi_type"] = m.epi_type.astype(object); m = add_sector(m)            # sector, country (text before the colon), year, serovar_clean
m["date_precision"] = np.select([m.collection_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}", na=False), m.collection_date.str.fullmatch(r"\d{4}-\d{2}", na=False), m.collection_date.str.fullmatch(r"\d{4}", na=False)], ["day", "month", "year"], "missing")
cd = pd.to_datetime(m.collection_date.where(m.date_precision == "day"), errors="coerce"); m["coll_day"] = cd
m["coll_month"] = pd.to_datetime(m.collection_date.where(m.date_precision.isin(["day", "month"])).str[:7] + "-01", errors="coerce")
m["release"] = pd.to_datetime(m.sra_release_date, errors="coerce"); m["creation"] = pd.to_datetime(m.target_creation_date, errors="coerce")
m["lag_days"] = (m.release - m.coll_day).dt.days
m["has_coords"] = m.lat_lon.notna() & ~m.lat_lon.str.lower().isin(["not collected", "missing", "not applicable", "not provided", "na", "unknown"]).fillna(False)
keep = ["target_acc", "sector", "country", "geo_loc_name", "year", "serovar_clean", "bioproject_center", "date_precision", "coll_day", "coll_month", "release", "creation", "lag_days", "has_coords", "outbreak", "host", "isolation_source", "source_type", "epi_type"]
m[keep].to_pickle("sal_surv.pkl"); print(f"{len(m):,} isolates; sectors {m.sector.value_counts().to_dict()}")
print("date precision (all):", m.date_precision.value_counts(normalize=True).round(3).to_dict()); print("date precision (human):", m[m.sector == "Human"].date_precision.value_counts(normalize=True).round(3).to_dict())
print("release before collection (day-level records):", int((m.lag_days < 0).sum()), "of", int(m.lag_days.notna().sum()))
