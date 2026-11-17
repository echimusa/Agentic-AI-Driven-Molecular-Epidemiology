"""Loader for the surveillance table: reads sal_surv.pkl if present, otherwise the portable sal_surv.csv.gz shipped with the chapter."""
import os
import pandas as pd
def load_surv():
    if os.path.exists("sal_surv.pkl"): return pd.read_pickle("sal_surv.pkl")
    m = pd.read_csv("sal_surv.csv.gz", dtype={"country": str, "year": float}, parse_dates=["coll_day", "coll_month", "release"], keep_default_na=True)
    m["has_coords"] = m.has_coords.astype(bool); return m
