"""
Chapter 8, Hands-on 3 (step 4) - prepare the time tree and a metadata table for Microreact (https://microreact.org).
Microreact takes a Newick tree and a table whose first column (id) matches the tip names; latitude and longitude (decimal degrees) place each genome on the map and
year, month and day columns drive the timeline. Locations are the country centroids (rounded to 0.1 degree), never the place of a patient.
Writes microreact_tree.nwk (591 tips, branch lengths in years) and microreact_metadata.csv, and checks that every tip has a row, a coordinate pair and a valid date.
"""
import json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from Bio import Phylo
from shapely.geometry import shape
import phylo8
st = phylo8.load_states(); out = set(pd.read_csv("zika_clock_outliers.csv").acc); tree = Phylo.read("zika_timetree.nwk", "newick")
for acc in out: tree.prune(acc)
tips = [t.name for t in tree.get_terminals()]; Phylo.write(tree, "microreact_tree.nwk", "newick")
geo = json.load(open("ne_50m_countries.geojson"))["features"]; polys = {f["properties"]["NAME"]: shape(f["geometry"]) for f in geo}
NAMES = {"USA": "United States of America", "Dominican Republic": "Dominican Rep.", "French Polynesia": "Fr. Polynesia", "Viet Nam": "Vietnam", "Virgin Islands": "U.S. Virgin Is.", "Cape Verde": "Cabo Verde"}
MANUAL = {"Martinique": (14.6, -61.0), "French Guiana": (3.9, -53.1)}          # territories that the Natural Earth country layer does not list; approximate centres typed by hand
def centroid(c):
    if c in MANUAL: return MANUAL[c]
    g = polys[NAMES.get(c, c)]; parts = list(g.geoms) if g.geom_type == "MultiPolygon" else [g]; big = max(parts, key=lambda p: p.area); return round(big.centroid.y, 1), round(big.centroid.x, 1)
meta = pd.read_csv("zika_selected_meta.csv", dtype=str).set_index("acc").loc[tips]; cen = {c: centroid(c) for c in meta.cty.unique()}
M = pd.DataFrame({"id": tips, "country": meta.cty.values, "place_group": st.loc[tips].state.values, "latitude": [cen[c][0] for c in meta.cty], "longitude": [cen[c][1] for c in meta.cty],
                  "year": meta.year.astype(float).astype(int).values, "month": meta.month.astype(float).astype(int).values, "day": [int(float(d)) if p == "day" else "" for d, p in zip(meta.day, meta.dprec)], "date_precision": meta.dprec.values})
M.to_csv("microreact_metadata.csv", index=False)
ok = (M.id.is_unique and set(M.id) == set(tips) and M.latitude.between(-90, 90).all() and M.longitude.between(-180, 180).all() and M.month.between(1, 12).all() and M.year.between(2006, 2018).all())
print(f"{len(M)} rows, {len(tips)} tips, {M.country.nunique()} countries, {M[['latitude', 'longitude']].drop_duplicates().shape[0]} distinct map positions; checks passed: {ok}")
