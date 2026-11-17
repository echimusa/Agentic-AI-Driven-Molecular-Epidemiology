"""
Chapter 12, Hands-on 7 - What free text costs: the vocabulary of the public Salmonella table.
For the fields that carry place, host and source, count how many distinct strings the 911,210 records use, how much a simple clean-up (lower case, spaces trimmed) removes, how concentrated the use is, and how many records
the Chapter 9 sector rules cannot classify although text is present. A controlled vocabulary would turn these long tails into a short list of terms.
Input: sal_surv.csv.gz. Output: table_standards.csv, table_placeholders.csv, table_unrecognised.csv, figs12/fig7_standards.png.
"""
import re
import numpy as np, pandas as pd
from surv12 import load_surv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
m = load_surv(); PLACE = re.compile(r"^(not provided|not collected|missing|unknown|not applicable|none|n/?a|restricted access|other)$", re.I)
rows = []; fields = {"geo_loc_name": "Place (geo_loc_name)", "host": "Host", "isolation_source": "Isolation source", "serovar_clean": "Serovar (cleaned)"}
for c, lab in fields.items():
    s = m[c].dropna().astype(str); vc = s.value_counts(); norm = s.str.lower().str.strip().str.replace(r"\s+", " ", regex=True)
    rows.append(dict(field=lab, records=len(s), distinct_raw=len(vc), distinct_lower_trimmed=norm.nunique(), top10_share_pct=100 * vc.head(10).sum() / len(s), top100_share_pct=100 * vc.head(100).sum() / len(s), used_once_pct_of_distinct=100 * (vc == 1).sum() / len(vc)))
cty = m.geo_loc_name.dropna().astype(str).str.split(":").str[0].str.strip(); rows.append(dict(field="Place reduced to the text before the colon", records=len(cty), distinct_raw=cty.nunique(), distinct_lower_trimmed=cty.str.lower().nunique(), top10_share_pct=100 * cty.value_counts().head(10).sum() / len(cty), top100_share_pct=100 * cty.value_counts().head(100).sum() / len(cty), used_once_pct_of_distinct=100 * (cty.value_counts() == 1).sum() / cty.nunique()))
T = pd.DataFrame(rows).round(1); T.to_csv("table_standards.csv", index=False); print(T.to_string(index=False))
ph = []
for c, lab in (("geo_loc_name", "Place"), ("host", "Host"), ("isolation_source", "Isolation source")):
    s = m[c].dropna().astype(str).str.split(":").str[0].str.strip() if c == "geo_loc_name" else m[c].dropna().astype(str).str.strip()
    for term, n in s[s.str.match(PLACE)].str.lower().value_counts().items(): ph.append(dict(field=lab, term=term, records=int(n), pct_of_field=round(100 * n / len(s), 2)))
P = pd.DataFrame(ph).sort_values(["field", "records"], ascending=[True, False]); P.to_csv("table_placeholders.csv", index=False); print(P.head(12).to_string(index=False))
hs = m.host.dropna().astype(str); hum = hs[hs.str.contains(r"homo|human|patient", case=False)]; print(f"host strings that name a human: {hum.nunique()} spellings over {len(hum):,} records; the commonest is {hum.value_counts().index[0]!r} ({hum.value_counts().iloc[0] / len(hum):.2%})")
u = m[(m.sector == "Unclassified")]; has_text = u.host.notna() | u.isolation_source.notna(); print(f"unclassified records: {len(u):,}; with host or source text: {has_text.sum():,} ({has_text.sum() / len(m):.1%} of all records); with neither: {(~has_text).sum():,} ({(~has_text).sum() / len(m):.1%})")
txt = u[has_text].assign(t=lambda d: d.isolation_source.fillna(d.host).astype(str).str.lower().str.strip()); U = txt.t.value_counts().head(12).rename_axis("text").reset_index(name="records"); U["pct_of_unclassified_with_text"] = (100 * U.records / has_text.sum()).round(1); U.to_csv("table_unrecognised.csv", index=False); print(U.to_string(index=False))
sec = pd.Series({"classified into a sector": int((m.sector != "Unclassified").sum()), "text present but not recognised": int(has_text.sum()), "no host or source text": int((~has_text).sum())}); print((100 * sec / len(m)).round(1).to_dict())
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.1, 1.0, 0.9]})
f = m.isolation_source.dropna().astype(str).value_counts().values; ax[0].loglog(np.arange(1, len(f) + 1), f, color="#2b6cb0", lw=2); ax[0].set(xlabel="Rank of the isolation-source string", ylabel="Records", title="A. A long tail: 9,248 strings, half used once")
ax[0].axvline(100, color="#c53030", ls="--", lw=1); ax[0].text(110, f.max() / 3, f"top 100 cover {T.loc[2, 'top100_share_pct']:.0f}%", color="#c53030", fontsize=9)
x = np.arange(4); ax[1].bar(x - 0.2, T.distinct_raw[:4], 0.4, color="#a0aec0", label="distinct strings"); ax[1].bar(x + 0.2, T.distinct_lower_trimmed[:4], 0.4, color="#2b6cb0", label="after lower case and trimming")
ax[1].set(yscale="log", ylim=(100, 20000), xticks=x, xticklabels=["Place", "Host", "Source", "Serovar"], ylabel="Distinct values (log scale)", title="B. What a simple clean-up removes"); ax[1].legend(frameon=False, fontsize=8.5)
ax[2].bar(range(3), 100 * sec.values / len(m), color=["#2b6cb0", "#d69e2e", "#a0aec0"]); ax[2].set(xticks=range(3), xticklabels=["classified", "text not\nrecognised", "no text"], ylabel="% of 911,210 records", title="C. Sector rules (Chapter 9)")
for i, v in enumerate(100 * sec.values / len(m)): ax[2].text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=9)
fig.tight_layout(); fig.savefig("figs12/fig7_standards.png", dpi=200); plt.close(fig); print("figure written")
