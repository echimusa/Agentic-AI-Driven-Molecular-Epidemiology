"""
Chapter 5 data preparation: 60 public influenza A (H3N2) hemagglutinin (HA, segment 4) complete coding sequences from NCBI, 12 per collection year 2014-2018, human isolates only
(strain names of the form place/identifier/year; names with a host such as swine/ or a duck species are excluded). The first 6,000 of the matching records returned by NCBI are used as the pool; the draw is seeded.
Output: h3n2_ha_raw.fasta, h3n2_meta.csv; then align with MAFFT and fit the tree with IQ-TREE (see run commands in the README).
"""
import json, re, time, random, urllib.request, urllib.parse, collections
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for a in range(5):
        try: return urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "molecular-epi-book"}), timeout=60).read().decode()
        except Exception: time.sleep(3 * (a + 1))
    raise RuntimeError(u)
q = '"Influenza A virus"[Organism] AND H3N2[All Fields] AND hemagglutinin[All Fields] AND segment 4[All Fields] AND complete cds[All Fields] AND 1650:1800[SLEN] AND ("2014/01/01"[PDAT]:"2019/12/31"[PDAT])'
r = json.loads(get(E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "term": q, "retmax": 6000, "retmode": "json"}))); ids = r["esearchresult"]["idlist"]; print("hits:", r["esearchresult"]["count"], "retrieved ids:", len(ids))
meta = {}; nonhuman = 0
for i in range(0, len(ids), 300):
    s = json.loads(get(E + "esummary.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids[i:i + 300]), "retmode": "json"})))["result"]
    for k in s["uids"]:
        t = s[k]["title"]; m = re.search(r"\(A/([^)]*?/(\d{4}))\(H3N2\)\)", t)
        if not m: continue
        if m.group(1).count("/") != 2: nonhuman += 1; continue                                              # host-prefixed names (swine/..., duck species) have an extra field
        meta[s[k]["accessionversion"]] = (m.group(1), int(m.group(2)), t)
    time.sleep(0.4)
print("human-style strain names:", len(meta), "| host-prefixed names excluded:", nonhuman); print(sorted(collections.Counter(v[1] for v in meta.values()).items()))
random.seed(5); pick = []
for y in range(2014, 2019):
    ys = sorted(a for a, v in meta.items() if v[1] == y); pick += random.sample(ys, min(12, len(ys)))
fa = get(E + "efetch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "id": ",".join(pick), "rettype": "fasta", "retmode": "text"})); open("h3n2_ha_raw.fasta", "w").write(fa); print("selected:", len(pick), "| fasta records:", fa.count(">"))
with open("h3n2_meta.csv", "w") as f:
    f.write("accession,strain,year\n"); [f.write(f"{a},\"{meta[a][0]}\",{meta[a][1]}\n") for a in pick]
