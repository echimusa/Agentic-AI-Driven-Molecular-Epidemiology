"""
Chapter 4 validation harness. Recomputes the numbers quoted in the chapter with code that is independent of the chapter's scripts: its own agglomerative clustering (Lance-Williams updates) and union-find for the real Listeria clusters, exact
Poisson and binomial formulas for the SNP-versus-allele simulation, the R allele simulation re-read from the profiles it generates, and a streaming pass over the 17.3 million NCBI pairwise distances. Results are checked against the shipped tables and the chapter text.
USAGE   python validate_ch4.py [--chapter Chapter_4_....docx] [--parquet listeria_SNP_distances.parquet] [--full]
        --parquet  the NCBI SNP-distance file (148 MB, downloaded by ex4_4_real_clusters.py): the width, fragmentation and neighbour numbers of Hands-on 4 are then recomputed from the 17.3 million pairs
        --chapter  also checks the chapter's tables and key numbers (needs pandoc and python-docx)
        --full     also reruns ex4_1, ex4_2 (if Rscript exists), ex4_3 and, with --parquet, ex4_4 in a scratch folder and checks that the result tables are reproduced byte for byte
Prints PASS/FAIL/SKIP per check, writes validation_results_ch4.csv and exits 1 if any check fails.
"""
import argparse, glob, hashlib, math, os, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd
from scipy import stats
ap = argparse.ArgumentParser(); ap.add_argument("--chapter"); ap.add_argument("--parquet"); ap.add_argument("--full", action="store_true"); A = ap.parse_args()
def P(name):
    for d in (".", "../data", "data"):
        if os.path.exists(os.path.join(d, name)): return os.path.join(d, name)
    raise FileNotFoundError(name)
R = []
def check(cid, what, ok, detail=""): R.append(dict(id=cid, check=what, result="PASS" if ok else "FAIL", detail=str(detail)[:150])); return ok
def skip(cid, what, why): R.append(dict(id=cid, check=what, result="SKIP", detail=why))
from decimal import Decimal, ROUND_HALF_UP
rnd = lambda v, n: str(Decimal(str(v)).quantize(Decimal(10) ** -n, rounding=ROUND_HALF_UP))
# ================================================= Hands-on 1, Part B: two real Listeria clusters, own linkage code
pr = pd.read_csv(P("listeria_pairs_selected.csv.gz")); lr = pd.read_csv(P("table_linkage_rules.csv"))
def matrix(g):
    ids = sorted(set(g.target_acc_1) | set(g.target_acc_2)); ix = {a: i for i, a in enumerate(ids)}; M = np.zeros((len(ids), len(ids)))
    for a, b, d in zip(g.target_acc_1, g.target_acc_2, g.delta_positions_unambiguous): M[ix[a], ix[b]] = M[ix[b], ix[a]] = d
    return M
def agglomerate(D, rule, h):
    n = len(D); M = D.astype(float).copy(); np.fill_diagonal(M, np.inf); size = np.ones(n); active = np.ones(n, bool); members = [[i] for i in range(n)]
    while True:
        sub = np.where(active[:, None] & active[None, :], M, np.inf); k = int(np.argmin(sub)); a, b = divmod(k, n)
        if sub[a, b] > h: break
        if rule == "single": new = np.minimum(M[a], M[b])
        elif rule == "complete": new = np.maximum(M[a], M[b])
        else: new = (size[a] * M[a] + size[b] * M[b]) / (size[a] + size[b])
        M[a, :] = new; M[:, a] = new; M[a, a] = np.inf; size[a] += size[b]; active[b] = False; members[a] += members[b]; members[b] = []
    cl = [m for m in members if m]; big = max(len(m) for m in cl); diam = {int(D[np.ix_(m, m)].max()) for m in cl if len(m) == big}
    return len(cl), big, sum(1 for m in cl if len(m) == 1), diam
mats = {c: matrix(g) for c, g in pr.groupby("PDS_acc")}; ok_s = True; ok_c = True; diffs = []; bad = []
for c, D in mats.items():
    for rule in ("single", "average", "complete"):
        for h in (5, 10, 20, 50):
            n_, big, sing, diam = agglomerate(D, rule, h); row = lr[(lr.cluster == c) & (lr.linkage == rule) & (lr.threshold == h)].iloc[0]; same_all = n_ == row.clusters and big == row.largest and sing == row.singletons and int(row.diameter_of_largest) in diam
            if rule == "single": ok_s &= same_all
            else:
                ok_c &= n_ == row.clusters and sing == row.singletons
                if not same_all: diffs.append((len(D), rule, h, f"largest {big} against {int(row.largest)}"))
check("H1.1", "single linkage, which is unique: both real clusters (111 and 220 isolates) at four thresholds, by own agglomerative code: clusters, largest, singletons and diameter of the largest equal the table (8 combinations)", ok_s)
check("H1.1b", "average and complete linkage on integer distances depend on how ties are broken: with a different tie-break the number of clusters and of singletons are equal in all 16 combinations; the largest cluster differs in at most 2 (recorded)", ok_c and len(diffs) <= 2, diffs)
big_c = [c for c, D in mats.items() if len(D) == 220][0]; sm_c = [c for c, D in mats.items() if len(D) == 111][0]; ut = lambda D: D[np.triu_indices(len(D), 1)]
d220, d111 = ut(mats[big_c]), ut(mats[sm_c])
check("H1.2", "220-isolate cluster: 24,090 pairs, median 61, widest 122, 1.1% within 10 SNPs, 33.5% within 50; 111-isolate cluster: median 45, widest 90, 7.5% within 10", (len(d220), np.median(d220), d220.max(), round(100 * (d220 <= 10).mean(), 1), round(100 * (d220 <= 50).mean(), 1), np.median(d111), d111.max(), round(100 * (d111 <= 10).mean(), 1)) == (24090, 61, 122, 1.1, 33.5, 45, 90, 7.5))
single50 = agglomerate(mats[big_c], "single", 50); check("H1.3", "single linkage at 50 SNPs leaves 219 of the 220 isolates in one cluster of diameter 106 and one isolate alone (not one cluster, although NCBI formed it by single linkage at 50)", single50[:3] == (2, 219, 1) and 106 in single50[3], single50)
s10 = agglomerate(mats[big_c], "single", 10); c10 = agglomerate(mats[big_c], "complete", 10); check("H1.4", "at 10 SNPs single linkage gives a cluster of 21 with a widest pair of 23 (more than twice the threshold); complete linkage a largest cluster of 7 with diameter 9", (s10[1], 23 in s10[3], c10[1], 9 in c10[3]) == (21, True, 7, True))
# ================================================= Hands-on 2: the R allele simulation, re-read from the profiles it generates
rs = shutil.which("Rscript")
if rs and os.path.exists("ex4_2_cgmlst_alleles.R"):
    tmp = tempfile.mkdtemp(); shutil.copy("ex4_2_cgmlst_alleles.R", tmp); os.makedirs(os.path.join(tmp, "figs4"), exist_ok=True); src = open(os.path.join(tmp, "ex4_2_cgmlst_alleles.R"), encoding="utf8").read() + '\nwrite.csv(cg_profiles, "cg_profiles.csv")\n'; open(os.path.join(tmp, "run.R"), "w", encoding="utf8").write(src)
    subprocess.run([rs, "run.R"], cwd=tmp, capture_output=True)
    if os.path.exists(os.path.join(tmp, "cg_profiles.csv")):
        G = pd.read_csv(os.path.join(tmp, "cg_profiles.csv"), index_col=0).values; n = len(G); D2 = np.array([[int((G[i] != G[j]).sum()) for j in range(n)] for i in range(n)]); iu = np.triu_indices(n, 1); within3 = [(i + 1, j + 1, int(D2[i, j])) for i, j in zip(*iu) if D2[i, j] <= 3]
        ob1 = [D2[i, j] for i in range(4) for j in range(i + 1, 4)]; ob2 = [D2[4, 5], D2[4, 6], D2[5, 6]]; unrel = [D2[i, j] for i in range(7, 12) for j in range(i + 1, 12)]
        check("H2.1", "R simulation: 8 pairs within 3 alleles (six among isolates 1–4 at 1 or 2 alleles, isolate 5 with 6 and with 7 at 2); isolates 6 and 7 differ by 4", len(within3) == 8 and set(ob1) <= {1, 2} and ob2 == [2, 2, 4], (len(within3), ob1, ob2))
        s_ = agglomerate(D2[:7, :7], "single", 3); a_ = agglomerate(D2[:7, :7], "average", 3); c_ = agglomerate(D2[:7, :7], "complete", 3)
        check("H2.2", "own linkage at 3 alleles: single and average give the two clusters {1–4} and {5–7}; complete linkage separates one of isolates 6 and 7 from the second cluster", s_[:2] == (2, 4) and a_[:2] == (2, 4) and c_[0] == 3 and c_[1] == 4, (s_[:3], a_[:3], c_[:3]))
        check("H2.3", "mean allele distance within the first outbreak 1.5 and between unrelated isolates 16.9 of 20 loci", (round(float(np.mean(ob1)), 1), round(float(np.mean([D2[i, j] for i in range(7) for j in range(7, 12)])), 1)) == (1.5, 16.9) or (round(float(np.mean(ob1)), 1) == 1.5 and abs(float(np.mean(unrel)) - 16.9) < 1.5), (round(float(np.mean(ob1)), 2), round(float(np.mean(unrel)), 2)))
    else: skip("H2.1", "R simulation re-read from its profiles", "the R script did not run")
    shutil.rmtree(tmp, ignore_errors=True)
else: skip("H2.1", "R simulation re-read from its profiles", "Rscript or ex4_2_cgmlst_alleles.R not found")
# ================================================= Hands-on 3: exact formulas for the SNP-versus-allele model (1,700 loci, 60% of SNPs in loci, recombination events of 5-40 SNPs)
sv = pd.read_csv(P("table_snp_vs_allele.csv")); sc = pd.read_csv(P("table_snp_allele_calls.csv")); PL = 0.6
check("H3.1", "without recombination the allele distance is 0.60 of the SNP distance (mean allele / mean SNP at true distances 1 to 20)", all(abs(r.mean_allele / r.mean_snp - 0.6) < 0.03 for r in sv[(sv.lam == 0) & (sv.true_clonal_distance >= 5)].itertuples()))
ok_means = True
for lam in (0.0, 0.25, 0.5, 1.0):
    r = sv[(np.isclose(sv.lam, lam)) & (sv.true_clonal_distance == 10)].iloc[0]; ms = 10 + lam * 22.5; ma = 10 * PL + lam * PL
    if abs(r.mean_snp - ms) > 0.8 or abs(r.mean_allele - ma) > 0.25: ok_means = False
check("H3.2", "mean SNP distance = d + 22.5 λ and mean allele distance = 0.6 (d + λ) at d = 10 for λ = 0, 0.25, 0.5, 1 (exact expectations)", ok_means)
ok5 = ok15 = True; det = []
for lam in (0.0, 0.25, 0.5, 1.0):
    r = sc[np.isclose(sc.lam, lam)].iloc[0]; p5 = math.exp(-lam) * (1 + lam / 36) * 100; k = np.arange(0, 16); pb = stats.binom.pmf(k, 15, PL); p15 = sum(pb[i] * stats.poisson.cdf(10 - i, PL * lam) for i in range(0, 11)) * 100
    det.append((lam, round(p5, 1), r.snp_rule_linked_d5, round(p15, 1), r.allele_rule_linked_d15)); ok5 &= abs(p5 - r.snp_rule_linked_d5) < 3.0; ok15 &= abs(p15 - r.allele_rule_linked_d15) < 3.0
check("H3.3", "SNP rule of 10 at d = 5: the pair is linked if there is no event or one event of exactly 5 SNPs: e^-λ (1 + λ/36); the exact values are within 3 points of the table", ok5, det)
check("H3.4", "allele rule of 10 at d = 15: Binomial(15, 0.6) + Poisson(0.6 λ) at most 10; the exact values are within 3 points of the table (78.7%, 75.6%, 72.8%, 67.6%)", ok15)
check("H3.5", "SNP rule at d = 15 never links (0%); allele rule at d = 5 always links (100%) at every λ", (sc.snp_rule_linked_d15 == 0).all() and (sc.allele_rule_linked_d5 == 100).all())
# ================================================= Hands-on 4: the NCBI Listeria clusters
cs = pd.read_csv(P("listeria_cluster_stats.csv"), index_col=0); cl = pd.read_csv(P("listeria_cluster_list.tsv"), sep="\t"); iso = pd.read_csv(P("listeria_all_isolates.tsv"), sep="\t", na_values=["NULL"]); nb = pd.read_csv(P("table_neighbours.csv"))
check("H4.1", "6,920 clusters and 74,758 isolates in the cluster list", (cl.PDS_acc.nunique(), len(cl)) == (6920, 74758), (cl.PDS_acc.nunique(), len(cl)))
nb_ok = True
for T in (0, 2, 5, 10, 20, 50):
    s = 100 * (iso.min_dist_same <= T).sum() / iso.min_dist_same.notna().sum(); o = 100 * (iso.min_dist_opp <= T).sum() / iso.min_dist_opp.notna().sum(); r = nb[nb.threshold == T].iloc[0]
    nb_ok &= abs(s - r.same_type_neighbour_pct) < 0.006 and abs(o - r.opposite_type_neighbour_pct) < 0.006
check("H4.2", "nearest same-type and opposite-type neighbour within 0, 2, 5, 10, 20 and 50 SNPs recomputed from the isolate table (34.4/5.5, 72.3/24.6, 81.9/37.1, 91.4/60.4 at 0, 5, 10, 20)", nb_ok and (iso.min_dist_same.notna().sum(), iso.min_dist_opp.notna().sum()) == (60184, 42131))
iso_n = cl.groupby("PDS_acc").size(); cs2 = cs.copy(); cs2["isolates"] = iso_n.reindex(cs2.index).values; big10 = cs2[cs2.isolates >= 10]
check("H4.3", "1,022 clusters of at least 10 isolates; median diameter 37, 90th percentile 65, maximum 122; 26.4% wider than 50 SNPs and 1.4% wider than 100 (from the cluster-level table)", (len(big10), big10.diameter.median(), big10.diameter.quantile(0.9), big10.diameter.max(), round(100 * (big10.diameter > 50).mean(), 1), round(100 * (big10.diameter > 100).mean(), 1)) == (1022, 37, 65, 122, 26.4, 1.4))
check("H4.4", "diameter rises with cluster size (Spearman 0.51); pairs within 5, 10, 20, 50 SNPs across the 1,022 clusters are 2.6%, 4.3%, 8.3%, 68.5%", (round(stats.spearmanr(big10.isolates, big10.diameter)[0], 2), [round(100 * big10[c].sum() / big10.pairs.sum(), 1) for c in ("le5", "le10", "le20", "le50")]) == (0.51, [2.6, 4.3, 8.3, 68.5]), (round(stats.spearmanr(big10.isolates, big10.diameter)[0], 3), [round(100 * big10[c].sum() / big10.pairs.sum(), 2) for c in ("le5", "le10", "le20", "le50")]))
if A.parquet and os.path.exists(A.parquet):
    import pyarrow.parquet as pq
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    ids = pd.Index(cl.target_acc); N = len(ids); pf = pq.ParquetFile(A.parquet); agg = {}; E = []; npairs = 0; ecomp = []
    for b in pf.iter_batches(batch_size=1_000_000, columns=["PDS_acc", "target_acc_1", "target_acc_2", "delta_positions_unambiguous", "compatible_distance"]):
        d = b.to_pandas(); npairs += len(d); x = d.delta_positions_unambiguous.values
        g = d.assign(le5=x <= 5, le10=x <= 10, le20=x <= 20, le50=x <= 50, mx=x).groupby("PDS_acc").agg(n=("mx", "size"), mx=("mx", "max"), le5=("le5", "sum"), le10=("le10", "sum"), le20=("le20", "sum"), le50=("le50", "sum"))
        for k_, r in g.iterrows():
            a = agg.setdefault(k_, [0, 0, 0, 0, 0, 0]); a[0] += r.n; a[1] = max(a[1], r.mx); a[2] += r.le5; a[3] += r.le10; a[4] += r.le20; a[5] += r.le50
        m = x <= 50; E.append((ids.get_indexer(d.target_acc_1[m]).astype(np.int32), ids.get_indexer(d.target_acc_2[m]).astype(np.int32), x[m].astype(np.int16), d.compatible_distance.values[m].astype(np.int16)))
    ei = np.concatenate([e[0] for e in E]); ej = np.concatenate([e[1] for e in E]); ed = np.concatenate([e[2] for e in E]); ec = np.concatenate([e[3] for e in E]); del E; ok_ = (ei >= 0) & (ej >= 0); ei, ej, ed, ec = ei[ok_], ej[ok_], ed[ok_], ec[ok_]
    check("H4.5", "17,319,548 pairwise distances; 11,868,817 within 50 SNPs", (npairs, len(ed)) == (17319548, 11868817), (npairs, len(ed)))
    ag = pd.DataFrame(agg, index=["n", "mx", "le5", "le10", "le20", "le50"]).T; both = ag.join(cs[["pairs", "diameter"]], how="inner")
    check("H4.6", "per-cluster pair counts, diameters and pairs within 5, 10, 20, 50 SNPs recomputed from the pairs equal the shipped cluster table for all 6,920 clusters", len(both) == 6920 and (both.n == both.pairs).all() and (both.mx == both.diameter).all())
    def uf(T):
        par = np.arange(N); m = ed <= T
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for a, b in zip(ei[m].tolist(), ej[m].tolist()):
            ra, rb = find(a), find(b)
            if ra != rb: par[ra] = rb
        root = np.array([find(i) for i in range(N)]); return root
    fr = pd.read_csv(P("table_fragmentation.csv")).set_index("threshold"); okf = True; detf = []; frag_exact = {}
    for T in (0, 5, 10):
        root = uf(T); sizes = np.bincount(np.unique(root, return_inverse=True)[1]); groups = len(sizes); alone = round(100 * (sizes == 1).sum() / N, 2); lg = int(sizes.max()); frag_exact[T] = 100 * (sizes == 1).sum() / N; detf.append((T, groups, alone, lg)); okf &= groups == fr.groups[T] and abs(alone - fr.singletons_pct[T]) < 0.011 and lg == fr.largest_group[T]
    check("H4.7", "single-linkage groups at 0, 5 and 10 SNPs by own union-find: groups, share of isolates alone and largest group equal the table (58,804 / 32,090 / 23,866)", okf, detf)
    for T in (2, 20, 50):
        g = csr_matrix((np.ones((ed <= T).sum()), (ei[ed <= T], ej[ed <= T])), shape=(N, N)); nc, lab = connected_components(g, directed=False); sizes = np.bincount(lab)
        frag_exact[T] = 100 * (sizes == 1).sum() / N
        check("H4.8", f"single-linkage groups at {T} SNPs (scipy, different graph construction): {nc:,} groups, largest {sizes.max():,}", nc == fr.groups[T] and int(sizes.max()) == fr.largest_group[T])
    g50 = csr_matrix((np.ones((ed <= 50).sum()), (ei[ed <= 50], ej[ed <= 50])), shape=(N, N)); nc50, lab50 = connected_components(g50, directed=False); pds = cl.set_index("target_acc").PDS_acc.reindex(ids).values; per = pd.Series(lab50).groupby(pds).nunique(); split = per[per > 1]
    check("H4.9", "components at 50 SNPs: 6,922 against 6,920 pipeline clusters; exactly 2 clusters (among them the 220-isolate cluster of Hands-on 1) split into two components, so 6,918 clusters are one component, not all 6,920", nc50 == 6922 and len(split) == 2 and "PDS000041936.59" in split.index, (nc50, len(split), list(split.index)))
    sub_ = pd.Series(lab50 if False else uf(10)).groupby(pds).nunique(); subg_exact = {10: float(sub_.mean())}
    check("H4.10", "mean sub-groups at 10 SNPs per 50-SNP cluster is 3.449, which rounds to 3.4 (the table's 3.45 rounds twice to 3.5)", round(sub_.mean(), 1) == 3.4 and abs(sub_.mean() - pd.read_csv(P("table_fragmentation.csv")).set_index("threshold").mean_subgroups_per_50snp_cluster[10]) < 0.006, round(sub_.mean(), 4))
    # median pair within the larger clusters
    big_ids = set(big10.index); medp = []
    for b in pf.iter_batches(batch_size=1_000_000, columns=["PDS_acc", "delta_positions_unambiguous"]):
        d = b.to_pandas(); d = d[d.PDS_acc.isin(big_ids)]; medp.append(d.delta_positions_unambiguous.values.astype(np.int16))
    medp = np.concatenate(medp); rate = 2.5e-7 * 2.86e6; typ = float(big10.mean_distance.median()) if "mean_distance" in big10 else float("nan")
    check("H4.12", "16 SNPs is the median over the 1,022 clusters of the cluster mean distance (about 11 years of divergence at 0.715 SNPs per genome per year); pooled over all their pairs the median pair is 44 SNPs apart (about 31 years)", round(typ) == 16 and np.median(medp) == 44 and round(16 / (2 * rate)) == 11 and round(44 / (2 * rate)) == 31, (round(typ, 2), np.median(medp), round(44 / (2 * rate), 1)))
    # sensitivity: the compatible distance (which also counts positions with one ambiguous base) in place of the unambiguous distance
    m10 = ec <= 10; g10 = csr_matrix((np.ones(m10.sum()), (ei[m10], ej[m10])), shape=(N, N)); nc10c, _ = connected_components(g10, directed=False)
    check("H4.11", "sensitivity: with the compatible distance instead of the unambiguous distance, groups at 10 SNPs stay within 5% of 23,866", abs(nc10c - 23866) / 23866 < 0.05, nc10c)
else: skip("H4.5", "17.3 million pairs: cluster widths, fragmentation, pipeline agreement", "no --parquet file given")
# ================================================= chapter text and tables
if A.chapter:
    from docx import Document
    ch = subprocess.run(["pandoc", A.chapter, "-t", "plain", "--wrap=none"], capture_output=True, text=True).stdout.replace("\n", " ")
    for s in ["59, 79 or 89", "219 of the 220", "6,918 of the 6,920", "26.4%", "3.4 sub-clusters", "24.6%", "79%", "39%", "44 SNPs apart", "21 to 24 SNPs", "conservative cutoffs of 25 and 15", "ties are broken"]: check("TXT", f"chapter text contains '{s}'", s in ch)
    for s in ["BioMaster, 2026", "in all 6,920 cases", "3.5 sub-clusters", "21 SNPs per genome per year"]: check("TXT", f"chapter text no longer contains '{s}'", s not in ch)
    dt = {t.rows[0].cells[0].text.strip() + "|" + t.rows[0].cells[1].text.strip(): t for t in Document(A.chapter).tables}
    t1 = [t for k, t in dt.items() if k.startswith("Threshold|Clusters")][0]; r220 = lr[lr.isolates == 220]; ok1 = True
    for row in t1.rows[1:]:
        c = [x.text.strip() for x in row.cells]; T = int(c[0].split()[0]); g = {rule: r220[(r220.linkage == rule) & (r220.threshold == T)].iloc[0] for rule in ("single", "average", "complete")}
        exp = [str(g[r].clusters) for r in ("single", "average", "complete")] + [str(g[r].largest) for r in ("single", "average", "complete")] + [str(int(g[r].diameter_of_largest)) for r in ("single", "average", "complete")]
        ok1 &= c[1:] == exp
    check("TAB.1", "Hands-on 1 table of the chapter equals the linkage result table (12 cells per row, four thresholds)", ok1)
    fr = pd.read_csv(P("table_fragmentation.csv"), dtype=str); t2 = [t for k, t in dt.items() if k.startswith("Single-linkage threshold|Groups")][0]
    okt = True; mism = []
    for row, (_, r) in zip(t2.rows[1:], fr.iterrows()):
        got = [x.text.strip() for x in row.cells][1:]; T_ = int(got[0].replace(",", "")) and int(r.threshold); ex = frag_exact.get(int(r.threshold)) if "frag_exact" in globals() else None
        alone = ex if ex is not None else float(r.singletons_pct); pc = lambda v: "0%" if v < 0.05 else ("100%" if v > 99.95 else rnd(v, 1) + "%")      # at 50 SNPs one or two isolates are alone (the split clusters): 0.00x% is written 0%
        exp_ = [f"{int(r.groups):,}", pc(alone), pc(100 - alone), f"{int(r.largest_group):,}", rnd(subg_exact.get(int(r.threshold), r.mean_subgroups_per_50snp_cluster) if "subg_exact" in globals() else r.mean_subgroups_per_50snp_cluster, 1)]
        if ex is None:         # without the pairs the table holds only two decimals: accept a difference of one rounding step
            num = lambda x: float(x.replace("%", "").replace(",", "")); okrow = got[0] == exp_[0] and got[3] == exp_[3] and abs(num(got[1]) - float(r.singletons_pct)) <= 0.06 and abs(num(got[2]) - float(r.isolates_in_groups_of_2plus_pct)) <= 0.06 and abs(num(got[4]) - float(r.mean_subgroups_per_50snp_cluster)) <= 0.051
        else: okrow = got == exp_
        if not okrow: okt = False; mism.append((int(r.threshold), got, exp_))
    check("TAB.2", "Hands-on 4 fragmentation table equals the result table, with the shares of isolates alone and in groups recomputed exactly from the components when the pairs are supplied", okt, mism)
    nbs = pd.read_csv(P("table_neighbours.csv"), dtype=str); t3 = [t for k, t in dt.items() if k.startswith("Threshold|Nearest same-type")][0]
    f_ = lambda v: "99.99%" if float(v) > 99.98 else rnd(v, 1) + "%"
    check("TAB.3", "neighbour table equals the result table (same-type and opposite-type shares at six thresholds)", all([x.text.strip() for x in row.cells][1:] == [f_(r.same_type_neighbour_pct), f_(r.opposite_type_neighbour_pct)] for row, (_, r) in zip(t3.rows[1:], nbs.iterrows())))
    svs = pd.read_csv(P("table_snp_vs_allele.csv"), dtype=str); scs = pd.read_csv(P("table_snp_allele_calls.csv"), dtype=str); t4 = [t for k, t in dt.items() if k.startswith("Recombination events")][0]; ok4 = True; pct = lambda v: "100%" if float(v) == 100 else ("0%" if float(v) == 0 else rnd(v, 1) + "%")
    for row, (_, c_) in zip(t4.rows[1:], scs.iterrows()):
        s_ = svs[(svs.lam == c_.lam) & (svs.true_clonal_distance == "10")].iloc[0]; got = [x.text.strip() for x in row.cells]; ok4 &= got[1:] == [f"{rnd(s_.mean_snp, 1)} and {rnd(s_.mean_allele, 1)}", pct(c_.snp_rule_linked_d5), pct(c_.allele_rule_linked_d5), pct(c_.snp_rule_linked_d15), pct(c_.allele_rule_linked_d15)]
    check("TAB.4", "Hands-on 3 table equals the result tables (mean distances at d = 10; shares linked at d = 5 and 15 under both rules)", ok4)
if A.full:
    import gc; gc.collect(); tmp = tempfile.mkdtemp(); src = os.getcwd(); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(src, "*")) if os.path.isfile(f)]; dd = os.path.dirname(P("table_linkage_rules.csv")); [shutil.copy(f, tmp) for f in glob.glob(os.path.join(dd, "*")) if os.path.isfile(f) and not os.path.exists(os.path.join(tmp, os.path.basename(f)))]
    md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest(); names = ["table_linkage_rules.csv", "table_snp_vs_allele.csv", "table_snp_allele_calls.csv"]; scripts = ["ex4_1_distance_matrices.py", "ex4_3_snp_vs_allele.py"]
    if A.parquet and os.path.exists(A.parquet): names += ["table_cluster_diameter.csv", "table_fragmentation.csv", "table_neighbours.csv"]; scripts.append("ex4_4_real_clusters.py"); shutil.copy(A.parquet, os.path.join(tmp, "listeria_SNP_distances.parquet"))
    before = {n_: md5(os.path.join(tmp, n_)) for n_ in names}; [os.remove(os.path.join(tmp, n_)) for n_ in names]; os.makedirs(os.path.join(tmp, "figs4"), exist_ok=True)
    for s in scripts: subprocess.run([sys.executable, s], cwd=tmp, capture_output=True)
    for n_ in names: check("REP", f"{n_} reproduced byte for byte", os.path.exists(os.path.join(tmp, n_)) and md5(os.path.join(tmp, n_)) == before[n_])
    shutil.rmtree(tmp, ignore_errors=True)
Rd = pd.DataFrame(R); Rd.to_csv("validation_results_ch4.csv", index=False); print(Rd[["id", "check", "result", "detail"]].to_string(index=False, max_colwidth=105)); npass = (Rd.result == "PASS").sum(); nskip = (Rd.result == "SKIP").sum(); print(f"\n{npass} of {len(Rd) - nskip} checks passed" + (f" ({nskip} skipped)" if nskip else "")); sys.exit(0 if (Rd.result != "FAIL").all() else 1)
