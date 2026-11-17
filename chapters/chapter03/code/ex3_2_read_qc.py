"""
Chapter 3, Hands-on 2 - Read-level quality control on public Illumina, Oxford Nanopore and PacBio runs of Salmonella enterica from the ENA.
For each run the first part of the compressed FASTQ file is downloaded (a byte range), parsed with a small FASTQ reader, and summarised: read length, Phred-based accuracy, per-base quality, GC content,
adapter content, duplication and the effect of the sliding-window trimming rule of the draft chapter (window 5, mean Q below 20). Output: table_read_qc.csv, per_base_quality.csv (for the R plot), Figure 2.
"""
import gzip, io, os, subprocess, zlib, urllib.request
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
RUNS = {"MiSeq (2014)": ("DRR090587", "ena_illumina.csv", 25_000_000, 100_000), "NovaSeq 6000 (2024)": ("SRR30474864", "ena2_illumina.csv", 20_000_000, 100_000), "Nanopore MinION": ("DRR1060924", "ena2_ont.csv", 40_000_000, 100_000), "PacBio Sequel IIe": ("SRR33434758", "ena2_pacbio.csv", 40_000_000, 100_000)}
ADAPTER = "AGATCGGAAGAGC"

def fetch(acc, table, nbytes):
    path = f"{acc}_prefix.fastq"
    if os.path.exists(path): return path
    t = pd.read_csv(table).set_index("run_accession"); url = "https://" + t.loc[acc, "fastq_ftp"].split(";")[0]
    gz = f"{acc}_prefix.gz"
    for attempt in range(6):                                                                                                   # the ENA mirror answers 503 under load, so retry with a pause
        r = subprocess.run(["curl", "-s", "-L", "-m", "240", "-r", f"0-{nbytes - 1}", "-o", gz, "-w", "%{http_code}", url], capture_output=True, text=True)
        if r.stdout.strip() in ("200", "206") and os.path.exists(gz) and os.path.getsize(gz) > 1000 and open(gz, "rb").read(2) == b"\x1f\x8b": break          # accept only a real gzip stream, not an error page
        import time; time.sleep(3 * (attempt + 1))
    if not (os.path.exists(gz) and open(gz, "rb").read(2) == b"\x1f\x8b"):                                                      # the mirror stayed unavailable
        demo = f"{acc}_demo.fastq.gz"
        if not os.path.exists(demo): raise RuntimeError(f"could not download {acc} from the ENA and no {demo} demo file is present")
        print(f"NOTE: {acc} could not be downloaded; using the shipped demo subset of the first reads, so the numbers will differ slightly from the chapter"); out = gzip.open(demo, "rb").read(); open(path, "wb").write(out); return path
    out = zlib.decompressobj(31).decompress(open(gz, "rb").read())                                                             # a truncated gzip stream still yields everything up to the cut
    open(path, "wb").write(out); return path

def read_fastq(path, limit):
    with open(path, "rb") as f: lines = f.read().split(b"\n")
    n = min((len(lines) - 1) // 4, limit); recs = []
    for i in range(n - 1):                                                                                                     # the last record may be cut short, so it is dropped
        seq, qual = lines[4 * i + 1], lines[4 * i + 3]
        if len(seq) == len(qual) and lines[4 * i + 2].startswith(b"+"): recs.append((seq.decode(), np.frombuffer(qual, dtype=np.uint8).astype(np.int16) - 33))
    return recs

def sliding_window_trim(scores, window_size=5, threshold=20):
    """The rule of the draft chapter: cut the read at the start of the first window whose mean quality falls below the threshold."""
    if len(scores) < window_size: return len(scores)
    c = np.convolve(scores, np.ones(window_size) / window_size, mode="valid"); bad = np.where(c < threshold)[0]; return int(bad[0]) if len(bad) else len(scores)

def n50(lengths):
    s = np.sort(lengths)[::-1]; return int(s[np.searchsorted(np.cumsum(s), s.sum() / 2)])

rows, perbase, store = [], [], {}
for name, (acc, table, nbytes, limit) in RUNS.items():
    path = fetch(acc, table, nbytes); recs = read_fastq(path, limit); L = np.array([len(s) for s, _ in recs]); Q = [q for _, q in recs]
    allq = np.concatenate(Q); perr = 10 ** (-allq / 10.0); mean_q_prob = -10 * np.log10(perr.mean())                           # accuracy averaged on the probability scale, not on the Phred scale
    read_mean_q = np.array([-10 * np.log10((10 ** (-q / 10.0)).mean()) for q in Q]); seqs = [s for s, _ in recs]
    gc = np.array([(s.count("G") + s.count("C")) / max(len(s), 1) for s in seqs]); adapter = np.mean([ADAPTER in s for s in seqs]) * 100
    dup = np.nan if L.mean() > 1000 else 100 * (1 - len({s[:40] for s in seqs}) / len(seqs)); trim_len = np.array([sliding_window_trim(q) for q in Q]); distinct_q = len(np.unique(allq))
    rows.append(dict(run=name, accession=acc, reads=len(recs), bases_mb=L.sum() / 1e6, mean_length=L.mean(), median_length=np.median(L), N50=n50(L), mean_Q_arith=allq.mean(), mean_Q_prob=mean_q_prob, error_per_kb=1000 * perr.mean(), q30_pct=100 * (allq >= 30).mean(), q20_pct=100 * (allq >= 20).mean(), distinct_q_values=distinct_q, gc_pct=100 * gc.mean(), gc_sd_pct=100 * gc.std(), adapter_pct=adapter, dup40_pct=dup, retained_after_trim_pct=100 * trim_len.sum() / L.sum(), reads_lost_pct=100 * (trim_len < 50).mean()))
    store[name] = dict(L=L, read_mean_q=read_mean_q, gc=gc)
    maxpos = 150 if L.mean() < 1000 else 300; M = np.full((len(Q), maxpos), np.nan)
    for i, q in enumerate(Q): k = min(len(q), maxpos); M[i, :k] = q[:k]
    for p in range(maxpos):
        col = M[:, p]; col = col[~np.isnan(col)]
        if len(col) >= 200: perbase.append(dict(run=name, position=p + 1, mean=col.mean(), q10=np.percentile(col, 10), q25=np.percentile(col, 25), q75=np.percentile(col, 75), q90=np.percentile(col, 90), n=len(col)))
T = pd.DataFrame(rows); T.round(3).to_csv("table_read_qc.csv", index=False); pd.DataFrame(perbase).round(3).to_csv("per_base_quality.csv", index=False); print(T.round(2).T.to_string())
PB = pd.DataFrame(perbase)
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [1.1, 1, 1]}); cols = {"MiSeq (2014)": "#2b6cb0", "NovaSeq 6000 (2024)": "#2f855a", "Nanopore MinION": "#dd6b20", "PacBio Sequel IIe": "#805ad5"}
ax[0].axhspan(28, 41, color="darkgreen", alpha=0.10); ax[0].axhspan(20, 28, color="orange", alpha=0.12); ax[0].axhspan(0, 20, color="darkred", alpha=0.08)
for r in ["MiSeq (2014)", "NovaSeq 6000 (2024)"]: x = PB[PB.run == r]; ax[0].plot(x.position, x["mean"], color=cols[r], lw=2, label=r); ax[0].fill_between(x.position, x.q25, x.q75, color=cols[r], alpha=0.18); ax[0].plot(x.position, x.q10, color=cols[r], ls="--", lw=1)
ax[0].set(xlabel="Position in read (bp)", ylabel="Phred score (Q)", title="Per-base quality of two Illumina runs", ylim=(0, 41)); ax[0].legend(frameon=False, fontsize=8.5, loc="lower left")
for r, c in cols.items():
    Ls = np.sort(store[r]["L"])[::-1]; ax[1].plot(Ls, np.cumsum(Ls) / Ls.sum() * 100, color=c, lw=2, label=f"{r} (N50 {n50(Ls):,})")
ax[1].set_xscale("log"); ax[1].invert_xaxis(); ax[1].set(xlabel="Read length (bp, log scale)", ylabel="Bases in reads at least this long (%)", title="Where the bases are"); ax[1].legend(frameon=False, fontsize=7.8, loc="lower right")
for r, c in cols.items(): ax[2].hist(store[r]["read_mean_q"], bins=np.arange(0, 45, 1), histtype="step", color=c, lw=1.8, density=True, label=r)
ax[2].set(xlabel="Mean quality of a read (error-probability scale)", ylabel="Density", title="Accuracy of whole reads"); fig.tight_layout(); fig.savefig("figs3/fig2_read_qc.png", dpi=200); plt.close(fig); print("figure written")
