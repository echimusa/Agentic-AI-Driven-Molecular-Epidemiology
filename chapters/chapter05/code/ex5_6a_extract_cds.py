"""
Chapter 5, Hands-on 5 preparation: cuts the coding region (1,701 nt, start motif ATGAAGAC) out of the H3N2 HA alignment and writes h3n2_cds.fasta and one file per codon position (h3n2_cp1-3.fasta).
Checks that every sequence has the start motif, that the block has no gaps, and that no in-frame stop codon occurs.
"""
import re
import numpy as np
from phylo_utils import read_aln
names, X = read_aln("h3n2_ha_aln.fasta"); n, L = X.shape; starts = []
for i in range(n):
    col = [j for j in range(L) if X[i, j] != "-"]; ung = "".join(X[i, col]); m = re.search("ATGAAGAC", ung); starts.append(col[m.start()] if m else None)
print("sequences with the start motif ATGAAGAC:", sum(s is not None for s in starts), "of", n); st = int(np.median([s for s in starts if s is not None])); cds = X[:, st:st + 1701]
print("CDS block:", cds.shape, "| gap fraction", float((cds == "-").mean()))
stops = sum(("".join(cds[i]).replace("-", "N")[k:k + 3] in ("TAA", "TAG", "TGA")) for i in range(n) for k in range(0, 1698, 3)); print("in-frame internal stop codons:", stops)
with open("h3n2_cds.fasta", "w") as f: [f.write(f">{nm}\n{''.join(cds[i])}\n") for i, nm in enumerate(names)]
for p in range(3):
    with open(f"h3n2_cp{p + 1}.fasta", "w") as f: [f.write(f">{nm}\n{''.join(cds[i, p::3])}\n") for i, nm in enumerate(names)]
