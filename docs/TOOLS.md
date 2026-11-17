# Tools and packages by chapter

Reference versions: Ubuntu 24.04 (Python 3.12.3, R 4.3.3, IQ-TREE 2.0.7, FastTree 2.1.11, MAFFT 7.505, bcftools 1.19, PLINK 1.90b7.2, pandoc 3.1.3). Python packages are pinned in `requirements.txt`; R packages are installed from CRAN by `r-packages.R`.

| Ch | Python packages (beyond numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml) | R packages | Needed for the quick check | Needed to rebuild inputs | Heavy (optional) |
|---|---|---|---|---|---|
| 1 | biopython | none | none | none | none |
| 2 | scikit-learn, statsmodels | ggplot2 | none | none | none |
| 3 | biopython, scikit-learn, statsmodels | ggplot2 | none | none | none |
| 4 | pyarrow, seaborn | ggplot2, gplots | none | none | none |
| 5 | biopython | ape, ggplot2 | iqtree2, FastTree | none | none |
| 6 | biopython, mappy, phylo-treetime | ggplot2 | none | FastTree | beast2 |
| 7 | biopython, mappy, networkx, pyreadr | outbreaker2, TransPhylo, outbreaks, coda | none | none | none |
| 8 | biopython, phylo-treetime, shapely | none | none | mafft, FastTree | none |
| 9 | biopython, networkx, statsmodels | ggplot2 | none | none | amrfinder, rgi |
| 10 | none | ggplot2 | none | bcftools, tabix, plink1.9 | none |
| 11 | none | none | none | bcftools, tabix | none |
| 12 | none | none | none | none | none |
| 13 | biopython, phylo-treetime | none | none | mafft, FastTree | none |

Three levels of need:

1. **Quick check and full rerun from the shipped data**: Python, the R packages above, pandoc, and for Chapter 5 only, IQ-TREE 2 and FastTree. In the test every other command-line tool was removed and the network blocked and all chapters but Chapter 5 still reproduced their tables; IQ-TREE 2 and FastTree were removed together, so it is not established that both are needed.
2. **Rebuilding inputs from public sources**: MAFFT and FastTree (alignments and trees of Chapters 6, 8 and 13), bcftools and tabix (1000 Genomes and summary-statistics extracts of Chapters 10 and 11), PLINK 1.9 (Chapter 10), plus a network connection. `scripts/setup.sh --rebuild` (or `scripts/setup.ps1 -Rebuild`) installs them on Ubuntu.
3. **Heavy analyses**: BEAST 2 (Chapter 6, Hands-on 5), AMRFinderPlus and RGI (Chapter 9, Hands-on 5 and 6). See [HEAVY-ANALYSES.md](HEAVY-ANALYSES.md).

Ubuntu package names: IQ-TREE `iqtree` (provides `iqtree2`), `fasttree`, `mafft`, `bcftools`, `tabix`, `plink1.9`, `pandoc`.
