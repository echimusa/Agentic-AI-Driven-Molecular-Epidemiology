# Chapter 4: Genomic Clustering: SNP and cgMLST Approaches

Code, data and checks for the exercises of Chapter 4 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 10 files, 8.0 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch4.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | distance matrices, thresholds and linkage rules | Python | `ex4_1_distance_matrices.py` |
| 2 | allele distances, pairs and clusters | R | `ex4_2_cgmlst_alleles.R` |
| 3 | SNP distance against allele distance | Python | `ex4_3_snp_vs_allele.py` |
| 4 | chaining, fragmentation and links in real clusters | Python | `ex4_4_real_clusters.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 4                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 4 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 4 --script ex4_1_distance_matrices.py   # run one exercise; the output stays in chapters/chapter04/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 4
./scripts/run-chapter.ps1 -Chapter 4 -Full
./scripts/run-chapter.ps1 -Chapter 4 -Script ex4_1_distance_matrices.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 4`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 40 of 40, in 10 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, pyarrow, seaborn
- R packages: ggplot2, gplots
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

10 files, 8.0 MB. Largest: `listeria_cluster_list.tsv` (4.4 MB), `listeria_all_isolates.tsv` (3.2 MB), `listeria_cluster_stats.csv` (0.3 MB).
The result tables are `table_snp_vs_allele.csv, table_fragmentation.csv, table_neighbours.csv, table_snp_allele_calls.csv, table_cluster_diameter.csv, table_linkage_rules.csv`.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
