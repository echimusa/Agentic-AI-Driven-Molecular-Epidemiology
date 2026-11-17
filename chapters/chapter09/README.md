# Chapter 9: Antimicrobial Resistance Genomics in a One Health Context

Code, data and checks for the exercises of Chapter 9 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 141 files, 37.6 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch9.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | from public isolates to a One Health sampling frame | Python | `ex9_1_sampling_frame.py`, `ex9_7_metadata_completeness.py` |
| 2 | resistance determinants and multidrug resistance by sector | Python and R | `ex9_2_resistance_by_sector.py`, `ex9_2b_glm.R` |
| 3 | genotype versus phenotype | Python | `ex9_3_genotype_phenotype.py`, `ex9_3b_discordance.py` |
| 4 | cross-sector clusters and genes that travel together | Python | `ex9_4_clusters_cooccurrence.py` |
| 5 | calling resistance genes from genomes with AMRFinderPlus | Python | `ex9_5_amrfinder_genomes.py` |
| 6 | two detection tools on the same genomes | Python and R | `ex9_5b_two_tools.py`, `ex9_5c_summarise.R` |
| 7 | how many cross-sector pairs survive each SNP threshold? | Python | `ex9_6_snp_thresholds.py` |
| 8 | reading three cross-sector studies and recomputing their numbers | Python | `ex9_8_study_audit.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 9                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 9 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 9 --script ex9_1_sampling_frame.py   # run one exercise; the output stays in chapters/chapter09/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 9
./scripts/run-chapter.ps1 -Chapter 9 -Full
./scripts/run-chapter.ps1 -Chapter 9 -Script ex9_1_sampling_frame.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 9`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 52 passed, 0 failed, in 95 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, networkx, statsmodels
- R packages: ggplot2
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): amrfinder, rgi
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

AMRFinderPlus and RGI (Hands-on 5 and 6) need the databases described in docs/HEAVY-ANALYSES.md; their raw outputs are shipped and the check reads them.

## Data

141 files, 37.6 MB. Largest: `salmonella_cluster_list.tsv` (20.8 MB), `sal_us_frame.csv.gz` (5.5 MB), `ReferenceGeneCatalog.txt` (2.4 MB).
The result tables are `table_metadata_completeness.csv, table_serovar_mix.csv, table_concordance.csv, table_snp_min_distance.csv, table_class_counts_r.csv, table_mixed_clusters.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
