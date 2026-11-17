# Chapter 10: Host Genetic Epidemiology and Genome-Wide Association

Code, data and checks for the exercises of Chapter 10 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 35 files, 2.5 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch10.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | what GWAS of infection have found, and in whom | Python | `ex10_1_catalog.py` |
| 2 | quality control and loci from public summary statistics | Python and R | `ex10_2_qc_manhattan.py`, `ex10_2b_r_check.R` |
| 3 | population stratification with a known truth | Python | `ex10_3_stratification.py` |
| 4 | fine-mapping the 3p21.31 locus | Python | `ex10_4_finemap.py` |
| 5 | the quality-control command, ancestry and Hardy–Weinberg | Python and PLINK | `ex10_5_qc_structure.py` |
| 6 | how many people? | Python | `ex10_6_power_lambda.py` |
| 7 | why a signal may not replicate, the 3p21.31 variants in five ancestries | Python | `ex10_7_ancestry_3p21.py` |
| 8 | admixture mapping with a known truth | Python | `ex10_8_admixture_mapping.py` |
| 9 | recomputing the worked examples, and a checklist for reading a GWAS paper | Python | `ex10_9_notes_audit.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 10                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 10 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 10 --script ex10_1_catalog.py   # run one exercise; the output stays in chapters/chapter10/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 10
./scripts/run-chapter.ps1 -Chapter 10 -Full
./scripts/run-chapter.ps1 -Chapter 10 -Script ex10_1_catalog.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 10`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 39 passed, 0 failed, in 30 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml
- R packages: ggplot2
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: bcftools, tabix, plink1.9
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

35 files, 2.5 MB. Largest: `chr22_group_geno_counts.tsv.gz` (1.2 MB), `chr3_region_eur.gt.gz` (0.7 MB), `chr3_region_all.gt.gz` (0.2 MB).
The result tables are `table_example_allelic.csv, table_loci.csv, table_notes_audit_summary.csv, table_catalog_ancestry.csv, table_lambda_scaling.csv, table_power_min_or.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
