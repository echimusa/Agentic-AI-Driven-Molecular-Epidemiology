# Chapter 11: Mendelian Randomisation and Causal Inference

Code, data and checks for the exercises of Chapter 11 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 29 files, 14.6 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch11.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | body mass index and hospitalised COVID-19 | Python | `ex11_1_bmi_covid.py` |
| 2 | robust estimators and diagnostics | Python and R | `ex11_1_bmi_covid.py`, `ex11_2b_r_check.R` |
| 3 | a known truth | Python | `ex11_3_simulation.py` |
| 4 | drug-target MR at IL6R | Python | `ex11_4_cis_mr.py` |
| 5 | why alleles must be harmonised | Python | `ex11_5_harmonise.py` |
| 6 | vitamin D and hospitalised COVID-19 | Python | `ex11_6_vitamin_d.py` |
| 7 | reconciling an observational association with a null MR | Python | `ex11_7_reconcile.py` |
| 8 | one sample, two samples and the winner’s curse | Python | `ex11_8_one_two_sample.py` |
| 9 | a binary outcome, testing against estimating | Python | `ex11_9_binary_outcome.py` |
| 10 | invalid instruments, the weighted mode and multivariable MR | Python | `ex11_10_mode_mvmr.py` |
| 11 | a within-family design against dynastic effects | Python | `ex11_11_family_design.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 11                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 11 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 11 --script ex11_1_bmi_covid.py   # run one exercise; the output stays in chapters/chapter11/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 11
./scripts/run-chapter.ps1 -Chapter 11 -Full
./scripts/run-chapter.ps1 -Chapter 11 -Script ex11_1_bmi_covid.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 11`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 40 passed, 0 failed, in 1.8 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml
- R packages: none
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: bcftools, tabix
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

29 files, 14.6 MB. Largest: `chr1_il6r_eur.gt` (14.1 MB), `bmi_assoc.csv` (0.2 MB), `table_mr_bmi_instruments.csv` (0.1 MB).
The result tables are `table_simulation.csv, table_mode_bmi.csv, table_vitd_methods.csv, table_vitd_summary.csv, table_cis_mr.csv, table_reconcile.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
