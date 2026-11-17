# Chapter 2: Epidemiological Study Design, Bias and Sampling for Molecular Studies

Code, data and checks for the exercises of Chapter 2 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 12 files, 2.3 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch2.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | nested case–control sampling | Python and R | `ex2_1_nested_case_control.py`, `ex2_1b_infert_clogit.R` |
| 2 | batch effects | R and Python | `ex2_2a_batch_pca.R`, `ex2_2_batch_effects.py` |
| 3 | population stratification | Python | `ex2_3_stratification.py` |
| 4 | sampling for molecular studies | Python | `ex2_4_sampling.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 2                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 2 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 2 --script ex2_1_nested_case_control.py   # run one exercise; the output stays in chapters/chapter02/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 2
./scripts/run-chapter.ps1 -Chapter 2 -Full
./scripts/run-chapter.ps1 -Chapter 2 -Script ex2_1_nested_case_control.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 2`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 61 of 61, in 4.2 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, scikit-learn, statsmodels
- R packages: ggplot2
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

12 files, 2.3 MB. Largest: `asm_batch_enteritidis.csv.gz` (1.2 MB), `owid_cases_2021_2022.csv` (0.7 MB), `lapis_daily_2021_2022.csv` (0.4 MB).
The result tables are `table_batch_design.csv, table_stratification_strength.csv, table_batch_labs.csv, table_sequencing_coverage.csv, table_sample_size.csv, table_infert_clogit.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
