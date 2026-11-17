# Chapter 12: Genomic Surveillance Systems, Ethics and Data Governance

Code, data and checks for the exercises of Chapter 12 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 20 files, 23.5 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch12.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | completeness and timeliness from public metadata | Python | `ex12_1_indicators.py` |
| 2 | representativeness and equity | Python and R | `ex12_2_equity.py`, `ex12_2b_gini.R` |
| 3 | privacy risk and utility | Python | `ex12_3_privacy.py` |
| 4 | how reporting delay postpones cluster detection | Python | `ex12_4_detection.py` |
| 5 | how many sequences a week does an aim need? | Python | `ex12_5_sizing.py` |
| 6 | where the sequences are taken decides when a variant is seen | Python | `ex12_6_allocation.py` |
| 7 | what free text costs | Python | `ex12_7_standards.py` |
| 8 | what an embargo costs | Python | `ex12_8_embargo.py` |
| 9 | check a surveillance protocol | Python | `ex12_9_protocol.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 12                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 12 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 12 --script ex12_1_indicators.py   # run one exercise; the output stays in chapters/chapter12/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 12
./scripts/run-chapter.ps1 -Chapter 12 -Full
./scripts/run-chapter.ps1 -Chapter 12 -Script ex12_1_indicators.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 12`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 50 passed, 0 failed, in 48 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml
- R packages: none
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

20 files, 23.5 MB. Largest: `cluster_list_us.tsv` (16.5 MB), `sal_surv.csv.gz` (6.8 MB), `wb_countries.json` (0.1 MB).
The result tables are `table_sizing.csv, table_protocol_check.csv, table_sizing_real.csv, table_income_groups.csv, table_embargo.csv, table_allocation.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
