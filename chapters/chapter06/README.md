# Chapter 6: Phylodynamics and Molecular Clocks

Code, data and checks for the exercises of Chapter 6 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 48 files, 34.2 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch6.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | temporal signal in real SARS-CoV-2 genomes | Python | `ex6_root_to_tip.py` |
| 2 | skyline plots, checked against a known truth | R | `ex6_skyline.R` |
| 3 | time-scaled trees with TreeTime | Python | `ex6_treetime.py` |
| 4 | time-dependent rates in a model | Python | `ex6_time_dependent.py` |
| 5 | Bayesian dating, BETS and Bayes factors | BEAST 2 | `ex6_beast.py`, `run_all.sh`, `ex6_beast_analysis.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 6                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 6 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 6 --script ex6_root_to_tip.py   # run one exercise; the output stays in chapters/chapter06/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 6
./scripts/run-chapter.ps1 -Chapter 6 -Full
./scripts/run-chapter.ps1 -Chapter 6 -Script ex6_root_to_tip.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 6`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 43 of 43, in 20 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, mappy, phylo-treetime
- R packages: ggplot2
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: FastTree
- Heavy analyses (optional): beast2
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

BEAST 2 (Hands-on 5) and the long TreeTime reruns (`-Full -Sim`, about 25 minutes) are not part of the quick check; the BEAST 2 results are shipped.

## Data

48 files, 34.2 MB. Largest: `sars2.fasta` (13.2 MB), `ebola_aln.fasta` (3.8 MB), `ps_het_relaxed.xml` (1.2 MB).
The result tables are `table_growth_R.csv, table_treetime_drt.csv, table_time_dependent.csv, table_time_dependent_dating.csv, table_skyline_ebola.csv, table_treetime.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
