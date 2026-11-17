# Chapter 7: Transmission Inference and Outbreak Investigation

Code, data and checks for the exercises of Chapter 7 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 25 files, 5.1 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch7.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | who infected whom in a measles outbreak | Python | `ex1_hagelloch.py` |
| 2 | what do genomes add? Ebola virus in Sierra Leone | Python | `ex2_ebola.py` |
| 3 | describe the outbreak in time, place and person, and size a tracing window | Python | `ex4_descriptive.py` |
| 4 | Bayesian reconstruction with outbreaker2 | R | `ex3_outbreaker2.R` |
| 5 | TransPhylo on simulated outbreaks with a known truth | R | `ex5_transphylo.R`, `ex5_transphylo_analysis.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 7                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 7 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 7 --script ex1_hagelloch.py   # run one exercise; the output stays in chapters/chapter07/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 7
./scripts/run-chapter.ps1 -Chapter 7 -Full
./scripts/run-chapter.ps1 -Chapter 7 -Script ex1_hagelloch.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 7`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 43 of 43, in 41 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, mappy, networkx, pyreadr
- R packages: outbreaker2, TransPhylo, outbreaks, coda
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

The 12 TransPhylo runs of Hands-on 5 take about 40 minutes on one processor; their results are shipped and the check reads them.

## Data

25 files, 5.1 MB. Largest: `ebola.fasta` (4.8 MB), `hagelloch_outbreaker2.rds` (0.2 MB), `transphylo_pairs.csv` (0.0 MB).
The result tables are `table_transphylo_short_calibration.csv, table_generations.csv, table_transphylo_calibration.csv, table_descriptive_summary.csv, table_transphylo_short_runs.csv, table_tracing.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
