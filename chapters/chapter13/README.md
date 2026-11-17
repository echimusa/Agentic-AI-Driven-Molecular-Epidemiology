# Chapter 13: Integrating Molecular Data with Transmission Models

Code, data and checks for the exercises of Chapter 13 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 27 files, 9.5 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch13.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | reproduction number and growth rate from public case data | Python | `ex13_1_ebola_rt.py`, `ex13_1b_r_check.R` |
| 2 | what the genomes show, and what they cannot | Python | `ex13_2_genome_growth.py` |
| 3 | a known truth, and what genomes add to case data | Python | `ex13_3_simulation.py` |
| 4 | variant replacement and a two-variant renewal forecast | Python | `ex13_4_variants.py` |
| 5 | the SIR model, and why R = exp(r T) is not enough | Python | `ex13_5_sir_gi.py` |
| 6 | different windows, different answers | Python | `ex13_6_windows.py` |
| 7 | a partly sampled outbreak with a host attribute | Python | `ex13_7_partial_sampling.py` |
| 8 | two databases, one outbreak | Python | `ex13_8_integration_audit.py` |
| 9 | check a capstone analysis plan and briefing | Python | `ex13_9_capstone_check.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 13                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 13 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 13 --script ex13_1_ebola_rt.py   # run one exercise; the output stays in chapters/chapter13/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 13
./scripts/run-chapter.ps1 -Chapter 13 -Full
./scripts/run-chapter.ps1 -Chapter 13 -Script ex13_1_ebola_rt.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 13`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 42 passed, 0 failed, in 24 s.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, phylo-treetime
- R packages: none
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: mafft, FastTree
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

27 files, 9.5 MB. Largest: `ebola.fasta` (4.8 MB), `ebola_aln.fasta` (3.8 MB), `ebola_sl_linelist.csv` (0.8 MB).
The result tables are `table_sir.csv, table_ebola_rt.csv, table_gi_shape.csv, table_simulation13.csv, table_variant_summary.csv, table_genome_growth.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
