# Chapter 1: Principles, Designs, and Applications of Molecular Epidemiology

Code, data and checks for the exercises of Chapter 1 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 20 files, 4.9 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch1.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | discrimination and agreement between typing methods | Python | `ex1_1_discrimination.py` |
| 2 | SNP distance, time and the molecular clock | Python and R | `ex1_2_snp_distance.py`, `ex1_2b_r_check.R` |
| 3 | design matters | Python | `ex1_3_design_simulation.py` |
| 4 | an application, lineage diversity and variant sweeps | Python | `ex1_4_lineage_diversity.py` |
| 5 | simulating the epidemic process behind a genomic data set | Python | `ex1_5_epidemic_process.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 1                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 1 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 1 --script ex1_1_discrimination.py   # run one exercise; the output stays in chapters/chapter01/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 1
./scripts/run-chapter.ps1 -Chapter 1 -Full
./scripts/run-chapter.ps1 -Chapter 1 -Script ex1_1_discrimination.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 1`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 93 of 93, in 3.1 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython (optional: anthropic)
- R packages: none
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.
- Note: Hands-on 6 and 7 of the agentic section run without a language model; ex1_7 uses the anthropic package and an API key only if you want a live model.

## Data

20 files, 4.9 MB. Largest: `ebola_aln.fasta` (3.8 MB), `lapis_lineage_monthly.csv` (1.0 MB), `typing_pfge_subset.csv.gz` (0.1 MB).
The result tables are `table_sim_sampling.csv, table_lineage_diversity.csv, table_snp_thresholds.csv, table_mock_agent_trace.csv, table_sim_steps.csv, table_design_simulation.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
