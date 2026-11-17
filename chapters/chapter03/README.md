# Chapter 3: Molecular Typing and Sequencing Methods and Data Quality

Code, data and checks for the exercises of Chapter 3 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 19 files, 27.4 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch3.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | assembly quality control on public genomes | Python | `ex3_1_assembly_qc.py` |
| 2 | read-level quality control on public runs | Python and R | `ex3_2_read_qc.py`, `ex3_2b_per_base_quality.R` |
| 3 | depth, error and SNP calling | Python | `ex3_3_snp_calling_simulation.py` |
| 4 | does assembly quality change the typing result? | Python | `ex3_4_quality_and_typing.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 3                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 3 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 3 --script ex3_1_assembly_qc.py   # run one exercise; the output stays in chapters/chapter03/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 3
./scripts/run-chapter.ps1 -Chapter 3 -Full
./scripts/run-chapter.ps1 -Chapter 3 -Script ex3_1_assembly_qc.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 3`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 56 of 56, in 2.9 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, scikit-learn, statsmodels
- R packages: ggplot2
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

19 files, 27.4 MB. Largest: `salmonella_assembly_qc.csv.gz` (14.3 MB), `ebola_aln.fasta` (3.8 MB), `DRR090587_demo.fastq.gz` (3.3 MB).
The result tables are `table_assembly_qc.csv, table_assembly_qc_labs.csv, table_quality_typing.csv, table_read_qc.csv, table_snp_calling.csv, table_snp_threshold.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
