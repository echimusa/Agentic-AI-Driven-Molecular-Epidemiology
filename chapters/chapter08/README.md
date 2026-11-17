# Chapter 8: Phylogeography and Visualisation of Spread

Code, data and checks for the exercises of Chapter 8 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 27 files, 17.3 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch8.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | from public Zika genomes to a map of introductions | Python | `ex8_1_dated_tree.py`, `ex8_2_migration.py` |
| 2 | do sampling and the clock rate change the answer? | Python | `ex8_1_dated_tree.py`, `ex8_3_sampling_bias.py`, `ex8_4_clock_sensitivity_figure.py`, `ex8_5_dating_reruns.py` |
| 3 | drawing where, how much and when | Python | `ex8_2_migration.py`, `ex8_10_microreact_export.py` |
| 4 | diversity, differentiation and isolation by distance | Python | `ex8_7_popgen.py` |
| 5 | which sampling strategy recovers the true source? | Python | `ex8_6_sampling_truth.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 8                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 8 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 8 --script ex8_1_dated_tree.py   # run one exercise; the output stays in chapters/chapter08/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 8
./scripts/run-chapter.ps1 -Chapter 8 -Full
./scripts/run-chapter.ps1 -Chapter 8 -Script ex8_1_dated_tree.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 8`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 47 of 47, in 1.9 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython, phylo-treetime, shapely
- R packages: none
- Command-line tools needed for the quick check and the full rerun: none
- Tools needed only to rebuild inputs from public sources: mafft, FastTree
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

27 files, 17.3 MB. Largest: `zika_aln.fasta` (7.0 MB), `zika_selected.fasta` (6.5 MB), `ne_50m_countries.geojson` (3.1 MB).
The result tables are `table_tb_chi2.csv, table_popgen_groups.csv, table_sampling_truth_summary.csv, table_dating_reruns_summary.csv, table_popgen_ibd.csv, table_popgen_regions.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
