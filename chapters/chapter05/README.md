# Chapter 5: Phylogenetic Inference

Code, data and checks for the exercises of Chapter 5 of *Agentic AI-driven Molecular Epidemiology*. Everything here runs from the files in `data/`; no network connection is needed for the quick check.

## Folders

| Folder | Contents |
|---|---|
| `code/` | the exercise scripts (Python, with R where a hands-on has an R version) |
| `data/` | inputs and the result tables (`table_*.csv`): 71 files, 4.9 MB |
| `figures/` | the figures of the chapter |
| `review/` | the chapter checks (`validate_ch5.py`) and their results |
| `PACKAGE_NOTES.txt` | the original notes of the package: provenance of the data, expected values, limits |

## Hands-on exercises

| # | Exercise | Language | Script |
|---|---|---|---|
| 1 | distances and distance trees | Python and R | `ex5_1_distances_nj.py`, `ex5_1b_nj_bootstrap_ape.R`, `ex5_1c_tree_likeness.py` |
| 2 | the likelihood from first principles | Python | `ex5_2_likelihood.py` |
| 3 | support and accuracy | Python | `ex5_3_bootstrap_accuracy.py`, `ex5_3d_ebola_full_support.py` |
| 4 | rooting and the temporal signal | Python | `ex5_4_root_to_tip.py` |
| 5 | modelling in detail | Python | `ex5_6a_extract_cds.py`, `ex5_6_modelling_detail.py` |
| 6 | trees from variable sites only | Python | `ex5_7_ascertainment.py` |
| 7 | recombination, branch lengths and the tree | Python | `ex5_8_recombination.py` |

## Run it

From the repository root, after `scripts/setup.sh` or `scripts/setup.ps1` (see `docs/INSTALL-Ubuntu.md` or `docs/INSTALL-Windows.md`). On Ubuntu, with Bash:

```bash
./scripts/run-chapter.sh 5                      # recompute every number of the chapter and compare it with the text
./scripts/run-chapter.sh 5 --full               # also rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 5 --script ex5_1_distances_nj.py   # run one exercise; the output stays in chapters/chapter05/.work
```

With PowerShell (Windows, or Ubuntu with `pwsh`):

```powershell
./scripts/run-chapter.ps1 -Chapter 5
./scripts/run-chapter.ps1 -Chapter 5 -Full
./scripts/run-chapter.ps1 -Chapter 5 -Script ex5_1_distances_nj.py
```

On Windows PowerShell 5.1 write `.\scripts\run-chapter.ps1`.

## What the check covers

- Quick check (`-Chapter 5`): recomputes the numbers of the chapter with code that does not import the chapter scripts and compares them with the results shipped in `data/`. The comparison with the text of the chapter needs your own copy of the Word file (`--docx FILE`; the Word files are not in this repository), otherwise those checks are skipped.
- Full rerun (`-Full`): copies code and data to a scratch folder, reruns every script offline and compares each result table with the shipped one byte for byte. Last run on the reference machine: 67 of 67, in 5.0 min.

## Tools

- Python packages: numpy, pandas, scipy, matplotlib, python-docx, python-pptx, lxml, biopython
- R packages: ape, ggplot2
- Command-line tools needed for the quick check and the full rerun: iqtree2, FastTree
- Tools needed only to rebuild inputs from public sources: none
- Heavy analyses (optional): none
- Also needed: pandoc, only if you compare the chapter text with `--docx`.

## Data

71 files, 4.9 MB. Largest: `ebola_aln.fasta` (3.8 MB), `h3n2_ha_aln.fasta` (0.1 MB), `h3n2_ha_raw.fasta` (0.1 MB).
The result tables are `table_codon_positions.csv, table_adequacy.csv, table_recombination.csv, table_recombination_raw.csv, table_modelfinder.csv, table_tree_scores.csv` and others.
Sources, access dates and terms of use are in `docs/DATA.md` and in `PACKAGE_NOTES.txt`. To rebuild an input from its public source, delete the file and run the exercise that creates it with a network connection; public data change, so counts may differ slightly.
