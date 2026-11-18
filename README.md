# Agentic AI-driven Molecular Epidemiology: code, data and teaching material

Code, data and checks for the 13 chapters of *Agentic AI-driven Molecular Epidemiology* (Emile Rugamika Chimusa). The book has 85 hands-on exercises on public data; every one runs from the files in this repository, and a check script recomputes the numbers printed in each chapter.

## Quick start

Ubuntu 22.04 or 24.04 (or WSL2 on Windows), with Bash:

```bash
# get git installed if it is not yet by sudo apt-get install git
git clone https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology.git
cd Agentic-AI-Driven-Molecular-Epidemiology
./scripts/setup.sh --chapters 12,13          # or leave out --chapters to prepare all 13
./scripts/check-install.sh --chapters 12,13
./scripts/run-chapter.sh 13
```

Every script also exists for PowerShell (Windows, or Ubuntu with `pwsh`): `scripts/setup.ps1 -Chapters 12,13`, `scripts/check-install.ps1`, `scripts/run-chapter.ps1 -Chapter 13`, `get git install on PowerShell windows with winget install --id Git.Git -e --source winget`.
If Powershell flag error on policies, to clear the flag on all of them at once `Get-ChildItem -Recurse -Filter *.ps1 | Unblock-File`
If Powershell still flag error do `powershell -ExecutionPolicy Bypass -File .\scripts\setup.ps1`
Windows without WSL2: see [docs/INSTALL-Windows.md](docs/INSTALL-Windows.md) (chapters 5, 6 and 7 need WSL2). Ubuntu in detail: [docs/INSTALL-Ubuntu.md](docs/INSTALL-Ubuntu.md).

## The chapters

| Ch | Title | Hands-on | Tools needed | Full rerun | Data |
|---|---|---|---|---|---|
| [1](chapters/chapter01/README.md) | Principles, Designs, and Applications of Molecular Epidemiology | 5 | none | 3.1 min | 5 MB |
| [2](chapters/chapter02/README.md) | Epidemiological Study Design, Bias and Sampling for Molecular Studies | 4 | none | 4.2 min | 2 MB |
| [3](chapters/chapter03/README.md) | Molecular Typing and Sequencing Methods and Data Quality | 4 | none | 2.9 min | 27 MB |
| [4](chapters/chapter04/README.md) | Genomic Clustering: SNP and cgMLST Approaches | 4 | none | 10 s | 8 MB |
| [5](chapters/chapter05/README.md) | Phylogenetic Inference | 7 | iqtree2, FastTree | 5.0 min | 5 MB |
| [6](chapters/chapter06/README.md) | Phylodynamics and Molecular Clocks | 5 | none | 20 s | 34 MB |
| [7](chapters/chapter07/README.md) | Transmission Inference and Outbreak Investigation | 5 | none | 41 s | 5 MB |
| [8](chapters/chapter08/README.md) | Phylogeography and Visualisation of Spread | 5 | none | 1.9 min | 17 MB |
| [9](chapters/chapter09/README.md) | Antimicrobial Resistance Genomics in a One Health Context | 8 | none | 95 s | 38 MB |
| [10](chapters/chapter10/README.md) | Host Genetic Epidemiology and Genome-Wide Association | 9 | none | 30 s | 2 MB |
| [11](chapters/chapter11/README.md) | Mendelian Randomisation and Causal Inference | 11 | none | 1.8 min | 15 MB |
| [12](chapters/chapter12/README.md) | Genomic Surveillance Systems, Ethics and Data Governance | 9 | none | 48 s | 23 MB |
| [13](chapters/chapter13/README.md) | Integrating Molecular Data with Transmission Models | 9 | none | 24 s | 10 MB |

"Full rerun" is the time of `run-chapter.sh N --full` (`run-chapter.ps1 -Chapter N -Full`) on the reference machine (one processor, no network): it reruns every script and compares every result table. The quick check without `-Full` recomputes the numbers of the chapter and is much shorter (seconds to under half a minute). "Tools needed" are the command-line tools beyond Python, R packages and pandoc.

## Layout

```
chapters/chapterNN/   code/ data/ figures/ review/ README.md PACKAGE_NOTES.txt
bibliography/         one reference list for the whole book, Harvard style (Markdown, CSV)
docs/                 INSTALL-Ubuntu.md  INSTALL-Windows.md  TOOLS.md  DATA.md  SOURCE-LICENCES.csv  CHAPTER-MAP.md  HEAVY-ANALYSES.md  TROUBLESHOOTING.md  PUBLISHING.md  SOURCE-LICENCES.csv
scripts/              setup  check-install  run-chapter  run-all  setup-heavy, each as .sh (Bash) and .ps1 (PowerShell)
tools/                run_chapter.py  needed.py  check_install.py (used by the scripts)  chapter_requirements.json
requirements.txt      pinned Python packages      r-packages.R   R packages
```

## What was tested, and what was not

- Reference machine: Ubuntu 24.04 (Python 3.12.3, R 4.3.3, IQ-TREE 2.0.7, FastTree 2.1.11, MAFFT 7.505, bcftools 1.19, PLINK 1.90b7.2, pandoc 3.1.3). All 13 chapter checks pass there. The Bash scripts were run on Ubuntu 24.04 (`check-install`, `run-chapter` for chapters and for single exercises, `run-all` with a passing and a failing chapter, `setup` as a dry run) and pass shellcheck. The PowerShell scripts were run with PowerShell 7.4 on the same machine (`check-install`, `run-chapter`, `run-all`, `setup` as a dry run).
- Installation: `scripts/setup.sh --chapters 12 --skip-system` was run into a fresh virtual environment (the pinned packages install, `check-install.sh` reports them and Chapter 12 passes with that environment). The `apt` and CRAN steps were only run as a dry run, so the first full `setup.sh` or `setup.ps1` run on a new machine may need the fixes in [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md).
- Not tested at all: Windows (native or WSL2). The Windows instructions follow the documentation of winget, WSL and the packages.
- Result tables are compared byte for byte; a different version of IQ-TREE, TreeTime or a numerical library can change the last digits. The reference versions are in `tools/chapter_requirements.json`.

## Licence and citation

- Code: MIT ([LICENSE](LICENSE)). Text, figures, tables and bibliography: CC BY 4.0 ([LICENSE-CC-BY-4.0.txt](LICENSE-CC-BY-4.0.txt)). Data keep the terms of their sources ([docs/DATA.md](docs/DATA.md)); the licences of every source are in [docs/SOURCE-LICENCES.csv](docs/SOURCE-LICENCES.csv). See [LICENSING.md](LICENSING.md).
- Cite as in `CITATION.cff`. Repository: https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology
- The repository is about 222 MB. No file is above 50 MB (the Chapter 10 genotype extract is stored compressed, 0.7 MB); the largest is `chapters/chapter11/data/chr1_il6r_eur.gt` (14 MB).

contact : 
Prof. Emile R. Chimusa
Email: emile.chimusa@northumbria.ac.uk or echimusa@gmail.com
WhatApp: +353 89 262 0356
