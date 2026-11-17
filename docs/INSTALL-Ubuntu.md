# Install on Ubuntu (22.04 or 24.04)

Tested on Ubuntu 24.04 with Python 3.12 and R 4.3.3. The same steps work in WSL2 on Windows (see [INSTALL-Windows.md](INSTALL-Windows.md)). Everything is installed once; afterwards each chapter runs from the files in this repository.

The scripts come in two versions that do the same job: **Bash** (`.sh`, used below; nothing else to install) and **PowerShell** (`.ps1`, in the last section, for people who work in PowerShell).

## What you need

- About 3 GB of free disk space (the repository is about 260 MB; the Python environment and R packages take the rest) and 4 GB of memory (the most demanding check, Chapter 9, peaked at about 2.5 GB).
- A user with `sudo` rights, for `apt`.
- A network connection for the installation. The exercises themselves run offline from the shipped data.

## 1. Get the repository

```bash
sudo apt-get update && sudo apt-get install -y git
git clone https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology.git
cd Agentic-AI-Driven-Molecular-Epidemiology
```

## 2. Install the tools

```bash
./scripts/setup.sh                      # all 13 chapters
./scripts/setup.sh --chapters 12,13     # only what these chapters need (ranges work too: 1-5)
./scripts/setup.sh --rebuild            # also MAFFT, FastTree, bcftools, tabix, PLINK, to rebuild inputs from public sources
./scripts/setup.sh --heavy              # also Java, BLAST and HMMER for BEAST 2 and AMRFinderPlus (then ./scripts/setup-heavy.sh)
./scripts/setup.sh --dry-run            # print every step and run none
./scripts/setup.sh --help
```

If a script says "Permission denied", run it as `bash scripts/setup.sh` or make the scripts executable with `chmod +x scripts/*.sh`.

The script installs with `apt` (Python, R, pandoc, compilers, and IQ-TREE 2 and FastTree for Chapter 5), creates a Python environment in `.venv` with the pinned versions of `requirements.txt`, and installs the R packages into `R/library` inside the repository, so nothing else on the machine is changed. R packages are compiled from CRAN and the first run takes several minutes.

If you already have Python 3.12, R and pandoc, add `--skip-system`.

## 3. Check the installation

```bash
./scripts/check-install.sh --chapters 12,13
```

It lists every package and tool with its version, marks any that differ from the tested version, and says which chapters are ready to run. It exits with a non-zero code if a chapter is not ready.

## 4. Run

```bash
./scripts/run-chapter.sh 13                                 # recompute the numbers of the chapter and compare them with the text
./scripts/run-chapter.sh 13 --full                          # rerun every script and compare each result table byte for byte
./scripts/run-chapter.sh 13 --script ex13_5_sir_gi.py       # one exercise; output in chapters/chapter13/.work
./scripts/run-chapter.sh 6 --full --sim                     # add the long TreeTime reruns of Chapter 6
./scripts/run-all.sh --chapters 1,2,4,12,13                 # several chapters and a summary table
```

Other options of `run-chapter.sh`: `--doi` (resolve every DOI in Crossref, needs the network) and `--offline` (block the network, to show that the run needs only the shipped data).

Typical times on one processor are in each chapter's README and in the root README. Chapters 1 to 4 and 6 to 13 need no command-line tools beyond pandoc; Chapter 5 needs IQ-TREE 2 and FastTree.

## 5. Rebuild from public sources (optional)

Each exercise reads its inputs from `data/`. To rebuild one from its public source, run `./scripts/setup.sh --rebuild`, delete the file, and run the exercise that creates it (the scripts say which in their first lines and in `PACKAGE_NOTES.txt`). Public data grow and change, so counts will differ a little from the book.

## 6. Remove everything

Delete the repository folder: the Python environment (`.venv`), the R packages (`R/library`) and the work folders (`.work`) are all inside it. The `apt` packages stay installed.

## PowerShell instead of Bash

The PowerShell scripts take the same options in PowerShell form. Install PowerShell 7 first:

```bash
sudo apt-get install -y wget apt-transport-https software-properties-common
wget -q "https://packages.microsoft.com/config/ubuntu/$(lsb_release -rs)/packages-microsoft-prod.deb"
sudo dpkg -i packages-microsoft-prod.deb && sudo apt-get update && sudo apt-get install -y powershell
```

(or `sudo snap install powershell --classic`), then start `pwsh` in the repository folder and use:

| Bash | PowerShell |
|---|---|
| `./scripts/setup.sh --chapters 12,13 --rebuild` | `./scripts/setup.ps1 -Chapters 12,13 -Rebuild` |
| `./scripts/setup.sh --skip-system --dry-run` | `./scripts/setup.ps1 -SkipSystem -DryRun` |
| `./scripts/check-install.sh --chapters 12,13` | `./scripts/check-install.ps1 -Chapters 12,13` |
| `./scripts/run-chapter.sh 13 --full` | `./scripts/run-chapter.ps1 -Chapter 13 -Full` |
| `./scripts/run-chapter.sh 13 --script ex13_5_sir_gi.py` | `./scripts/run-chapter.ps1 -Chapter 13 -Script ex13_5_sir_gi.py` |
| `./scripts/run-all.sh --chapters 1,2,4,12,13` | `./scripts/run-all.ps1 -Chapters 1,2,4,12,13` |
| `./scripts/setup-heavy.sh` | `./scripts/setup-heavy.ps1` |
