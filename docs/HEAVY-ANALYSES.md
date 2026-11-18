# Heavy analyses and large downloads

The quick checks and the full reruns read the shipped results of these analyses; they do not repeat them. To repeat them yourself, use Ubuntu or WSL2.

| Chapter | Analysis | Needs | Time (one processor) |
|---|---|---|---|
| 6 | Hands-on 5: BEAST 2 chains and path sampling for six models on 60 Ebola genomes | BEAST 2.7.7 with the MODEL_SELECTION package 1.6.2 (`scripts/setup.sh --heavy`; BEAST 2 brings its own Java 17) | about 90 minutes (`code/run_all.sh`) |
| 6 | Hands-on 3: TreeTime reruns for Ebola virus and SARS-CoV-2 | `phylo-treetime`, FastTree | about 25 minutes (`run-chapter.ps1 -Chapter 6 -Full -Sim`) |
| 7 | Hands-on 5: 12 TransPhylo runs | R packages `TransPhylo`, `coda` | about 40 minutes |
| 9 | Hands-on 5: AMRFinderPlus on eight assemblies | AMRFinderPlus 4.2.7, its database (release 4.2/2026-10-29.1, 43 files), BLAST+, HMMER | about 1.5 minutes per genome |
| 9 | Hands-on 6: RGI (CARD) on the same eight assemblies | RGI 6.0.8, DIAMOND, Prodigal | not timed |
| 9 | Hands-on 1: the NCBI Pathogen Detection Salmonella table | a download of 883 MB; the table keeps growing, so counts change | download time |
| 10 | Hands-on 2 and 3: harmonised COVID-19 HGI file (336 MB) and a chromosome 22 VCF (205 MB) | downloads from EBI | download time |

## BEAST 2, Java and AMRFinderPlus: one command

```bash
./scripts/setup.sh --chapters 6,9 --heavy
```

(PowerShell: `./scripts/setup.ps1 -Chapters 6,9 -Heavy`; on Windows this part runs inside WSL2 Ubuntu.)

After the usual setup, this runs `scripts/setup-heavy.sh`, which checks each piece and installs only what is missing, so it is safe to run again:

| Piece | Checked how | If missing |
|---|---|---|
| BEAST 2.7.7 (Chapter 6) | `tools/linux/bin/beast` or `beast` on your `PATH` reports version 2.7.7 | downloads the Linux build (95 MB) into `tools/linux/` |
| Java | BEAST 2 includes its own Java 17 (`tools/linux/beast/jre`) | only if that is absent: a system Java 17 or newer is checked, and `openjdk-17-jre-headless` is installed with `apt` |
| MODEL_SELECTION 1.6.2 (Chapter 6, path sampling) | `~/.beast/2.7/MODEL_SELECTION/version.xml` | downloads the package into `~/.beast/2.7` |
| AMRFinderPlus 4.2.7 (Chapter 9) | `amrf/amrfinder --version` | downloads the binaries (6 MB) into `amrf/` |
| BLAST+ and HMMER | `makeblastdb` and `hmmpress` exist | installed with `apt` (they build the database index) |
| AMRFinderPlus database 4.2/2026-10-29.1 | all 43 files in `amrf/db`, `version.txt` says 2026-10-29.1, and the index was built | downloads the missing files from NCBI (132 MB) and runs `amrf/amrfinder_index amrf/db` |

`scripts/setup-heavy.sh --chapters 6` (or `9`) prepares one chapter only; `--force` installs again; `--dry-run` only reports. It needs about 360 MB of free disk space for BEAST 2 and about 300 MB for AMRFinderPlus, an x86-64 Linux system, and a network connection. `scripts/check-install.sh --chapters 6,9` lists what is installed. The other scripts add `tools/linux/bin` to the `PATH` themselves; add it to your own shell with `export PATH=$PWD/tools/linux/bin:$PATH`. `run-chapter.sh 9 --script ex9_5_amrfinder_genomes.py` finds `amrf/` for you.

The AMRFinderPlus database is pinned to release 2026-10-29.1 so that the results of the chapter can be reproduced; the built-in updater (`amrfinder -u`) would fetch a newer release and the counts could change.

For BEAST 2 in Chapter 6:

```bash
cd chapters/chapter06/code
python ex6_beast.py make --mcmc-length 6000000 --ps-steps 20 --ps-length 250000 --ps-pre 100000
bash run_all.sh
python ex6_beast_analysis.py
```

Chains with the same seeds should give the same logs on the same BEAST 2 version; other versions differ.

## RGI

RGI is not on PyPI under a name that installs. Clone `https://github.com/arpcard/rgi` into a separate virtual environment (RGI 6.0.8 worked with pandas 2.2 and numpy 2.2), install `diamond-aligner` and `prodigal` with `apt`, load the CARD data as the RGI documentation describes, and run `chapters/chapter09/code/ex9_5b_two_tools.py`.

## Large downloads

The scripts that download data say so in their first lines. Run them from `chapters/chapterNN/.work` (created by `run-chapter.sh --script`), with a network connection, after deleting the cached input they replace.
