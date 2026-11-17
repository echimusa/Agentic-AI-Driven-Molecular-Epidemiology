# Heavy analyses and large downloads

The quick checks and the full reruns read the shipped results of these analyses; they do not repeat them. To repeat them yourself, use Ubuntu or WSL2.

| Chapter | Analysis | Needs | Time (one processor) |
|---|---|---|---|
| 6 | Hands-on 5: BEAST 2 chains and path sampling for six models on 60 Ebola genomes | BEAST 2.7.7 with the MODEL_SELECTION package 1.6.2 (`scripts/setup-heavy.ps1`), Java | about 90 minutes (`code/run_all.sh`) |
| 6 | Hands-on 3: TreeTime reruns for Ebola virus and SARS-CoV-2 | `phylo-treetime`, FastTree | about 25 minutes (`run-chapter.ps1 -Chapter 6 -Full -Sim`) |
| 7 | Hands-on 5: 12 TransPhylo runs | R packages `TransPhylo`, `coda` | about 40 minutes |
| 9 | Hands-on 5: AMRFinderPlus on eight assemblies | AMRFinderPlus 4.2.7, its database (release 4.2/2026-10-29.1, 43 files), BLAST+, HMMER | about 1.5 minutes per genome |
| 9 | Hands-on 6: RGI (CARD) on the same eight assemblies | RGI 6.0.8, DIAMOND, Prodigal | not timed |
| 9 | Hands-on 1: the NCBI Pathogen Detection Salmonella table | a download of 883 MB; the table keeps growing, so counts change | download time |
| 10 | Hands-on 2 and 3: harmonised COVID-19 HGI file (336 MB) and a chromosome 22 VCF (205 MB) | downloads from EBI | download time |

## BEAST 2 and AMRFinderPlus

```bash
./scripts/setup.sh --chapters 6,9 --heavy      # Java, BLAST, HMMER and the other tools
./scripts/setup-heavy.sh                       # BEAST 2.7.7, its MODEL_SELECTION package, and the AMRFinderPlus binaries
```

(PowerShell: `./scripts/setup.ps1 -Chapters 6,9 -Heavy` and `./scripts/setup-heavy.ps1`.)

`setup-heavy.sh` downloads BEAST 2 into `tools/linux/` and AMRFinderPlus into `amrf/`. The other scripts add `tools/linux/bin` to the `PATH` themselves; add it to your own shell with `export PATH=$PWD/tools/linux/bin:$PATH`. Then download the 43 files of the AMRFinderPlus database release 4.2/2026-10-29.1 from `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/` into `amrf/db` and run `amrf/amrfinder_index amrf/db`. The built-in updater failed in the environment used to prepare the chapter, which is why the files were fetched by hand.

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
