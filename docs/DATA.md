# Data: sources, terms and files

The exercises use public data. Each chapter's `data/` folder holds the inputs the scripts need to run offline (some trimmed to the columns and rows used) and the result tables (`table_*.csv`). This page says where the data come from and what terms I know of. **It is not legal advice, and the redistribution status of each file has not been reviewed by a lawyer.** The licences of this repository (MIT for code, CC BY 4.0 for text) do not cover the data; see `LICENSING.md`. Check the source terms before you publish the repository, and remove a file if its terms forbid redistribution; the scripts can rebuild most inputs from the source (see each chapter's `PACKAGE_NOTES.txt`).

| Data | Source | Used in | Terms as I understand them |
|---|---|---|---|
| Ebola virus genomes, Sierra Leone 2014 (249 records, 195 analysed) | NCBI BioProject PRJNA257197 (GenBank) | Chapters 1, 3, 5, 6, 7, 13 | NCBI places no restrictions on use of GenBank data; submitters may claim rights (NCBI policies page) |
| Ebola line list, Sierra Leone (11,903 cases) | R package `outbreaks` 1.9.0, `ebola_sierraleone_2014` | Chapter 13 | package licence GPL (>= 2) |
| Hagelloch measles outbreak 1861 (188 cases) | R package `outbreaks` / the `outbreaker2` data | Chapter 7 | package licences: `outbreaks` GPL (>= 2), `outbreaker2` MIT |
| Salmonella and Listeria isolate metadata, assembly statistics and SNP clusters | NCBI Pathogen Detection (PDG000000002.4290 and others; PDG000000001.4985 for Listeria) | Chapters 1, 2, 3, 4, 9, 12 | as NCBI above |
| SARS-CoV-2 genomes for temporal signal; Zika, H3N2 and other virus genomes | NCBI Virus and GenBank | Chapters 5, 6, 8 | as NCBI above |
| Sequencing reads of four public runs (first reads only) | ENA / SRA | Chapter 3 | INSDC open data |
| SARS-CoV-2 sequence counts by lineage and clade | CoV-Spectrum open API (LAPIS) | Chapters 1, 2, 13 | aggregated counts; check the API's current terms, since the underlying sequences have their own |
| England daily cases by specimen date | UKHSA data dashboard API | Chapter 13 | Open Government Licence v3.0 (check the dashboard's current terms) |
| 1000 Genomes high-coverage genotypes, European samples | EBI / IGSR | Chapters 10, 11 | open; cite the 1000 Genomes Project |
| GWAS summary statistics and catalogue tables; COVID-19 Host Genetics Initiative release 5 | GWAS Catalog (EBI); COVID-19 HGI | Chapters 10, 11 | freely available for research with citation; check the current terms |
| AMRFinderPlus database 4.2/2026-10-29.1 and raw outputs for eight assemblies | NCBI | Chapter 9 | NCBI public data; the program is public domain |
| RGI (CARD) raw outputs for 12 assemblies (`chapters/chapter09/data/rgi_outputs`) | CARD, McMaster University | Chapter 9 | **Decision needed.** CARD allows free non-commercial use by individuals at academic institutions only; non-academic organisations need written permission and commercial use needs a licence (LICENSE file of the RGI repository, read 16 November 2026). A public repository cannot meet those terms for every downloader: remove these 12 files from the public repository or ask McMaster University |
| World Bank country classification and population (API data of 8 October 2026) | World Bank Data Catalog | Chapter 12 | the catalogue lists CC BY 4.0 among its standard licences; the licence of each indicator is on its own page (not checked indicator by indicator) |
| Natural Earth 1:50 million country boundaries | Natural Earth | Chapter 8 | public domain |
| Infertility case-control data (`infert`) | R `datasets` package | Chapter 2 | R is GPL |
| Numbers typed from published papers and teaching notes (audits of Li et al., Klaas et al., Kitamura et al., Khanna et al., the NARMS slide set, and worked examples in teaching notes) | the papers and documents named in the chapters | Chapters 8, 9, 10, 13 | only numbers are recomputed; no figures or text are copied. Some of the papers are CC BY-NC-ND or CC BY-NC-SA |

Population figures, income groups and similar reference tables used in Chapter 12 are described in that chapter and in `chapters/chapter12/PACKAGE_NOTES.txt`.

## Files larger than 1 MB, by chapter


- **Chapter 1**: `ebola_aln.fasta` (3.8 MB).
- **Chapter 2**: `asm_batch_enteritidis.csv.gz` (1.2 MB).
- **Chapter 3**: `salmonella_assembly_qc.csv.gz` (14.3 MB), `ebola_aln.fasta` (3.8 MB), `DRR090587_demo.fastq.gz` (3.3 MB), `SRR33434758_demo.fastq.gz` (2.8 MB), `DRR1060924_demo.fastq.gz` (1.8 MB), `SRR30474864_demo.fastq.gz` (1.3 MB).
- **Chapter 4**: `listeria_cluster_list.tsv` (4.4 MB), `listeria_all_isolates.tsv` (3.2 MB).
- **Chapter 5**: `ebola_aln.fasta` (3.8 MB).
- **Chapter 6**: `sars2.fasta` (13.2 MB), `ebola_aln.fasta` (3.8 MB).
- **Chapter 7**: `ebola.fasta` (4.8 MB).
- **Chapter 8**: `zika_aln.fasta` (7.0 MB), `zika_selected.fasta` (6.5 MB), `ne_50m_countries.geojson` (3.1 MB).
- **Chapter 9**: `salmonella_cluster_list.tsv` (20.8 MB), `sal_us_frame.csv.gz` (5.5 MB), `frame_for_r.csv.gz` (2.2 MB).
- **Chapter 10**: `chr22_group_geno_counts.tsv.gz` (1.2 MB).
- **Chapter 11**: `chr1_il6r_eur.gt` (14.1 MB).
- **Chapter 12**: `cluster_list_us.tsv` (16.5 MB), `sal_surv.csv.gz` (6.8 MB).
- **Chapter 13**: `ebola.fasta` (4.8 MB), `ebola_aln.fasta` (3.8 MB).
