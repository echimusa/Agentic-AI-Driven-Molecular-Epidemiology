# Chapter map: hands-on exercises and their scripts

Each chapter's README gives the same table with the commands to run it.

## Chapter 1: Principles, Designs, and Applications of Molecular Epidemiology

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | discrimination and agreement between typing methods | `ex1_1_discrimination.py` |
| 2 (Python and R) | SNP distance, time and the molecular clock | `ex1_2_snp_distance.py`, `ex1_2b_r_check.R` |
| 3 (Python) | design matters | `ex1_3_design_simulation.py` |
| 4 (Python) | an application, lineage diversity and variant sweeps | `ex1_4_lineage_diversity.py` |
| 5 (Python) | simulating the epidemic process behind a genomic data set | `ex1_5_epidemic_process.py` |

## Chapter 2: Epidemiological Study Design, Bias and Sampling for Molecular Studies

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python and R) | nested case–control sampling | `ex2_1_nested_case_control.py`, `ex2_1b_infert_clogit.R` |
| 2 (R and Python) | batch effects | `ex2_2a_batch_pca.R`, `ex2_2_batch_effects.py` |
| 3 (Python) | population stratification | `ex2_3_stratification.py` |
| 4 (Python) | sampling for molecular studies | `ex2_4_sampling.py` |

## Chapter 3: Molecular Typing and Sequencing Methods and Data Quality

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | assembly quality control on public genomes | `ex3_1_assembly_qc.py` |
| 2 (Python and R) | read-level quality control on public runs | `ex3_2_read_qc.py`, `ex3_2b_per_base_quality.R` |
| 3 (Python) | depth, error and SNP calling | `ex3_3_snp_calling_simulation.py` |
| 4 (Python) | does assembly quality change the typing result? | `ex3_4_quality_and_typing.py` |

## Chapter 4: Genomic Clustering: SNP and cgMLST Approaches

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | distance matrices, thresholds and linkage rules | `ex4_1_distance_matrices.py` |
| 2 (R) | allele distances, pairs and clusters | `ex4_2_cgmlst_alleles.R` |
| 3 (Python) | SNP distance against allele distance | `ex4_3_snp_vs_allele.py` |
| 4 (Python) | chaining, fragmentation and links in real clusters | `ex4_4_real_clusters.py` |

## Chapter 5: Phylogenetic Inference

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python and R) | distances and distance trees | `ex5_1_distances_nj.py`, `ex5_1b_nj_bootstrap_ape.R`, `ex5_1c_tree_likeness.py` |
| 2 (Python) | the likelihood from first principles | `ex5_2_likelihood.py` |
| 3 (Python) | support and accuracy | `ex5_3_bootstrap_accuracy.py`, `ex5_3d_ebola_full_support.py` |
| 4 (Python) | rooting and the temporal signal | `ex5_4_root_to_tip.py` |
| 5 (Python) | modelling in detail | `ex5_6a_extract_cds.py`, `ex5_6_modelling_detail.py` |
| 6 (Python) | trees from variable sites only | `ex5_7_ascertainment.py` |
| 7 (Python) | recombination, branch lengths and the tree | `ex5_8_recombination.py` |

## Chapter 6: Phylodynamics and Molecular Clocks

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | temporal signal in real SARS-CoV-2 genomes | `ex6_root_to_tip.py` |
| 2 (R) | skyline plots, checked against a known truth | `ex6_skyline.R` |
| 3 (Python) | time-scaled trees with TreeTime | `ex6_treetime.py` |
| 4 (Python) | time-dependent rates in a model | `ex6_time_dependent.py` |
| 5 (BEAST 2) | Bayesian dating, BETS and Bayes factors | `ex6_beast.py`, `run_all.sh`, `ex6_beast_analysis.py` |

## Chapter 7: Transmission Inference and Outbreak Investigation

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | who infected whom in a measles outbreak | `ex1_hagelloch.py` |
| 2 (Python) | what do genomes add? Ebola virus in Sierra Leone | `ex2_ebola.py` |
| 3 (Python) | describe the outbreak in time, place and person, and size a tracing window | `ex4_descriptive.py` |
| 4 (R) | Bayesian reconstruction with outbreaker2 | `ex3_outbreaker2.R` |
| 5 (R) | TransPhylo on simulated outbreaks with a known truth | `ex5_transphylo.R`, `ex5_transphylo_analysis.py` |

## Chapter 8: Phylogeography and Visualisation of Spread

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | from public Zika genomes to a map of introductions | `ex8_1_dated_tree.py`, `ex8_2_migration.py` |
| 2 (Python) | do sampling and the clock rate change the answer? | `ex8_1_dated_tree.py`, `ex8_3_sampling_bias.py`, `ex8_4_clock_sensitivity_figure.py`, `ex8_5_dating_reruns.py` |
| 3 (Python) | drawing where, how much and when | `ex8_2_migration.py`, `ex8_10_microreact_export.py` |
| 4 (Python) | diversity, differentiation and isolation by distance | `ex8_7_popgen.py` |
| 5 (Python) | which sampling strategy recovers the true source? | `ex8_6_sampling_truth.py` |

## Chapter 9: Antimicrobial Resistance Genomics in a One Health Context

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | from public isolates to a One Health sampling frame | `ex9_1_sampling_frame.py`, `ex9_7_metadata_completeness.py` |
| 2 (Python and R) | resistance determinants and multidrug resistance by sector | `ex9_2_resistance_by_sector.py`, `ex9_2b_glm.R` |
| 3 (Python) | genotype versus phenotype | `ex9_3_genotype_phenotype.py`, `ex9_3b_discordance.py` |
| 4 (Python) | cross-sector clusters and genes that travel together | `ex9_4_clusters_cooccurrence.py` |
| 5 (Python) | calling resistance genes from genomes with AMRFinderPlus | `ex9_5_amrfinder_genomes.py` |
| 6 (Python and R) | two detection tools on the same genomes | `ex9_5b_two_tools.py`, `ex9_5c_summarise.R` |
| 7 (Python) | how many cross-sector pairs survive each SNP threshold? | `ex9_6_snp_thresholds.py` |
| 8 (Python) | reading three cross-sector studies and recomputing their numbers | `ex9_8_study_audit.py` |

## Chapter 10: Host Genetic Epidemiology and Genome-Wide Association

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | what GWAS of infection have found, and in whom | `ex10_1_catalog.py` |
| 2 (Python and R) | quality control and loci from public summary statistics | `ex10_2_qc_manhattan.py`, `ex10_2b_r_check.R` |
| 3 (Python) | population stratification with a known truth | `ex10_3_stratification.py` |
| 4 (Python) | fine-mapping the 3p21.31 locus | `ex10_4_finemap.py` |
| 5 (Python and PLINK) | the quality-control command, ancestry and Hardy–Weinberg | `ex10_5_qc_structure.py` |
| 6 (Python) | how many people? | `ex10_6_power_lambda.py` |
| 7 (Python) | why a signal may not replicate, the 3p21.31 variants in five ancestries | `ex10_7_ancestry_3p21.py` |
| 8 (Python) | admixture mapping with a known truth | `ex10_8_admixture_mapping.py` |
| 9 (Python) | recomputing the worked examples, and a checklist for reading a GWAS paper | `ex10_9_notes_audit.py` |

## Chapter 11: Mendelian Randomisation and Causal Inference

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | body mass index and hospitalised COVID-19 | `ex11_1_bmi_covid.py` |
| 2 (Python and R) | robust estimators and diagnostics | `ex11_1_bmi_covid.py`, `ex11_2b_r_check.R` |
| 3 (Python) | a known truth | `ex11_3_simulation.py` |
| 4 (Python) | drug-target MR at IL6R | `ex11_4_cis_mr.py` |
| 5 (Python) | why alleles must be harmonised | `ex11_5_harmonise.py` |
| 6 (Python) | vitamin D and hospitalised COVID-19 | `ex11_6_vitamin_d.py` |
| 7 (Python) | reconciling an observational association with a null MR | `ex11_7_reconcile.py` |
| 8 (Python) | one sample, two samples and the winner’s curse | `ex11_8_one_two_sample.py` |
| 9 (Python) | a binary outcome, testing against estimating | `ex11_9_binary_outcome.py` |
| 10 (Python) | invalid instruments, the weighted mode and multivariable MR | `ex11_10_mode_mvmr.py` |
| 11 (Python) | a within-family design against dynastic effects | `ex11_11_family_design.py` |

## Chapter 12: Genomic Surveillance Systems, Ethics and Data Governance

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | completeness and timeliness from public metadata | `ex12_1_indicators.py` |
| 2 (Python and R) | representativeness and equity | `ex12_2_equity.py`, `ex12_2b_gini.R` |
| 3 (Python) | privacy risk and utility | `ex12_3_privacy.py` |
| 4 (Python) | how reporting delay postpones cluster detection | `ex12_4_detection.py` |
| 5 (Python) | how many sequences a week does an aim need? | `ex12_5_sizing.py` |
| 6 (Python) | where the sequences are taken decides when a variant is seen | `ex12_6_allocation.py` |
| 7 (Python) | what free text costs | `ex12_7_standards.py` |
| 8 (Python) | what an embargo costs | `ex12_8_embargo.py` |
| 9 (Python) | check a surveillance protocol | `ex12_9_protocol.py` |

## Chapter 13: Integrating Molecular Data with Transmission Models

| Hands-on | Exercise | Script |
|---|---|---|
| 1 (Python) | reproduction number and growth rate from public case data | `ex13_1_ebola_rt.py`, `ex13_1b_r_check.R` |
| 2 (Python) | what the genomes show, and what they cannot | `ex13_2_genome_growth.py` |
| 3 (Python) | a known truth, and what genomes add to case data | `ex13_3_simulation.py` |
| 4 (Python) | variant replacement and a two-variant renewal forecast | `ex13_4_variants.py` |
| 5 (Python) | the SIR model, and why R = exp(r T) is not enough | `ex13_5_sir_gi.py` |
| 6 (Python) | different windows, different answers | `ex13_6_windows.py` |
| 7 (Python) | a partly sampled outbreak with a host attribute | `ex13_7_partial_sampling.py` |
| 8 (Python) | two databases, one outbreak | `ex13_8_integration_audit.py` |
| 9 (Python) | check a capstone analysis plan and briefing | `ex13_9_capstone_check.py` |
