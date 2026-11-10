# 03. *R. australiae*: population genomics

Genetic diversity, differentiation, structure and isolation by distance and environment in the bottlenose wedgefish (*Rhynchobatus australiae*, RA), analysed separately for the four SNP datasets from [stage 2](02_seascape_outliers.md): neutral, adaptive (any method), high confidence and environmental (RDA).

Everything is run in R from one script, [`scripts/03_population_genomics.R`](../scripts/03_population_genomics.R), which starts by loading the objects saved at the end of stage 2 (`RA_seascape_FINAL_checkpoint.RData`).

## Missing genotypes

HO, HE, FIS and FST are calculated from genotypes **without imputation**, because `hierfstat` handles missing data itself, and genotypes are passed to `hierfstat` in its two-digit format (11, 12, 22). The sNMF analysis (section 7) uses a matrix in which missing genotypes are replaced by the most common genotype at each locus, as in stage 2.

## Populations

MZ (Mozambique, n = 28), TZ (Tanzania, 8), MD (Madagascar, 2), RI (Réunion Island, 9), SY (Seychelles, 12), MA (Malaysia, 10), AU (Australia, 2). Madagascar and Australia are excluded from population-level diversity statistics because of their small sample sizes.

## Sections

| Section | Analysis | Main settings | Main outputs |
|---|---|---|---|
| 1 | Diversity statistics: HO, HE, FIS | `hierfstat` without imputation; 95% confidence intervals for FIS from 500 bootstrap resamples of loci | `RA_set1_*_diversity_table.csv`, `RA_set1_diversity_results.RData` |
| 2 | Allelic richness, private alleles and weighted means | richness rarefied to the smallest sample (`hierfstat::allelic.richness`); private alleles counted as loci with `poppr::private_alleles`; weighted means of HO, HE, FIS and AR use sample size as weights | `RA_set1_AR_PA_all_datasets.csv`, `RA_set1_AR_PA_results.RData`, `RA_set1_diversity_weighted_means.csv` |
| 3 | Pairwise and global FST | Weir and Cockerham (`hierfstat`); permutation tests (9 999 permutations of individuals between the two populations, or all possible arrangements when fewer exist); Benjamini-Hochberg adjustment within each dataset; significance at alpha = 0.05 | `RA_set1_FST_results.RData`, `RA_set1_*_FST_permutation_table.csv`, `RA_set1_FST_permutation_all_datasets.csv` |
| 4 | FST heatmaps | neutral below and adaptive above the diagonal (main text); high confidence below and environmental above (supplement); an asterisk marks a significant pair and negative values are shown as zero | `RA_set1_FST_heatmap*.pdf/png` |
| 5 | DAPC | PCs explaining 80% of variance; `find.clusters` for K = 1 to 10 (BIC); `xvalDapc` (100 replicates, 90% training set) and `optim.a.score` (100 simulations); with and without population priors | `RA_set1_DAPC_summary_table.csv`, one results file per dataset |
| 6 | Adaptive DAPC without Australia | Australia (n = 2) removed because its extreme scores compress the rest of the plot | `RA_set1_adaptive_dapc_noAU_results.RData` |
| 7 | sNMF | `LEA`; K = 1 to 10, 10 replicates, alpha = 100, K chosen by cross-entropy; barplots with `pophelper` and colours fixed to geographic groups | `RA_set1_snmf_summary.csv`, barplots |
| 8 | Isolation by distance and environment | see below | `RA_set1_IBD_IBE_results.RData`, `RA_set1_IBD_IBE_results.csv`, plots |
| 9 | Diversity statistics figure (neutral against adaptive) | built from the section 1 and 2 CSV files | `diversity_neutral_vs_adaptive.pdf/png` |
| 10 | Sampling map | `sf` and `rnaturalearth` | `RA_IndoPacific_zoomed_out.png/pdf` |

## Input

Objects from stage 2 (`RA_seascape_FINAL_checkpoint.RData`) and the stage 1 VCF (`RApop_set1.snps.vcf`).
