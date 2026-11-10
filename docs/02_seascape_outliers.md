# 02. *R. australiae*: seascape genomics and outlier detection

Detection of loci under selection in the bottlenose wedgefish (*Rhynchobatus australiae*, RA) with structure-based methods (OutFLANK, pcadapt) and genotype-environment association methods (partial RDA, LFMM), followed by compilation of four SNP datasets used in the population genomic analyses ([stage 3](03_population_genomics.md)).

The analysis is one R script, [`scripts/02_seascape_outliers.R`](../scripts/02_seascape_outliers.R), run in sections. Two steps ran on a cluster (**[HPC]**) because they are computationally heavy, and the rest in R (**[R]**). Parameters are defined once at the top of the script.

## Input

| File | Contents |
|---|---|
| `RApop_set1.snps.vcf` | RA_set1: 71 individuals, 27 562 SNPs (output of [stage 1](01_assembly_and_filtering.md)) |
| `RApopmap.txt` | two columns: individual ID, population code |
| `RA_coordinates.csv` | ID, latitude, longitude |

Populations: MZ (Mozambique, n = 28), TZ (Tanzania, 8), MD (Madagascar, 2), RI (Réunion Island, 9), SY (Seychelles, 12), MA (Malaysia, 10), AU (Australia, 2).

## Sections

| Section | What it does | Where |
|---|---|---|
| 1 | Load VCF, popmap and coordinates; build genotype matrix (0/1/2), replace missing genotypes with the most common genotype at each locus, build genind | R |
| 2 | OutFLANK | R |
| 3 | pcadapt (K from scree plot and genomic inflation factor), Storey q-values | R |
| 4 | Extract 11 environmental variables (Bio-ORACLE v3, MARSPEC) at sampling coordinates | R |
| 5 | Marine least-cost distances and distance-based Moran's eigenvector maps (dbMEMs) | R |
| 6 | Prepare input and run forward-backward selection of environmental variables (`ordistep`) | R, then **HPC** |
| 7 | Variance inflation factor (VIF) filtering of environmental variables and MEMs | R |
| 8 | Partial RDA (environment conditioned on space), 999-permutation ANOVA, variance partitioning | **HPC**, then R |
| 9 | RDA outlier loci (±3 SD on significant axes), correlation of candidate loci with each variable | R |
| 10 | LFMM (`lfmm2`) on environmental PCs, FDR < 0.01 | R |
| 11 | Compile the four datasets | R |
| 12 | Venn diagram of outlier sets | R |
| 13 | RDA triplots | R |
| 14 | Top 10 candidate loci | R |

### HPC steps

Sections 6 and 8 write input files, run a job on the cluster (R 4.3.2) and load the result back into R. The cluster scripts are in [`scripts/`](../scripts/):

| Step | Script | Inputs | Output |
|---|---|---|---|
| Variable selection (section 6) | `set1_RDA_ordistep11.R` | `G_imp.rds`, `locus_names_full.rds`, `RA_coordinates.csv`, `RA_environmental_data_11vars.csv` | `RA_set1_ordistep_results_11.rds` |
| Partial RDA (section 8) | `set1_partialRDA.R` | `G_imp.rds`, `locus_names_full.rds`, `RA_coordinates.csv`, `RApopmap.txt`, `RA_dbmems.csv`, `RA_env_final_set1_11vars.csv` | `RA_set1_partialRDA_results_11.rds` |

`RA_env_final_set1_11vars.csv` is written in section 7 after VIF filtering, so the partial RDA is submitted after that section.

## Key parameters

| Parameter | Value |
|---|---|
| OutFLANK | trim fractions 0.05, q < 0.05, Hmin 0.1 (no outliers detected, so not used further) |
| pcadapt | K = 4; Storey q < 0.05 |
| RDA | partial RDA on 4 retained variables, conditioned on 2 dbMEMs; candidates beyond 3 SD; 999 permutations |
| LFMM | K = 4 latent factors; 3 environmental PCs (broken stick); FDR < 0.01 |

## The four datasets

Section 11 builds four SNP sets from the three outlier methods (pcadapt, RDA, LFMM) and saves them in `RA_set1_four_datasets_v2.RData`:

| Dataset | Definition | SNPs |
|---|---|---|
| Neutral | all SNPs not detected by any method | 24 415 |
| Adaptive | detected by any of the three methods | 3 147 |
| High confidence | detected by two or more methods | 720 |
| Environmental | RDA outliers | 846 |

## Output

Main objects are saved at the end of the script in `RA_seascape_FINAL_checkpoint.RData`, and the four datasets in `RA_set1_four_datasets_v2.RData`. These are loaded in stage 3.
