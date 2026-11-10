# HPC step (section 8 of 02_seascape_outliers.R)
# Partial RDA of genotypes on the four retained environmental variables, conditioned on two
# dbMEMs, with 999-permutation ANOVA tests and variance partitioning.
# Inputs (same directory): G_imp.rds, locus_names_full.rds, RA_coordinates.csv, RApopmap.txt,
#                          RA_dbmems.csv, RA_env_final_set1_11vars.csv
# Output: RA_set1_partialRDA_results_11.rds

# .libPaths("<path to R library>")
library(vegan)
library(readxl)

# setwd("<working directory containing the input files>")

cat("Loading data...\n")

# ---- Load objects ----
G_imp            <- readRDS("G_imp.rds")
locus_names_full <- readRDS("locus_names_full.rds")

coords  <- read.csv("RA_coordinates.csv")
popmap  <- read.table("RApopmap.txt", header=FALSE)
colnames(popmap) <- c("ID", "pop")

dbmem_df <- read.csv("RA_dbmems.csv")
env_df   <- read.csv("RA_env_final_set1_11vars.csv")

cat("G_imp dimensions:", dim(G_imp), "\n")
cat("Alignment check:", all(rownames(G_imp) == coords$ID), "\n")

# ---- Prepare env and spatial data ----
env_df   <- env_df[match(coords$ID, env_df$ID), ]
rownames(env_df) <- env_df$ID
env_only_final   <- env_df[, !colnames(env_df) %in% "ID"]

mem_df <- data.frame(
  MEM1 = dbmem_df$MEM1[match(coords$ID, dbmem_df$ID)],
  MEM2 = dbmem_df$MEM2[match(coords$ID, dbmem_df$ID)],
  row.names = coords$ID
)

cat("Env variables:", colnames(env_only_final), "\n")
cat("Start time:", format(Sys.time()), "\n")

# ---- Partial RDA ----
cat("Running partial RDA...\n")
rda_partial <- rda(G_imp ~ . + Condition(as.matrix(mem_df)),
                   data = env_only_final)

# ---- Significance testing ----
cat("Running anova...\n")
anova_partial      <- anova(rda_partial, permutations=999)
anova_partial_axes <- anova(rda_partial, by="axis", permutations=999)

cat("Partial RDA results:\n")
print(anova_partial)
print(anova_partial_axes)

# ---- Variance partitioning ----
cat("Running varpart...\n")
vp <- varpart(G_imp, env_only_final, mem_df)
print(vp)

# Additional models for full variance partitioning reporting
rda_env_only <- rda(G_imp ~ ., data = env_only_final)
anova_env    <- anova(rda_env_only, permutations = 999)
cat("Env-only adj R2:",
    RsquareAdj(rda_env_only)$adj.r.squared, "\n")
print(anova_env)

rda_sp_only <- rda(G_imp ~ ., data = mem_df)
anova_sp    <- anova(rda_sp_only, permutations = 999)
cat("Space-only adj R2:",
    RsquareAdj(rda_sp_only)$adj.r.squared, "\n")
print(anova_sp)

rda_full   <- rda(G_imp ~ ., data = data.frame(env_only_final, mem_df))
anova_full <- anova(rda_full, permutations = 999)
cat("Full model adj R2:",
    RsquareAdj(rda_full)$adj.r.squared, "\n")
print(anova_full)

cat("End time:", format(Sys.time()), "\n")

# Save all results
saveRDS(list(
  rda_partial        = rda_partial,
  anova_partial      = anova_partial,
  anova_partial_axes = anova_partial_axes,
  vp                 = vp,
  env_only_final     = env_only_final,
  mem_df             = mem_df,
  rda_env_only       = rda_env_only,
  anova_env          = anova_env,
  rda_sp_only        = rda_sp_only,
  anova_sp           = anova_sp,
  rda_full           = rda_full,
  anova_full         = anova_full
), file = "RA_set1_partialRDA_results_11.rds")

cat("Done — results saved to RA_set1_partialRDA_results_11.rds\n")