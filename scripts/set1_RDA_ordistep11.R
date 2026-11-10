# HPC step (section 6 of 02_seascape_outliers.R)
# Bidirectional stepwise selection of environmental variables with vegan::ordistep
# (11 starting variables, 500 permutations).
# Inputs (same directory): G_imp.rds, locus_names_full.rds, RA_coordinates.csv,
#                          RA_environmental_data_11vars.csv
# Output: RA_set1_ordistep_results_11.rds

# .libPaths("<path to R library>")
library(vegan)

# setwd("<working directory containing the input files>")
cat("Loading data...\n")

G_imp            <- readRDS("G_imp.rds")
locus_names_full <- readRDS("locus_names_full.rds")
coords <- read.csv("RA_coordinates.csv")
env_df <- read.csv("RA_environmental_data_11vars.csv")

cat("G_imp dimensions:", dim(G_imp), "\n")
cat("Alignment check:", all(rownames(G_imp) == coords$ID), "\n")

env_df <- env_df[match(coords$ID, env_df$ID), ]
rownames(env_df) <- env_df$ID
env_vars <- env_df[, !colnames(env_df) %in% "ID"]
rownames(env_vars) <- coords$ID

cat("Number of env variables:", ncol(env_vars), "\n")
cat("Env variables:", colnames(env_vars), "\n")
cat("Start time:", format(Sys.time()), "\n")

cat("Running ordistep on 11 env variables...\n")
rda_null <- rda(G_imp ~ 1, data=as.data.frame(env_vars))
rda_full <- rda(G_imp ~ ., data=as.data.frame(env_vars))

rda_step <- ordistep(
  rda_null,
  scope        = formula(rda_full),
  direction    = "both",
  permutations = 500,
  trace        = TRUE
)

selected_env <- attr(rda_step$terms, "term.labels")
cat("\nSelected env variables:\n")
print(selected_env)
cat("Number selected:", length(selected_env), "\n")
cat("End time:", format(Sys.time()), "\n")

saveRDS(list(
  selected_env = selected_env,
  rda_step     = rda_step
), file="RA_set1_ordistep_results_11.rds")

cat("Done — results saved to RA_set1_ordistep_results_11.rds\n")
