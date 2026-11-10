# ======================================================================
# Population genomics: diversity, differentiation, structure and IBD/IBE
# Species: Rhynchobatus australiae (RA set1: 71 individuals, 27 562 SNPs)
# Four SNP datasets from stage 2: neutral, adaptive (any method), high confidence,
# environmental (RDA). All steps run in R.
#
# Missing genotypes: HO, HE, FIS and FST are calculated from genotypes WITHOUT
# imputation (hierfstat handles missing data itself), passed to hierfstat in its
# two-digit format (11, 12, 22). The sNMF analysis uses a matrix in which missing
# genotypes are replaced by the most common genotype at each locus.
# ======================================================================

# ── Setup: objects created in stage 2 (02_seascape_outliers.R) ─────────────
# Provides popmap, coords, G_imp, locus_names_full, pop_factor, gind_full, marine_dist,
# the outlier sets (outliers_neutral, outliers_any, outliers_highconf, outliers_rda,
# outliers_enviro) and the population colours (pop_cols_RA)
load("RA_seascape_FINAL_checkpoint.RData")

library(vcfR)
library(hierfstat)
library(parallel)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

pop_levels <- c("MZ","TZ","MD","RI","SY","MA","AU")   # mainland, islands, Indo-Pacific
pop_cols   <- pop_cols_RA

dataset_list <- list(
  neutral  = outliers_neutral,
  adaptive = outliers_any,
  highconf = outliers_highconf,
  enviro   = outliers_rda
)


# ======================================================================
# 1. Genotype matrices and diversity statistics (HO, HE, FIS)
# ======================================================================


# reload genotype matrix

obj.vcfR_s1 <- read.vcfR("RApop_set1.snps.vcf")
geno_s1     <- extract.gt(obj.vcfR_s1)
position_s1   <- getPOS(obj.vcfR_s1)
chromosome_s1 <- getCHROM(obj.vcfR_s1)
G_s1 <- matrix(NA, nrow=nrow(geno_s1), ncol=ncol(geno_s1))
G_s1[geno_s1 %in% c("0/0","0|0")] <- 0
G_s1[geno_s1 %in% c("0/1","1/0","1|0","0|1")] <- 1
G_s1[geno_s1 %in% c("1/1","1|1")] <- 2
G_imp_s1 <- apply(G_s1, 1, function(x){
  x[is.na(x)] <- as.numeric(names(sort(table(x),decreasing=TRUE)[1])); x})
if(nrow(G_imp_s1) != nrow(coords)) G_imp_s1 <- t(G_imp_s1)
locus_names_s1 <- paste(chromosome_s1, position_s1, sep="_")
colnames(G_imp_s1) <- locus_names_s1
rownames(G_imp_s1) <- colnames(geno_s1)
cat("Set1 reloaded:", dim(G_imp_s1), "\n")

# Genotypes without imputation (individuals x loci), used for HO, HE, FIS and FST
G_na <- t(G_s1)
colnames(G_na) <- locus_names_s1
rownames(G_na) <- colnames(geno_s1)
cat("Missing genotypes:", round(100 * mean(is.na(G_na)), 2), "%\n")

# hierfstat expects diploid genotypes as two-digit allele codes: 0/1/2 -> 11/12/22
to_hf <- function(G) {
  stopifnot(all(G %in% c(0, 1, 2, NA)))
  out <- G
  out[G == 0] <- 11; out[G == 1] <- 12; out[G == 2] <- 22
  out
}

# Check: neutral global FST without imputation (0.3222219 for this dataset)
G_check <- G_na[, colnames(G_na) %in% outliers_neutral]
cat("Neutral global FST:",
    round(wc(cbind(pop = as.numeric(factor(pop_factor)), as.data.frame(to_hf(G_check))))$FST, 7), "\n")

## Population diversity: HO, HE, FIS for the four datasets, MD and AU excluded (n = 2 each)
# The FIS bootstrap resamples loci: the per-locus FIS values from basic.stats are
# resampled, which gives the same result as recomputing basic.stats on resampled loci.

exclude_pops <- c("MD","AU")
keep_idx     <- which(!pop_factor %in% exclude_pops)
pop_filt     <- factor(pop_factor[keep_idx])
pop_num      <- as.numeric(pop_filt)
pop_names    <- levels(pop_filt)
cat("Populations:", pop_names, "\n")
cat("N per pop:", table(pop_filt), "\n")

nboot       <- 500
div_results <- list()

for(ds_name in names(dataset_list)){
  cat("\n===", ds_name, "===\n")
  loci  <- dataset_list[[ds_name]]
  G_div <- G_na[keep_idx, colnames(G_na) %in% loci]
  cat("SNPs:", ncol(G_div), "| Individuals:", nrow(G_div), "\n")

  hf <- cbind(pop = pop_num, as.data.frame(to_hf(G_div)))
  bs <- basic.stats(hf)

  ho_pop  <- colMeans(bs$Ho,  na.rm = TRUE)
  he_pop  <- colMeans(bs$Hs,  na.rm = TRUE)
  fis_pop <- colMeans(bs$Fis, na.rm = TRUE)

  set.seed(42)
  fis_boot <- t(replicate(nboot, colMeans(
    bs$Fis[sample(nrow(bs$Fis), replace = TRUE), , drop = FALSE], na.rm = TRUE)))
  fis_ll <- apply(fis_boot, 2, quantile, 0.025, na.rm = TRUE)
  fis_ul <- apply(fis_boot, 2, quantile, 0.975, na.rm = TRUE)

  div_table <- data.frame(
    Population = pop_names,
    n          = as.numeric(table(pop_filt)[pop_names]),
    HO         = round(ho_pop, 4),
    HE         = round(he_pop, 4),
    FIS        = round(fis_pop, 4),
    FIS_lower  = round(fis_ll, 4),
    FIS_upper  = round(fis_ul, 4),
    row.names  = NULL
  )
  print(div_table)
  write.csv(div_table, paste0("RA_set1_", ds_name, "_diversity_table.csv"), row.names = FALSE)
  div_results[[ds_name]] <- div_table
}

save(div_results, file = "RA_set1_diversity_results.RData")
cat("\nDone\n")


# ======================================================================
# 2. Private alleles (poppr) and allelic richness (hierfstat)
# ======================================================================

# Madagascar (n = 2) and Australia (n = 2) are excluded from population-level diversity
pops_for_div <- c("MZ","TZ","RI","SY","MA")
gind_div <- gind_full[which(as.character(pop(gind_full)) %in% pops_for_div), ]
pop(gind_div) <- factor(as.character(pop(gind_div)), levels = pops_for_div)

# ============================================================
# Private Alleles (poppr) + Allelic Richness (hierfstat)
# RA set1 — four datasets
# ============================================================

library(poppr)
library(hierfstat)

dataset_list <- list(
  neutral  = outliers_neutral,
  adaptive = outliers_any,
  highconf = outliers_highconf,
  enviro   = outliers_rda
)

ar_pa_results <- list()

for (ds_name in names(dataset_list)) {
  cat("\n===", ds_name, "===\n")
  
  loci_present <- dataset_list[[ds_name]][dataset_list[[ds_name]] %in% locNames(gind_div)]
  gind_sub     <- gind_div[, loc = loci_present]
  pop_order    <- levels(pop(gind_sub))
  
  # Allelic richness — rarefied to smallest pop (n=8, TZ)
  hf     <- genind2hierfstat(gind_sub)
  min_n  <- min(table(pop(gind_sub)))
  ar     <- allelic.richness(hf, min.n = min_n)
  ar_pop <- colMeans(ar$Ar, na.rm = TRUE)
  cat("  Rarefaction n =", min_n, "\n")
  
  # Private alleles — number of loci (not allele copies)
  pa_mat <- private_alleles(gind_sub, count.alleles = FALSE, drop = TRUE)
  pa_pop <- rowSums(pa_mat)
  
  # Summary table
  ar_pa_table <- data.frame(
    Dataset    = ds_name,
    Population = pop_order,
    n          = as.integer(table(pop(gind_sub))[pop_order]),
    AR         = round(ar_pop[pop_order], 4),
    PA         = as.integer(pa_pop[pop_order]),
    row.names  = NULL
  )
  
  print(ar_pa_table)
  ar_pa_results[[ds_name]] <- ar_pa_table
}

ar_pa_all <- do.call(rbind, ar_pa_results)
write.csv(ar_pa_all, "RA_set1_AR_PA_all_datasets.csv", row.names = FALSE)
save(ar_pa_results, file = "RA_set1_AR_PA_results.RData")
cat("\nDone.\n")


# ======================================================================
# 2b. Weighted means of HO, HE, FIS and AR across the five populations (weights = sample size)
# sum(n_i * statistic_i) / sum(n_i), with n = 67 (MD and AU excluded)
# Private alleles are not averaged; per-population counts are in RA_set1_AR_PA_all_datasets.csv
# ======================================================================

wm <- do.call(rbind, lapply(names(dataset_list), function(ds) {
  d1 <- div_results[[ds]]
  d2 <- ar_pa_results[[ds]]
  m  <- merge(d1[, c("Population", "n", "HO", "HE", "FIS")],
              d2[, c("Population", "AR")], by = "Population")
  data.frame(Dataset = ds,
             N       = sum(m$n),
             HO      = weighted.mean(m$HO,  m$n),
             HE      = weighted.mean(m$HE,  m$n),
             FIS     = weighted.mean(m$FIS, m$n),
             AR      = weighted.mean(m$AR,  m$n))
}))
wm[, -(1:2)] <- round(wm[, -(1:2)], 4)
print(wm)
write.csv(wm, "RA_set1_diversity_weighted_means.csv", row.names = FALSE)



# ======================================================================
# 3. Pairwise and global FST (hierfstat) with permutation tests
#    9 999 permutations of individuals between the two populations (all possible
#    arrangements when fewer exist), Benjamini-Hochberg adjustment within each dataset
# ======================================================================

nperm     <- 9999
n_cores   <- max(1, detectCores() - 1)
seed_base <- 123

# Weir and Cockerham (1984) theta for two populations, locus by locus, from the
# number of typed individuals (n), the allele dosage sums (s) and the numbers of
# heterozygotes (h) in each population; ratio of sums over loci, as in hierfstat.
# Missing data are handled as in hierfstat: each variance component (among
# populations a, among individuals b, within individuals c) is summed over the
# loci where it is defined. A locus typed in only one population contributes b
# and c but not a, and a locus typed in a single individual per population
# contributes c only. (Validated against hierfstat to about 1e-16.)
wc_theta <- function(n1, n2, s1, s2, h1, h2) {
  N    <- n1 + n2
  nbar <- N / 2
  nc   <- N - (n1^2 + n2^2) / N
  p1   <- s1 / (2 * n1);  p2 <- s2 / (2 * n2)
  pbar <- (n1 * p1 + n2 * p2) / N
  s2v  <- (n1 * (p1 - pbar)^2 + n2 * (p2 - pbar)^2) / nbar
  hbar <- (h1 + h2) / N
  a  <- (nbar / nc) * (s2v - (1 / (nbar - 1)) * (pbar * (1 - pbar) - s2v / 2 - hbar / 4))
  b  <- (nbar / (nbar - 1)) * (pbar * (1 - pbar) - s2v / 2 - (2 * nbar - 1) * hbar / (4 * nbar))
  cc <- hbar / 2
  one <- (n1 == 0 | n2 == 0) & N >= 1          # typed in one population only
  if (any(one)) {
    nn <- N[one]; pp <- (s1[one] + s2[one]) / (2 * nn); hh <- (h1[one] + h2[one]) / nn
    a[one]  <- NA
    b[one]  <- (nn / (nn - 1)) * (pp * (1 - pp) - (2 * nn - 1) * hh / (4 * nn))
    cc[one] <- hh / 2
  }
  a[!is.finite(a)] <- NA; b[!is.finite(b)] <- NA; cc[!is.finite(cc)] <- NA
  sum(a, na.rm = TRUE) / (sum(a, na.rm = TRUE) + sum(b, na.rm = TRUE) + sum(cc, na.rm = TRUE))
}

# Permutation test for one pair. Gp: individuals of the two populations x loci
# (0/1/2, NA = missing); the first n1 rows are population 1. Individuals are
# shuffled between the two populations (missing genotypes stay with their
# individual). If there are no more than nperm distinct arrangements, all of
# them are enumerated (exact p); otherwise nperm random shuffles are used with
# p = (1 + number >= observed) / (nperm + 1).
perm_pair <- function(Gp, n1, nperm, seed) {
  N  <- nrow(Gp)
  M  <- (!is.na(Gp)) * 1
  G0 <- Gp; G0[is.na(G0)] <- 0
  H0 <- (G0 == 1) * 1
  totM <- colSums(M); totG <- colSums(G0); totH <- colSums(H0)
  stat <- function(idx) {
    ind1 <- numeric(N); ind1[idx] <- 1
    m1 <- drop(ind1 %*% M); g1 <- drop(ind1 %*% G0); k1 <- drop(ind1 %*% H0)
    wc_theta(m1, totM - m1, g1, totG - g1, k1, totH - k1)
  }
  obs <- stat(seq_len(n1))
  if (choose(N, n1) <= nperm) {
    th <- apply(combn(N, n1), 2, stat)
    p      <- mean(th >= obs - 1e-12)
    p_min  <- mean(abs(th - max(th)) < 1e-12)
    method <- "exact"; B <- length(th)
  } else {
    set.seed(seed)
    th <- replicate(nperm, stat(sample.int(N, n1)))
    p      <- (1 + sum(th >= obs - 1e-12)) / (nperm + 1)
    p_min  <- 1 / (nperm + 1)
    method <- "permutation"; B <- nperm
  }
  list(obs = obs, p = p, p_min = p_min, method = method, B = B)
}

pop_full       <- factor(pop_factor)
pop_names_full <- levels(pop_full)
pairs          <- t(combn(pop_names_full, 2))

fst_results <- list()
fst_tables  <- list()

for (ds_name in names(dataset_list)) {
  cat("\n===", ds_name, "===\n")
  G_fst <- G_na[, colnames(G_na) %in% dataset_list[[ds_name]]]
  cat("SNPs:", ncol(G_fst), "| Individuals:", nrow(G_fst), "\n")

  # observed global and pairwise FST (hierfstat)
  hf_full      <- cbind(pop = as.numeric(pop_full), as.data.frame(to_hf(G_fst)))
  global_fst   <- wc(hf_full)$FST
  pairwise_fst <- pairwise.WCfst(hf_full)
  rownames(pairwise_fst) <- colnames(pairwise_fst) <- pop_names_full
  cat("Global FST:", round(global_fst, 4), "\n")

  # permutation tests for all pairs
  run_pair <- function(k) {
    idx_i <- which(pop_full == pairs[k, 1]); idx_j <- which(pop_full == pairs[k, 2])
    perm_pair(G_fst[c(idx_i, idx_j), , drop = FALSE], length(idx_i), nperm, seed_base + k)
  }
  cat("Permutation tests (", nrow(pairs), "pairs ) ...\n")
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("G_fst", "pop_full", "pairs", "nperm", "seed_base",
                        "perm_pair", "wc_theta"), envir = environment())
    res <- parLapply(cl, seq_len(nrow(pairs)), run_pair)
    stopCluster(cl)
  } else {
    res <- lapply(seq_len(nrow(pairs)), run_pair)
  }

  tab <- data.frame(
    Pop1 = pairs[, 1], Pop2 = pairs[, 2],
    n1 = as.numeric(table(pop_full)[pairs[, 1]]), n2 = as.numeric(table(pop_full)[pairs[, 2]]),
    FST_raw = sapply(res, `[[`, "obs"),
    p_perm  = sapply(res, `[[`, "p"),
    p_min   = sapply(res, `[[`, "p_min"),
    method  = sapply(res, `[[`, "method"),
    n_arrangements = sapply(res, `[[`, "B"))
  tab$p_BH     <- p.adjust(tab$p_perm, method = "BH")
  tab$FST      <- pmax(tab$FST_raw, 0)
  tab$testable <- tab$p_min < 0.05      # FALSE when the smallest attainable p is 0.05 or more

  # validation: own estimator must match hierfstat for every pair
  hier <- mapply(function(a, b) pairwise_fst[a, b], tab$Pop1, tab$Pop2)
  max_diff <- max(abs(hier - tab$FST_raw))
  cat("Largest difference between own and hierfstat pairwise FST:", signif(max_diff, 3), "\n")
  stopifnot(max_diff < 1e-6)

  print(tab[, c("Pop1","Pop2","n1","n2","FST","p_perm","p_BH","testable")])
  write.csv(tab, paste0("RA_set1_", ds_name, "_FST_permutation_table.csv"), row.names = FALSE)

  # matrices used by the heatmap and IBD code (p-values are BH-adjusted)
  pmat <- matrix(NA, length(pop_names_full), length(pop_names_full),
                 dimnames = list(pop_names_full, pop_names_full))
  for (k in seq_len(nrow(tab))) {
    pmat[tab$Pop1[k], tab$Pop2[k]] <- tab$p_BH[k]
    pmat[tab$Pop2[k], tab$Pop1[k]] <- tab$p_BH[k]
  }
  fst_results[[ds_name]] <- list(fst = pairwise_fst, pvals = pmat, global = global_fst)
  fst_tables[[ds_name]]  <- cbind(Dataset = ds_name, tab)
}

save(fst_results, fst_tables, file = "RA_set1_FST_results.RData")
write.csv(do.call(rbind, fst_tables), "RA_set1_FST_permutation_all_datasets.csv", row.names = FALSE)



# ======================================================================
# 4. FST heatmaps
# ======================================================================

pop_levels_RA <- c("MZ","TZ","MD","RI","SY","MA","AU")
pop_display_RA <- c(
  "MZ" = "Mozambique",
  "TZ" = "Tanzania",
  "MD" = "Madagascar",
  "RI" = "Réunion Island",
  "SY" = "Seychelles",
  "MA" = "Malaysia",
  "AU" = "Australia"
)



# Significance: * = Benjamini-Hochberg adjusted permutation p-value below sig_alpha.
# Negative FST values are shown as 0. Cells without an asterisk are not significant.

sig_alpha <- 0.05

make_fst_heatmap_RA <- function(fst_upper, fst_lower,
                                pval_upper, pval_lower,
                                filename_pdf, filename_png) {
  n <- length(pop_levels_RA)

  fst_upper  <- pmax(fst_upper[pop_levels_RA,  pop_levels_RA], 0)
  fst_lower  <- pmax(fst_lower[pop_levels_RA,  pop_levels_RA], 0)
  pval_upper <- pval_upper[pop_levels_RA, pop_levels_RA]
  pval_lower <- pval_lower[pop_levels_RA, pop_levels_RA]

  fst_comb  <- matrix(NA, n, n, dimnames = list(pop_levels_RA, pop_levels_RA))
  pval_comb <- matrix(NA, n, n, dimnames = list(pop_levels_RA, pop_levels_RA))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) { fst_comb[i,j] <- fst_upper[i,j]; pval_comb[i,j] <- pval_upper[i,j] }
      if (i > j) { fst_comb[i,j] <- fst_lower[i,j]; pval_comb[i,j] <- pval_lower[i,j] }
    }
  }

  fst_m  <- melt(fst_comb,  varnames = c("Pop1","Pop2"), value.name = "FST")
  pval_m <- melt(pval_comb, varnames = c("Pop1","Pop2"), value.name = "pval")
  df     <- merge(fst_m, pval_m, by = c("Pop1","Pop2"))

  df$sig <- ifelse(!is.na(df$pval) & df$pval < sig_alpha, "*", "")
  df$label <- ifelse(
    df$Pop1 == df$Pop2, "",
    ifelse(!is.na(df$FST),
           ifelse(df$sig != "",
                  paste0(round(df$FST, 3), "\n", df$sig),
                  as.character(round(df$FST, 3))),
           ""))

  df$Pop1 <- factor(df$Pop1, levels = pop_levels_RA)
  df$Pop2 <- factor(df$Pop2, levels = rev(pop_levels_RA))

  p <- ggplot(df, aes(x = Pop1, y = Pop2, fill = FST)) +
    geom_tile(color = NA) +
    geom_tile(data = subset(df, as.character(Pop1) == as.character(Pop2)),
              fill = "white", color = NA) +
    geom_text(aes(label = label), size = 2.8, color = "black") +
    scale_fill_gradientn(
      colors   = c("#FFFFFF","#C6DBEF","#6BAED6","#2171B5","#08306B"),
      na.value = "white",
      name     = expression(F[ST]),
      limits   = c(0, max(fst_comb, na.rm = TRUE))) +
    scale_x_discrete(labels = pop_display_RA) +
    scale_y_discrete(labels = pop_display_RA) +
    theme_bw(base_size = 11) +
    theme(axis.title      = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1),
          panel.grid      = element_blank(),
          panel.border    = element_blank(),
          axis.ticks      = element_blank(),
          legend.position = "right")

  pdf(filename_pdf, width = 6, height = 5.5); print(p); dev.off()
  cat("Saved:", filename_pdf, "\n")
  png(filename_png, width = 6, height = 5.5, units = "in", res = 300); print(p); dev.off()
  cat("Saved:", filename_png, "\n")
  return(p)
}

# Main: neutral in the matrix upper triangle (drawn below the diagonal because the y axis
# is reversed), adaptive in the matrix lower triangle (drawn above the diagonal)
p_main <- make_fst_heatmap_RA(
  fst_upper  = as.matrix(fst_results$neutral$fst),
  fst_lower  = as.matrix(fst_results$adaptive$fst),
  pval_upper = as.matrix(fst_results$neutral$pvals),
  pval_lower = as.matrix(fst_results$adaptive$pvals),
  filename_pdf = "RA_set1_FST_heatmap_neutral_adaptive.pdf",
  filename_png = "RA_set1_FST_heatmap_neutral_adaptive.png"
)

# Supplementary: high confidence in the matrix upper triangle (drawn below the diagonal),
# environmental in the matrix lower triangle (drawn above the diagonal)
p_supp <- make_fst_heatmap_RA(
  fst_upper  = as.matrix(fst_results$highconf$fst),
  fst_lower  = as.matrix(fst_results$enviro$fst),
  pval_upper = as.matrix(fst_results$highconf$pvals),
  pval_lower = as.matrix(fst_results$enviro$pvals),
  filename_pdf = "RA_set1_FST_heatmap_highconf_enviro.pdf",
  filename_png = "RA_set1_FST_heatmap_highconf_enviro.png"
)

# Global FST of each dataset, for the captions
cat("\nGlobal FST:\n")
print(round(sapply(fst_results, function(x) x$global), 4))



# ======================================================================
# 5. DAPC
# ======================================================================

## DAPC — Population Structure Analysis
## Replicable pipeline for RA set1 and RD
## Requires: VCF file, popmap, outlier sets, population colours

library(vcfR)
library(adegenet)
library(ggplot2)
library(patchwork)
library(cowplot)

## ── 1. Load VCF and build genlight ──────────────────────────────────────────

vcf <- read.vcfR("RApop_set1.snps.vcf")
gl  <- vcfR2genlight(vcf)

# Fix locus names — vcfR2genlight uses "CHR:POS" (0-based)
# Convert to "CHR_POS" (1-based) to match matrix locus names
loc_fixed <- gsub(":", "_", locNames(gl))
loc_split <- strsplit(loc_fixed, "_")
loc_fixed <- sapply(loc_split, function(x){
  paste0(x[1], "_", as.integer(x[2])+1)
})
locNames(gl) <- loc_fixed

# Assign population factor
pop_factor <- popmap$pop[match(indNames(gl), popmap$ID)]
pop_full   <- factor(pop_factor, levels=pop_levels)
pop(gl)    <- pop_full

cat("Genlight object:", nInd(gl), "individuals,", nLoc(gl), "loci\n")
cat("Populations:", levels(pop_full), "\n")

# Population colours — must match pop_levels order
pop_cols_dapc <- pop_cols[levels(pop_full)]

## ── 2. Define datasets ──────────────────────────────────────────────────────

# outlier sets were loaded in the setup block

dataset_list <- list(
  neutral  = outliers_neutral,
  adaptive = outliers_any,
  highconf = outliers_highconf,
  enviro   = outliers_rda
)

## ── 3. DAPC plotting function ────────────────────────────────────────────────

make_dapc_plot <- function(dapc_obj, filename){

  scores_df     <- as.data.frame(dapc_obj$ind.coord)
  scores_df$pop <- dapc_obj$grp

  eig_df <- data.frame(
    axis = 1:length(dapc_obj$eig),
    eig  = dapc_obj$eig
  )

  xlab <- paste0("DA1 (", round(dapc_obj$eig[1]/sum(dapc_obj$eig)*100, 1), "%)")
  ylab <- paste0("DA2 (", round(dapc_obj$eig[2]/sum(dapc_obj$eig)*100, 1), "%)")

  pop_labels <- c("TZ"="Tanzania",      "MD"="Madagascar",
                  "SY"="Seychelles",    "AU"="Australia",
                  "RI"="Réunion Island","MZ"="Mozambique",
                  "MA"="Malaysia")

  p_main <- ggplot(scores_df, aes(x=LD1, y=LD2, fill=pop)) +
    geom_point(shape=21, size=3.5, stroke=0.4, color="grey30") +
    scale_fill_manual(values=pop_cols_dapc,
                      labels=pop_labels[levels(scores_df$pop)],
                      name="Population") +
    geom_hline(yintercept=0, linetype="dashed",
               color="grey70", linewidth=0.4) +
    geom_vline(xintercept=0, linetype="dashed",
               color="grey70", linewidth=0.4) +
    xlab(xlab) + ylab(ylab) +
    theme_bw(base_size=11) +
    theme(panel.grid=element_blank(),
          legend.position="right",
          legend.title=element_text(size=9, face="bold"),
          legend.text=element_text(size=8),
          legend.key.size=unit(0.4, "cm"))

  p_eig <- ggplot(eig_df, aes(x=factor(axis), y=eig)) +
    geom_bar(stat="identity", fill="grey40", width=0.6) +
    ggtitle("DA eigenvalues") +
    theme_bw(base_size=8) +
    theme(panel.grid=element_blank(),
          axis.title=element_blank(),
          axis.text=element_text(size=7),
          plot.title=element_text(size=8, hjust=0.5),
          plot.margin=margin(2,2,2,2))

  p_main_noleg <- p_main + theme(legend.position="none")
  leg          <- cowplot::get_legend(p_main)
  right_col    <- plot_grid(leg, p_eig, ncol=1, rel_heights=c(0.6, 0.4))
  p_final      <- plot_grid(p_main_noleg, right_col,
                             ncol=2, rel_widths=c(0.78, 0.22))

  pdf(filename, width=8, height=6)
    print(p_final)
  dev.off()

  return(p_final)
}

## ── 4. DAPC loop — all four datasets ────────────────────────────────────────

dapc_summary <- data.frame()

for(ds_name in names(dataset_list)){
  cat("\n===", ds_name, "===\n")
  loci   <- dataset_list[[ds_name]]
  gl_sub <- gl[, locNames(gl) %in% loci]
  cat("SNPs:", nLoc(gl_sub), "| Individuals:", nInd(gl_sub), "\n")

  # Matrix for xvalDapc
  geno_mat <- tab(gl_sub, NA.method="mean")

  # PCA axes explaining 80% variance
  # NOTE: returnDotProduct argument removed — not supported in all versions
  pca_sub <- glPca(gl_sub, nf=min(nInd(gl_sub)-1, nLoc(gl_sub)))
  var_exp <- cumsum(pca_sub$eig / sum(pca_sub$eig))
  n_pc_80 <- which(var_exp >= 0.80)[1]
  cat("PCA axes (80% variance):", n_pc_80, "\n")

  # find.clusters — 10 replicates, K 1 to 10
  cat("Running find.clusters (10 replicates)...\n")
  best_k_vec <- integer(10)
  for(rep in 1:10){
    fc <- find.clusters(gl_sub,
                        max.n.clust    = 10,
                        n.pca          = n_pc_80,
                        choose.n.clust = FALSE,
                        criterion      = "diffNgroup",
                        n.iter         = 1e6,
                        n.start        = 10)
    best_k_vec[rep] <- length(levels(fc$grp))
    cat("  Rep", rep, "— K =", best_k_vec[rep], "\n")
  }
  kmeans_k <- as.integer(names(sort(table(best_k_vec), decreasing=TRUE)[1]))
  stable_k <- kmeans_k
  cat("K-means optimal K:", kmeans_k, "\n")
  cat("Stable K:",          stable_k,  "\n")

  # xvalDapc — cross-validation, 100 reps, 90% training set
  # NOTE: correct argument is n.pca not n.pca.range
  cat("Running xvalDapc (100 reps)...\n")
  set.seed(42)
  xval <- xvalDapc(geno_mat,
                   pop(gl_sub),
                   n.pca        = 1:min(n_pc_80*2, 50),
                   training.set = 0.9,
                   result       = "groupMean",
                   center       = TRUE,
                   scale        = FALSE,
                   n.rep        = 100,
                   xval.plot    = FALSE)
  n_pc_xval <- as.integer(xval$`Number of PCs Achieving Lowest MSE`)
  cat("Cross-val optimal PCs:", n_pc_xval, "\n")

  # optim.a.score — 100 simulations
  cat("Running optim.a.score (100 simulations)...\n")
  dapc_tmp <- dapc(gl_sub, pop(gl_sub),
                   n.pca = n_pc_80,
                   n.da  = nlevels(pop(gl_sub))-1)
  ascore      <- optim.a.score(dapc_tmp, n.sim=100, plot=FALSE)
  n_pc_ascore <- ascore$best
  cat("A-score optimal PCs:", n_pc_ascore, "\n")

  # DAPC with priors — cross-val PC retention (main text)
  dapc_prior <- dapc(gl_sub, pop(gl_sub),
                     n.pca = n_pc_xval,
                     n.da  = nlevels(pop(gl_sub))-1)

  # DAPC without priors — stable K from find.clusters (supplementary)
  fc_final     <- find.clusters(gl_sub,
                                 n.clust        = stable_k,
                                 n.pca          = n_pc_80,
                                 choose.n.clust = FALSE,
                                 n.iter         = 1e6,
                                 n.start        = 10)
  dapc_noprior <- dapc(gl_sub, fc_final$grp,
                        n.pca = n_pc_xval,
                        n.da  = stable_k-1)

  # BIC plot
  pdf(paste0("RA_set1_", ds_name, "_BIC.pdf"), width=6, height=5)
    plot(fc$Kstat, type="b", pch=19,
         xlab="Number of clusters (K)", ylab="BIC",
         cex.lab=1.1, las=1, bty="l")
    abline(v=stable_k, lty=2, col="red")
  dev.off()

  # DAPC plots
  # With priors — main text
  p_prior <- make_dapc_plot(dapc_prior,
                             paste0("RA_set1_", ds_name, "_dapc_prior.pdf"))

  # Without priors — supplementary
  pdf(paste0("RA_set1_", ds_name, "_dapc_noprior.pdf"), width=8, height=6)
    scatter(dapc_noprior, bg="white", solid=0.7, cex=1.5,
            clab=0, leg=TRUE, posi.leg="topright")
    title(main=paste0(ds_name, " — without priors (K=", stable_k, ")"),
          cex.main=0.9)
  dev.off()

  # Save DAPC objects
  save(dapc_prior, dapc_noprior, fc_final, xval, ascore,
       file=paste0("RA_set1_", ds_name, "_dapc_results.RData"))

  # Save plot object for editing
  saveRDS(p_prior, paste0("RA_set1_", ds_name, "_dapc_plot.rds"))

  # Summary table row
  dapc_summary <- rbind(dapc_summary, data.frame(
    Dataset          = ds_name,
    N                = nLoc(gl_sub),
    PCA_axes_80pct   = n_pc_80,
    Kmeans_optimal_K = kmeans_k,
    Stable_K         = stable_k,
    Ascore_PCs       = n_pc_ascore,
    Crossval_PCs     = n_pc_xval
  ))

  cat("===", ds_name, "complete ===\n")
}

## ── 5. Summary table

print(dapc_summary)
write.csv(dapc_summary, "RA_set1_DAPC_summary_table.csv", row.names=FALSE)
cat("\nAll DAPC analyses complete\n")

# ======================================================================
# 6. Adaptive DAPC excluding Australia (n = 2)
# ======================================================================

# genlight object for all 71 individuals (called gl_RA because gl is also a base R function)
gl_RA <- vcfR2genlight(read.vcfR("RApop_set1.snps.vcf", verbose = FALSE))
pop(gl_RA) <- factor(pop_factor, levels = c("MZ","TZ","MD","RI","SY","MA","AU"))
locNames(gl_RA) <- gsub(":", "_", locNames(gl_RA))
stopifnot(nInd(gl_RA) == 71, nLoc(gl_RA) == 27562)

# adaptive dapc without AU
# vcfR2genlight uses 0-based POS — shift by +1 to match locus_names_full
# locus_names_full uses 1-based POS from getPOS(obj.vcfR)
# Fix: rebuild locus names from gl_RA using +1 offset

# Option A — use locus_names_full directly (confirmed correct)
stopifnot(length(locus_names_full) == nLoc(gl_RA))
locNames(gl_RA) <- locus_names_full
cat("Locus names corrected from locus_names_full\n")
cat("First 3:", locNames(gl_RA)[1:3], "\n")
# Expected: 1_20  2_7  5_243

# Verify subset works
test_sub <- gl_RA[, locNames(gl_RA) %in% outliers_any]
cat("Adaptive loci found:", nLoc(test_sub), "\n")
# Expected: 3147
stopifnot(nLoc(test_sub) == length(outliers_any))
cat("Locus name fix verified ✓\n")

# ============================================================
# Adaptive DAPC — excluding Australia (n=2)
# Australia excluded: extreme discriminant scores (DA2≈-200)
# compress WIO + Malaysia structure in full dataset plot
# ============================================================

library(adegenet)

# ── Subset to adaptive loci excluding AU ─────────────────────────────────────
idx_noAU     <- which(as.character(pop_factor) != "AU")
gl_adap_noAU <- gl_RA[idx_noAU, ]
gl_adap_noAU <- gl_adap_noAU[,
  locNames(gl_adap_noAU) %in% outliers_any]
pop_noAU     <- droplevels(
  factor(as.character(pop_factor[idx_noAU]),
         levels = c("MZ","TZ","MD","RI","SY","MA")))
pop(gl_adap_noAU) <- pop_noAU

cat("No-AU adaptive genlight:\n")
cat("Individuals:", nInd(gl_adap_noAU), "\n")
cat("Loci:", nLoc(gl_adap_noAU), "\n")
cat("Populations:", table(pop_noAU), "\n")
# Expected: 69 individuals, 3147 loci

# ── Step 1: PCA axes explaining 80% variance ─────────────────────────────────
pca_adap <- glPca(gl_adap_noAU,
                   nf = min(nInd(gl_adap_noAU) - 1,
                             nLoc(gl_adap_noAU)))
var_exp  <- cumsum(pca_adap$eig / sum(pca_adap$eig))
n_pc_80  <- which(var_exp >= 0.80)[1]
cat("PCs explaining 80% variance:", n_pc_80, "\n")

# ── Step 2: find.clusters — 10 replicates ────────────────────────────────────
cat("Running find.clusters (10 replicates)...\n")
best_k_vec_noAU <- integer(10)
for (rep in 1:10) {
  fc <- find.clusters(gl_adap_noAU,
                       max.n.clust    = 8,
                       n.pca          = n_pc_80,
                       choose.n.clust = FALSE,
                       criterion      = "diffNgroup",
                       n.iter         = 1e6,
                       n.start        = 10)
  best_k_vec_noAU[rep] <- length(levels(fc$grp))
  cat("  Rep", rep, "— K =", best_k_vec_noAU[rep], "\n")
}
k_noAU <- as.integer(
  names(sort(table(best_k_vec_noAU),
             decreasing = TRUE)[1]))
cat("Optimal K (no AU):", k_noAU, "\n")

# ── Step 3: xvalDapc ──────────────────────────────────────────────────────────
cat("Running xvalDapc (100 reps)...\n")
geno_mat_noAU <- tab(gl_adap_noAU, NA.method = "mean")
set.seed(123)
xval_noAU <- xvalDapc(geno_mat_noAU,
                       grp        = pop_noAU,
                       n.pca.max  = n_pc_80,
                       training.set = 0.9,
                       result     = "groupMean",
                       center     = TRUE,
                       scale      = FALSE,
                       n.rep      = 100,
                       xval.plot  = FALSE)
n_pc_xval_noAU <- as.integer(
  xval_noAU$`Number of PCs Achieving Lowest MSE`)
cat("PCs retained (xvalDapc):", n_pc_xval_noAU, "\n")

# ── Step 4: Final DAPC ────────────────────────────────────────────────────────
set.seed(123)
dapc_adap_noAU <- dapc(gl_adap_noAU,
                        pop        = pop_noAU,
                        n.pca      = n_pc_xval_noAU,
                        n.da       = k_noAU - 1)

cat("DAPC complete\n")
cat("DA1 variance:", round(
  dapc_adap_noAU$eig[1]/sum(dapc_adap_noAU$eig)*100, 1), "%\n")
cat("DA2 variance:", round(
  dapc_adap_noAU$eig[2]/sum(dapc_adap_noAU$eig)*100, 1), "%\n")

cat("Cluster assignments:\n")
print(table(dapc_adap_noAU$assign, pop_noAU))

# ── Step 5: Plot ──────────────────────────────────────────────────────────────
pop_cols_RA <- c(
  "MZ" = "#0072B2", "TZ" = "#E69F00", "MD" = "#56B4E9",
  "RI" = "#D55E00", "SY" = "#009E73", "MA" = "#FF69B4")

pt_cols_noAU <- pop_cols_RA[as.character(pop_noAU)]

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf(
    "RA_set1_adaptive_dapc_noAU.pdf",
    width = 7, height = 6)
  else png("RA_set1_adaptive_dapc_noAU.png",
           width = 7, height = 6,
           units = "in", res = 300)

  scatter(dapc_adap_noAU,
          col        = pop_cols_RA[levels(pop_noAU)],
          bg         = "white",
          scree.da   = TRUE,
          posi.da    = "bottomright",
          cex        = 1.5,
          solid      = 0.7,
          pch        = 19,
          leg        = TRUE,
          txt.leg    = levels(pop_noAU),
          posi.leg   = "topleft",
          clab       = 0,
          clabel     = 0)
  dev.off()
}
cat("Adaptive DAPC no-AU plots saved\n")

# ── Save ──────────────────────────────────────────────────────────────────────
save(dapc_adap_noAU, k_noAU, n_pc_xval_noAU,
     best_k_vec_noAU,
     file = "RA_set1_adaptive_dapc_noAU_results.RData")
cat("Results saved\n")


# ======================================================================
# 7. sNMF: ancestry coefficients
# ======================================================================

  ## sNMF Population Structure — RA set1, all four datasets
## LEA snmf, K=1-10, 10 replicates, alpha=100
## Visualisation via pophelper v2.3.1

library(LEA)
library(pophelper)
library(RColorBrewer)

## ── 1. Individual metadata ordered by population

ind_meta <- data.frame(
  ID  = rownames(G_imp_s1),
  pop = pop_factor,
  stringsAsFactors = FALSE
)
ind_meta$pop <- factor(ind_meta$pop, levels=pop_levels)
ind_meta     <- ind_meta[order(ind_meta$pop),]

grp_labels   <- data.frame(pop=as.character(ind_meta$pop),
                             stringsAsFactors=FALSE)

## ── 2. Dataset list

dataset_list <- list(
  neutral  = outliers_neutral,
  adaptive = outliers_any,
  highconf = outliers_highconf,
  enviro   = outliers_rda
)

## ── 3. Run sNMF — K=1-10, 10 replicates, alpha=100

for(ds_name in names(dataset_list)){
  cat("\n===", ds_name, "===\n")
  loci  <- dataset_list[[ds_name]]
  G_sub <- G_imp_s1[, colnames(G_imp_s1) %in% loci]
  cat("SNPs:", ncol(G_sub), "| Individuals:", nrow(G_sub), "\n")

  # Reorder individuals by population
  G_sub <- G_sub[ind_meta$ID,]

  # Write .geno file
  # Format: rows = SNPs, cols = individuals, values 0/1/2, NA -> 9
  G_t <- t(G_sub)
  G_t[is.na(G_t)] <- 9
  geno_file <- paste0("RA_set1_", ds_name, ".geno")
  write.table(G_t, file=geno_file,
              row.names=FALSE, col.names=FALSE,
              sep="", quote=FALSE)
  cat("Geno file written:", geno_file, "\n")

  # Remove existing project if present
  proj_file <- paste0("RA_set1_", ds_name, ".snmfProject")
  if(file.exists(proj_file)) remove.snmfProject(proj_file)

  # Run sNMF
  snmf_proj <- snmf(geno_file,
                    K       = 1:10,
                    rep     = 10,
                    alpha   = 100,
                    entropy = TRUE,
                    ploidy  = 2,
                    project = "new")

  # Cross-entropy — mean across replicates per K
  ce_vals  <- sapply(1:10, function(k) mean(cross.entropy(snmf_proj, K=k)))
  best_k   <- which.min(ce_vals)
  best_run <- which.min(cross.entropy(snmf_proj, K=best_k))
  cat("Optimal K:", best_k, "\n")
  cat("Best run:", best_run, "\n")

  # Cross-entropy plot
  pdf(paste0("RA_set1_", ds_name, "_snmf_crossentropy.pdf"), width=6, height=5)
    plot(1:10, ce_vals, type="b", pch=19,
         xlab="Number of ancestral populations (K)",
         ylab="Cross-entropy",
         cex.lab=1.2, cex.axis=1.1, las=1, bty="l")
    abline(v=best_k, lty=2, col="red")
    legend("topright", legend=paste0("Optimal K = ", best_k),
           bty="n", cex=0.9)
  dev.off()

  # Save Q matrix
  qmat <- as.data.frame(Q(snmf_proj, K=best_k, run=best_run))
  colnames(qmat)  <- paste0("Cluster", 1:best_k)
  rownames(qmat)  <- rownames(G_imp_s1)
  qmat_ordered    <- qmat[ind_meta$ID,]
  write.csv(qmat_ordered,
            paste0("RA_set1_", ds_name, "_snmf_Qmatrix_K", best_k, ".csv"),
            row.names=TRUE)

  cat("===", ds_name, "complete ===\n")
}

## ── 4. Barplots — all four datasets

dark2_cols <- brewer.pal(8, "Dark2")

snmf_info <- list(
  neutral  = list(proj_file="RA_set1_neutral.snmfProject",  best_k=3),
  adaptive = list(proj_file="RA_set1_adaptive.snmfProject", best_k=5),
  highconf = list(proj_file="RA_set1_highconf.snmfProject", best_k=3),
  enviro   = list(proj_file="RA_set1_enviro.snmfProject",   best_k=3)
)

snmf_summary <- data.frame()

for(ds_name in names(snmf_info)){
  cat("\n===", ds_name, "===\n")
  info   <- snmf_info[[ds_name]]
  best_k <- info$best_k

  proj     <- load.snmfProject(info$proj_file)
  best_run <- which.min(cross.entropy(proj, K=best_k))
  cat("Best run:", best_run, "\n")

  # Extract Q matrix with correct sample IDs
  qmat <- as.data.frame(Q(proj, K=best_k, run=best_run))
  colnames(qmat)  <- paste0("Cluster", 1:best_k)
  rownames(qmat)  <- rownames(G_imp_s1)
  qmat            <- qmat[ind_meta$ID,]
  cat("NAs:", sum(is.na(qmat)), "| Dimensions:", dim(qmat), "\n")

  # pophelper qlist
  qlist <- as.qlist(list(qmat))

  # Barplot
  plotQ(qlist,
        clustercol      = dark2_cols[1:best_k],
        grplab          = grp_labels,
        selgrp          = "pop",
        ordergrp        = TRUE,
        showgrplab      = TRUE,
        grplabsize      = 4,
        grplabangle     = 0,
        grplabjust      = 0.5,
        grplabcol       = "grey30",
        showindlab      = FALSE,
        sortind         = "all",
        barbordercolour = "black",
        barbordersize   = 0.1,
        divtype         = 1,
        divcol          = "white",
        divsize         = 0.5,
        panelratio      = c(8,1),
        grplabheight    = 0.8,
        splab           = paste0("K = ", best_k),
        splabsize       = 4,
        showsp          = TRUE,
        outputfilename  = paste0("RA_set1_", ds_name,
                                  "_snmf_barplot_K", best_k),
        imgtype         = "pdf",
        exportplot      = TRUE,
        exportpath      = getwd(),
        height          = 6,
        width           = 40,
        units           = "cm",
        dpi             = 300)

  # Summary table row
  snmf_summary <- rbind(snmf_summary, data.frame(
    Dataset  = ds_name,
    N_SNPs   = length(dataset_list[[ds_name]]),
    Optimal_K = best_k,
    Best_run  = best_run
  ))

  cat(ds_name, "barplot saved\n")
}

## ── 5. Summary table

print(snmf_summary)
write.csv(snmf_summary, "RA_set1_snmf_summary.csv", row.names=FALSE)
cat("\nAll sNMF analyses complete\n")


# ── Final barplots with colours fixed to geographic groups ──────────────────
# ============================================================
# RA Set1 — sNMF barplots with consistent colours
# Population order: MZ, TZ, MD, RI, SY, MA, AU
# Dark2 palette — colours fixed to geographic groups
# ============================================================

# ============================================================
# RA Set1 — sNMF barplots with consistent colours
# Population order: MZ, TZ, MD, RI, SY, MA, AU
# ============================================================

library(LEA)
library(pophelper)
library(RColorBrewer)

# ── Already loaded from checkpoint ───────────────────────────────────────────
# G_imp, pop_factor, locus_names_full, dark2_cols all present

cat("G_imp:", dim(G_imp), "\n")
cat("pop_factor:", table(pop_factor), "\n")

# ── Population order ──────────────────────────────────────────────────────────
pop_levels_RA <- c("MZ","TZ","MD","RI","SY","MA","AU")

pop_display_RA <- c(
  "MZ" = "Mozambique",
  "TZ" = "Tanzania",
  "MD" = "Madagascar",
  "RI" = "Réunion Island",
  "SY" = "Seychelles",
  "MA" = "Malaysia",
  "AU" = "Australia"
)

# ── Reorder individuals by new pop order ──────────────────────────────────────
ind_meta <- data.frame(
  ID  = rownames(G_imp),
  pop = pop_factor,
  stringsAsFactors = FALSE
)
ind_meta$pop <- factor(ind_meta$pop, levels = pop_levels_RA)
ind_meta     <- ind_meta[order(ind_meta$pop), ]

grp_labels <- data.frame(
  pop = as.character(ind_meta$pop),
  stringsAsFactors = FALSE)

cat("Individual order — population counts:\n")
print(table(ind_meta$pop))

# ── Load sNMF projects ────────────────────────────────────────────────────────
proj_neutral  <- load.snmfProject("RA_set1_neutral.snmfProject")
proj_adaptive <- load.snmfProject("RA_set1_adaptive.snmfProject")
proj_highconf <- load.snmfProject("RA_set1_highconf.snmfProject")
proj_enviro   <- load.snmfProject("RA_set1_enviro.snmfProject")

# ── Best runs ─────────────────────────────────────────────────────────────────
best_run_neutral  <- which.min(cross.entropy(proj_neutral,  K = 3))
best_run_adaptive <- which.min(cross.entropy(proj_adaptive, K = 5))
best_run_highconf <- which.min(cross.entropy(proj_highconf, K = 3))
best_run_enviro   <- which.min(cross.entropy(proj_enviro,   K = 3))

cat("Best runs — neutral:", best_run_neutral,
    "| adaptive:", best_run_adaptive,
    "| highconf:", best_run_highconf,
    "| enviro:", best_run_enviro, "\n")

# ── Extract Q matrices in new individual order ────────────────────────────────
get_qmat <- function(proj, K, run) {
  qmat <- as.data.frame(Q(proj, K = K, run = run))
  colnames(qmat) <- paste0("Cluster", 1:K)
  rownames(qmat) <- rownames(G_imp)
  qmat[ind_meta$ID, ]
}

qmat_neutral  <- get_qmat(proj_neutral,  3, best_run_neutral)
qmat_adaptive <- get_qmat(proj_adaptive, 5, best_run_adaptive)
qmat_highconf <- get_qmat(proj_highconf, 3, best_run_highconf)
qmat_enviro   <- get_qmat(proj_enviro,   3, best_run_enviro)

# ── Check cluster identity ────────────────────────────────────────────────────
check_clusters <- function(qmat, ds_name) {
  cat("\n===", ds_name, "— Population mean Q values ===\n")
  K <- ncol(qmat)
  for (pop in pop_levels_RA) {
    idx <- which(ind_meta$pop == pop)
    if (length(idx) == 0) next
    means <- round(colMeans(qmat[idx, , drop=FALSE]), 3)
    dominant <- which.max(means)
    cat(sprintf("%-4s → dominant Cluster%d | %s\n",
                pop, dominant,
                paste(paste0("C", 1:K, "=", means), collapse=" ")))
  }
}

check_clusters(qmat_neutral,  "Neutral")
check_clusters(qmat_adaptive, "Adaptive")
check_clusters(qmat_highconf, "Highconf")
check_clusters(qmat_enviro,   "Enviro")

# Biological colour mapping
# Indo-Pacific (MA + AU) = Dark2 colour 1 (teal)
# Islands (MD + RI + SY) = Dark2 colour 2 (orange)  
# Mainland Africa (MZ + TZ) = Dark2 colour 3 (purple)
# MZ-specific (adaptive only) = Dark2 colour 4 (pink)
# TZ-specific (adaptive only) = Dark2 colour 5 (green)

dark2_cols <- RColorBrewer::brewer.pal(5, "Dark2")

# Colour vectors per dataset — reorder Dark2 to match cluster numbers
# Neutral K=3: C1=IndoPacific, C2=Mainland, C3=Islands
cols_neutral <- c(dark2_cols[1],   # C1 = MA/AU = teal
                  dark2_cols[3],   # C2 = MZ/TZ = purple
                  dark2_cols[2])   # C3 = MD/RI/SY = orange

# Adaptive K=5: C1=Islands, C2=MA, C3=AU, C4=MZ, C5=TZ
cols_adaptive <- c(dark2_cols[2],  # C1 = MD/RI/SY = orange
                   dark2_cols[1],  # C2 = MA = teal
                   dark2_cols[1],  # C3 = AU = teal (same group as MA)
                   dark2_cols[3],  # C4 = MZ = purple
                   dark2_cols[4])  # C5 = TZ = pink

# Highconf K=3: C1=IndoPacific, C2=Islands, C3=Mainland
cols_highconf <- c(dark2_cols[1],  # C1 = MA/AU = teal
                   dark2_cols[2],  # C2 = MD/RI/SY = orange
                   dark2_cols[3])  # C3 = MZ/TZ = purple

# Enviro K=3: C1=MA/AU/TZ, C2=MZ, C3=Islands
# Note: in the environmental dataset TZ clusters with the Indo-Pacific group, not the mainland
cols_enviro <- c(dark2_cols[1],    # C1 = TZ/MA/AU = teal
                 dark2_cols[3],    # C2 = MZ = purple
                 dark2_cols[2])    # C3 = MD/RI/SY = orange

# ── sNMF barplot function ──────────────────────────────────────────────────────
plot_snmf <- function(qmat, cols, K, ds_name, filename) {

  qlist <- as.qlist(list(qmat))

  for (ext in c("pdf","png")) {
    outfile <- paste0(filename, ".", ext)
    plotQ(qlist,
          clustercol      = cols,
          grplab          = grp_labels,
          selgrp          = "pop",
          ordergrp        = FALSE,
          showgrplab      = TRUE,
          grplabsize      = 3.5,
          grplabangle     = 45,
          grplabjust      = 1,
          grplabcol       = "grey30",
          showindlab      = FALSE,
          sortind         = "all",
          barbordercolour = "black",
          barbordersize   = 0.1,
          divtype         = 1,
          divcol          = "white",
          divsize         = 0.5,
          panelratio      = c(8,1),
          grplabheight    = 0.8,
          splab           = paste0("K = ", K),
          splabsize       = 4,
          showsp          = TRUE,
          outputfilename  = paste0(filename, "_", ext),
          imgtype         = ext,
          exportplot      = TRUE,
          exportpath      = getwd(),
          height          = 6,
          width           = 50,
          units           = "cm",
          dpi             = 300)
    cat(ds_name, ext, "saved\n")
  }
}

#updated AU colour for adaptive
cols_adaptive <- c(dark2_cols[2],  # C1 = MD/RI/SY = orange
                   dark2_cols[1],  # C2 = MA = teal
                   dark2_cols[5],  # C3 = AU = green
                   dark2_cols[3],  # C4 = MZ = purple
                   dark2_cols[4])  # C5 = TZ = pink

cat("Adaptive colours:\n")
print(cols_adaptive)

# Rerun adaptive barplot only
plot_snmf(qmat_adaptive, cols_adaptive, 5, "Adaptive",
          "RA_set1_adaptive_snmf_K5")
cat("Adaptive sNMF barplot regenerated\n")

# ── Generate all four barplots ─────────────────────────────────────────────────
plot_snmf(qmat_neutral,  cols_neutral,  3, "Neutral",
          "RA_set1_neutral_snmf_K3")

plot_snmf(qmat_adaptive, cols_adaptive, 5, "Adaptive",
          "RA_set1_adaptive_snmf_K5")

plot_snmf(qmat_highconf, cols_highconf, 3, "Highconf",
          "RA_set1_highconf_snmf_K3")

plot_snmf(qmat_enviro,   cols_enviro,   3, "Enviro",
          "RA_set1_enviro_snmf_K3")

cat("All sNMF barplots complete\n")


# ======================================================================
# 8. IBD and IBE
# ======================================================================

## IBD/IBE Analysis — RA set1, all four datasets
## GLS with permutation testing (10,000 permutations)
## Geographic: marine least-cost distances averaged to population level
## Environmental: Bray-Curtis, Manhattan, Canberra — best selected by AIC
## Canberra selected for final plots (best for 3/4 datasets)
## Genetic: linearized FST (FST/1-FST)

library(nlme)
library(vegan)
library(ggplot2)

## ── 1. Population-level geographic distance matrix ───────────────────────────

diag(geo_dist_pop) <- 0

# Reload if necessary
# coords$pop <- popmap$pop[match(coords$ID, popmap$ID)]
# pop_coords <- aggregate(cbind(as.numeric(lat), as.numeric(lon)) ~ pop,
#                          data=coords, FUN=mean)
# colnames(pop_coords) <- c("pop","lat","lon")
# pop_coords <- pop_coords[match(pop_levels, pop_coords$pop),]
# rownames(pop_coords) <- pop_coords$pop

gc_dist <- marine_dist
# Average gc_dist within population pairs
pops_all <- pop_levels
n_pops   <- length(pops_all)

geo_dist_pop <- matrix(NA, n_pops, n_pops,
                        dimnames=list(pops_all, pops_all))
for(i in 1:n_pops){
  for(j in 1:n_pops){
    idx_i <- which(pop_factor == pops_all[i])
    idx_j <- which(pop_factor == pops_all[j])
    geo_dist_pop[i,j] <- mean(gc_dist[idx_i, idx_j])
  }
}
diag(geo_dist_pop) <- 0
cat("Population geographic distance matrix (km):\n")
print(round(geo_dist_pop, 0))

## ── 2. Environmental distance matrices ───────────────────────────────────────

env_data_v2 <- read.csv("RA_environmental_data_v2.csv")
env_data_v2 <- env_data_v2[match(coords$ID, env_data_v2$ID),]

# Final RDA variables — set1
env_vars <- c("ph_mean","thetao_range","phyc_range","sws_mean")
env_sub  <- env_data_v2[, c("ID", env_vars)]

# Population mean environmental values
env_pop <- matrix(NA, length(pop_levels), length(env_vars),
                   dimnames=list(pop_levels, env_vars))
for(p in pop_levels){
  idx <- which(pop_factor == p)
  env_pop[p,] <- colMeans(env_sub[match(coords$ID[idx],
                                         env_sub$ID), env_vars],
                           na.rm=TRUE)
}

# Standardise (mean=0, SD=1)
env_pop_scaled <- scale(env_pop)

# Shift for Bray-Curtis (requires non-negative values)
env_pop_shifted <- env_pop_scaled - min(env_pop_scaled) + 0.001

# Three distance metrics
env_dist <- list(
  BrayCurtis = as.matrix(vegdist(env_pop_shifted, method="bray")),
  Manhattan  = as.matrix(vegdist(env_pop_scaled,  method="manhattan")),
  Canberra   = as.matrix(vegdist(env_pop_scaled,  method="canberra"))
)
cat("Environmental distance matrices computed\n")

## ── 3. Linearized FST matrices ───────────────────────────────────────────────

if (!exists("fst_results")) load("RA_set1_FST_results.RData")

linearize_fst <- function(fst_mat){
  fst_mat[fst_mat < 0] <- 0
  return(fst_mat / (1 - fst_mat))
}

fst_lin <- lapply(fst_results, function(x){
  linearize_fst(as.matrix(x$fst))
})

## ── 4. Extract pairwise vectors ───────────────────────────────────────────────

pairs_idx <- which(lower.tri(geo_dist_pop), arr.ind=TRUE)
n_pairs   <- nrow(pairs_idx)
cat("Number of pairwise comparisons:", n_pairs, "\n")

geo_vec  <- geo_dist_pop[pairs_idx]
env_vecs <- lapply(env_dist, function(d) d[pairs_idx])

## ── 5. Permutation test function ─────────────────────────────────────────────

perm_r2 <- function(y, x, n_perm=10000){
  obs   <- cor(y, x, use="complete.obs")^2
  perms <- replicate(n_perm, cor(sample(y), x,
                                  use="complete.obs")^2)
  p_val <- mean(perms >= obs)
  return(list(r2=obs, p=p_val))
}

## ── 6. GLS models — all datasets, all env metrics ────────────────────────────

ibd_ibe_results <- data.frame()
best_env_metric <- list()

for(ds_name in names(fst_lin)){
  cat("\n===", ds_name, "===\n")
  gen_vec <- fst_lin[[ds_name]][pairs_idx]

  # IBD model
  ibd_mod  <- gls(gen_vec ~ geo_vec, method="ML")
  ibd_perm <- perm_r2(gen_vec, geo_vec)
  cat("IBD — R²:", round(ibd_perm$r2,3),
      "p:", round(ibd_perm$p,4), "\n")

  # IBE models — three metrics, select best by AIC
  ibe_aic <- numeric(3)
  names(ibe_aic) <- names(env_dist)
  for(metric in names(env_dist)){
    env_v          <- env_vecs[[metric]]
    ibe_mod        <- gls(gen_vec ~ env_v, method="ML")
    ibe_aic[metric] <- AIC(ibe_mod)
  }
  best_metric             <- names(which.min(ibe_aic))
  best_env_metric[[ds_name]] <- best_metric
  cat("Best environmental metric:", best_metric,
      "(AIC:", round(min(ibe_aic),1), ")\n")

  # IBE with best metric
  env_best <- env_vecs[[best_metric]]
  ibe_mod  <- gls(gen_vec ~ env_best, method="ML")
  ibe_perm <- perm_r2(gen_vec, env_best)
  cat("IBE — R²:", round(ibe_perm$r2,3),
      "p:", round(ibe_perm$p,4), "\n")

  # Combined IBD+IBE
  comb_mod <- gls(gen_vec ~ geo_vec + env_best, method="ML")

  # AIC comparison
  aic_ibd  <- AIC(ibd_mod)
  aic_ibe  <- AIC(ibe_mod)
  aic_comb <- AIC(comb_mod)
  best_mod <- c("IBD","IBE","Combined")[which.min(c(aic_ibd,
                                                     aic_ibe,
                                                     aic_comb))]
  cat("Best overall model:", best_mod, "\n")
  cat("ΔAIC IBD:", round(aic_ibd-min(c(aic_ibd,aic_ibe,aic_comb)),1),
      "IBE:", round(aic_ibe-min(c(aic_ibd,aic_ibe,aic_comb)),1),
      "Combined:", round(aic_comb-min(c(aic_ibd,aic_ibe,aic_comb)),1),
      "\n")

  ibd_ibe_results <- rbind(ibd_ibe_results, data.frame(
    Dataset      = ds_name,
    IBD_R2       = round(ibd_perm$r2, 3),
    IBD_p        = round(ibd_perm$p,  4),
    IBE_R2       = round(ibe_perm$r2, 3),
    IBE_p        = round(ibe_perm$p,  4),
    Best_env     = best_metric,
    AIC_IBD      = round(aic_ibd,  1),
    AIC_IBE      = round(aic_ibe,  1),
    AIC_Combined = round(aic_comb, 1),
    ΔAIC_IBD     = round(aic_ibd  - min(c(aic_ibd,aic_ibe,aic_comb)), 1),
    ΔAIC_IBE     = round(aic_ibe  - min(c(aic_ibd,aic_ibe,aic_comb)), 1),
    ΔAIC_Comb    = round(aic_comb - min(c(aic_ibd,aic_ibe,aic_comb)), 1),
    Best_model   = best_mod
  ))
}

# Print full results table
cat("\n── Full IBD/IBE Results Table ──\n")
print(ibd_ibe_results)
cat("\n── AIC Model Comparison ──\n")
print(ibd_ibe_results[, c("Dataset","AIC_IBD","AIC_IBE",
                            "AIC_Combined","ΔAIC_IBD",
                            "ΔAIC_IBE","ΔAIC_Comb","Best_model")])

write.csv(ibd_ibe_results, "RA_set1_IBD_IBE_results.csv",
          row.names=FALSE)

## ── 7. Rerun IBE with Canberra for all datasets (consistent metric) ──────────

env_vec_canberra <- env_dist[["Canberra"]][pairs_idx]

ibd_ibe_canberra <- data.frame()
for(ds_name in names(fst_lin)){
  gen_vec  <- fst_lin[[ds_name]][pairs_idx]
  ibd_perm <- perm_r2(gen_vec, geo_vec)
  ibe_perm <- perm_r2(gen_vec, env_vec_canberra)
  ibd_ibe_canberra <- rbind(ibd_ibe_canberra, data.frame(
    Dataset = ds_name,
    IBD_R2  = round(ibd_perm$r2, 3),
    IBD_p   = round(ibd_perm$p,  4),
    IBE_R2  = round(ibe_perm$r2, 3),
    IBE_p   = round(ibe_perm$p,  4)
  ))
}
cat("\n── Canberra IBD/IBE Results ──\n")
print(ibd_ibe_canberra)
write.csv(ibd_ibe_canberra,
          "RA_set1_IBD_IBE_results_canberra.csv",
          row.names=FALSE)

## ── 8. Plots ──────────────────────────────────────────────────────────────────

ds_levels <- c("neutral","adaptive","highconf","enviro")
ds_labels <- c("Neutral Loci","All Adaptive","High Confidence","Environmental")
ds_cols   <- c("#4472C4","#ED7D31","#70AD47","#FF0000")
ds_shapes <- c(16, 15, 17, 4)

# Build plot dataframes
plot_ibd <- data.frame()
plot_ibe <- data.frame()

for(ds_name in names(fst_lin)){
  gen_vec <- fst_lin[[ds_name]][pairs_idx]
  plot_ibd <- rbind(plot_ibd, data.frame(
    Dataset = ds_name,
    Geo     = geo_vec,
    FST     = gen_vec
  ))
  plot_ibe <- rbind(plot_ibe, data.frame(
    Dataset = ds_name,
    Env     = env_vec_canberra,
    FST     = gen_vec
  ))
}

plot_ibd$Dataset <- factor(plot_ibd$Dataset,
                             levels=ds_levels, labels=ds_labels)
plot_ibe$Dataset <- factor(plot_ibe$Dataset,
                             levels=ds_levels, labels=ds_labels)

# R² subtitle labels
ibd_sub <- paste0("IBD R² values: ",
  paste(paste0(ds_labels, " = ",
               ibd_ibe_canberra$IBD_R2), collapse=", "))
ibe_sub <- paste0("IBE R² values: ",
  paste(paste0(ds_labels, " = ",
               ibd_ibe_canberra$IBE_R2), collapse=", "))

# IBD plot
p_ibd <- ggplot(plot_ibd, aes(x=Geo, y=FST,
                               colour=Dataset, shape=Dataset)) +
  geom_point(size=2.5, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8,
              aes(fill=Dataset), alpha=0.15) +
  scale_colour_manual(values=ds_cols, labels=ds_labels,
                       name="Dataset") +
  scale_fill_manual(values=ds_cols,   labels=ds_labels,
                     name="Dataset") +
  scale_shape_manual(values=ds_shapes,labels=ds_labels,
                      name="Dataset") +
  xlab("Geographic Distance (km)") +
  ylab(expression("Genetic Distance (F"[ST]*"/1-F"[ST]*")")) +
  labs(title="Isolation by Distance across genomic categories",
       subtitle=ibd_sub) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.title       = element_text(size=11, face="bold"),
        plot.subtitle    = element_text(size=7,  colour="grey30"))

# IBE plot
p_ibe <- ggplot(plot_ibe, aes(x=Env, y=FST,
                               colour=Dataset, shape=Dataset)) +
  geom_point(size=2.5, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8,
              aes(fill=Dataset), alpha=0.15) +
  scale_colour_manual(values=ds_cols, labels=ds_labels,
                       name="Dataset") +
  scale_fill_manual(values=ds_cols,   labels=ds_labels,
                     name="Dataset") +
  scale_shape_manual(values=ds_shapes,labels=ds_labels,
                      name="Dataset") +
  xlab("Environmental Distance (Canberra)") +
  ylab(expression("Genetic Distance (F"[ST]*"/1-F"[ST]*")")) +
  labs(title="Isolation by Environment across genomic categories",
       subtitle=ibe_sub) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.title       = element_text(size=11, face="bold"),
        plot.subtitle    = element_text(size=7,  colour="grey30"))

pdf("RA_set1_IBD_plot.pdf", width=9, height=7)
  print(p_ibd)
dev.off()

pdf("RA_set1_IBE_plot.pdf", width=9, height=7)
  print(p_ibe)
dev.off()

saveRDS(p_ibd, "RA_set1_IBD_plot.rds")
saveRDS(p_ibe, "RA_set1_IBE_plot.rds")

save(geo_dist_pop, env_dist, fst_lin,
     ibd_ibe_results, ibd_ibe_canberra,
     best_env_metric,
     file="RA_set1_IBD_IBE_results.RData")

cat("\nAll IBD/IBE analyses complete\n")

# Updated IBD and IBE plots

ds_labels_new <- c("Neutral","Adaptive","High confidence","Environmental")

plot_ibd$Dataset <- factor(as.character(plot_ibd$Dataset),
                            levels=c("Neutral Loci","All Adaptive",
                                     "High Confidence","Environmental"),
                            labels=ds_labels_new)
plot_ibe$Dataset <- factor(as.character(plot_ibe$Dataset),
                            levels=c("Neutral Loci","All Adaptive",
                                     "High Confidence","Environmental"),
                            labels=ds_labels_new)

# Updated R² subtitle
ibd_sub <- paste0("IBD R² values: ",
  paste(paste0(ds_labels_new, " = ",
               ibd_ibe_canberra$IBD_R2), collapse=", "))
ibe_sub <- paste0("IBE R² values: ",
  paste(paste0(ds_labels_new, " = ",
               ibd_ibe_canberra$IBE_R2), collapse=", "))

# IBD plot
p_ibd <- ggplot(plot_ibd, aes(x=Geo, y=FST,
                               colour=Dataset, shape=Dataset)) +
  geom_point(size=2.5, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8,
              aes(fill=Dataset), alpha=0.15) +
  scale_colour_manual(values=ds_cols, labels=ds_labels_new,
                       name="Dataset") +
  scale_fill_manual(values=ds_cols,   labels=ds_labels_new,
                     name="Dataset") +
  scale_shape_manual(values=ds_shapes,labels=ds_labels_new,
                      name="Dataset") +
  xlab("Geographic distance (km)") +
  ylab(expression("Genetic distance (F"[ST]*")")) +
  labs(title=expression("Isolation by distance across genomic categories"~
                         "-"~italic("R. australiae")),
       subtitle=ibd_sub) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor  = element_blank(),
        legend.position   = "right",
        plot.title        = element_text(size=11, face="bold"),
        plot.subtitle     = element_text(size=9, colour="grey30"))

# IBE plot
p_ibe <- ggplot(plot_ibe, aes(x=Env, y=FST,
                               colour=Dataset, shape=Dataset)) +
  geom_point(size=2.5, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8,
              aes(fill=Dataset), alpha=0.15) +
  scale_colour_manual(values=ds_cols, labels=ds_labels_new,
                       name="Dataset") +
  scale_fill_manual(values=ds_cols,   labels=ds_labels_new,
                     name="Dataset") +
  scale_shape_manual(values=ds_shapes,labels=ds_labels_new,
                      name="Dataset") +
  xlab("Environmental distance") +
  ylab(expression("Genetic distance (F"[ST]*")")) +
  labs(title=expression("Isolation by environment across genomic categories"~
                         "-"~italic("R. australiae")),
       subtitle=ibe_sub) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor  = element_blank(),
        legend.position   = "right",
        plot.title        = element_text(size=11, face="bold"),
        plot.subtitle     = element_text(size=9, colour="grey30"))

# Save PDF and PNG
for(plot_obj in list(list(p=p_ibd, name="RA_set1_IBD"),
                      list(p=p_ibe, name="RA_set1_IBE"))){
  pdf(paste0(plot_obj$name, ".pdf"), width=9, height=7)
    print(plot_obj$p)
  dev.off()
  png(paste0(plot_obj$name, ".png"), width=9, height=7,
      units="in", res=300)
    print(plot_obj$p)
  dev.off()
}

saveRDS(p_ibd, "RA_set1_IBD_plot.rds")
saveRDS(p_ibe, "RA_set1_IBE_plot.rds")


# ======================================================================
# 9. Diversity statistics figure (neutral vs adaptive)
# Built from the section 1 and 2 outputs (CSV files)
# ======================================================================

library(dplyr)
library(tidyr)



n_snps <- c(Neutral = 24415, Adaptive = 3147)     # SNPs per dataset (private alleles as % of SNPs)

div <- bind_rows(
  read.csv("RA_set1_neutral_diversity_table.csv"),
  read.csv("RA_set1_adaptive_diversity_table.csv")) %>%
  select(Dataset, Population, N = n, Ho = HO, He = HE, FIS)

arpa <- read.csv("RA_set1_AR_PA_all_datasets.csv") %>%
  filter(Dataset %in% c("neutral", "adaptive")) %>%
  select(Dataset, Population, AR, PA)

div <- div %>%
  left_join(arpa, by = c("Dataset", "Population")) %>%
  mutate(Dataset  = factor(ifelse(Dataset == "neutral", "Neutral", "Adaptive"),
                           levels = c("Neutral", "Adaptive")),
         Location = factor(pop_display_RA[Population],
                           levels = c("Malaysia", "Mozambique", "Tanzania",
                                      "Réunion Island", "Seychelles")),
         NP       = PA / n_snps[as.character(Dataset)] * 100)

div_long <- div %>%
  select(Location, Dataset, Ho, He, AR, FIS, NP) %>%
  pivot_longer(c(Ho, He, AR, FIS, NP), names_to = "Statistic", values_to = "Value") %>%
  mutate(Statistic = factor(Statistic,
                            levels = c("Ho", "He", "AR", "FIS", "NP"),
                            labels = c("Observed heterozygosity (Ho)",
                                       "Expected heterozygosity (HE)",
                                       "Allelic richness (AR)",
                                       "FIS",
                                       "Private alleles (% of SNPs)")))

# Version 1: colours by dataset
p1 <- ggplot(div_long, aes(Location, Value, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Statistic, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Neutral" = "blue4", "Adaptive" = "brown2"), name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")
ggsave("diversity_neutral_vs_adaptive_datasetcolours.png", p1, width = 11, height = 7, dpi = 300)

# Version 2: population colours, adaptive in full colour and neutral in a lighter tint
pop_cols <- c("Malaysia"       = "#E75480",
              "Mozambique"     = "#0B6DB7",
              "Tanzania"       = "#E69F00",
              "Réunion Island" = "#D55E00",
              "Seychelles"     = "#009E73")
lighten <- function(col, amount = 0.55) {
  m <- col2rgb(col) / 255
  m <- m + (1 - m) * amount
  rgb(m[1, ], m[2, ], m[3, ])
}
fill_cols <- c(setNames(lighten(pop_cols), paste(names(pop_cols), "Neutral")),
               setNames(pop_cols,          paste(names(pop_cols), "Adaptive")))

p2 <- ggplot(div_long, aes(Location, Value, group = Dataset,
                           fill = paste(Location, Dataset))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Statistic, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())
ggsave("diversity_neutral_vs_adaptive.pdf", p2, width = 11, height = 7, device = cairo_pdf)
ggsave("diversity_neutral_vs_adaptive.png", p2, width = 11, height = 7, dpi = 300)

# ======================================================================
# 10. Sampling map
# ======================================================================

#Figure 1 map
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# ── Data ─────────────────────────────────────────────────────
world <- ne_countries(scale = "medium", returnclass = "sf")

# ── Population Colours ────────────────────────────────────────
pop_cols_RA <- c(
  "MZ" = "#0072B2",
  "TZ" = "#E69F00",
  "MD" = "#56B4E9",
  "RI" = "#D55E00",
  "SY" = "#009E73",
  "MA" = "#FF69B4",
  "AU" = "#8B008B"
)

# population code for each individual (used for the point colours)
coords$pop <- popmap$pop[match(coords$ID, popmap$ID)]

# ── Main map ──────────────────────────────────────────────────
focal_map <- ggplot() +
    geom_sf(data = world,
            fill = "gray95", color = "gray70", linewidth = 0.2) +
    
    # ── Sampling Dots Layer ───────────────────────────────────
    geom_point(data = coords, 
               aes(x = lon, y = lat, fill = pop), 
               shape = 21, color = "black", size = 3.5, stroke = 0.5) +
    
    # ── Custom Colour Scale ───────────────────────────────────
    scale_fill_manual(values = pop_cols_RA) +
    
    # ── Zoomed-Out Indo-Pacific Bounding Box ──────────────────
    coord_sf(xlim = c(20, 160), ylim = c(-50, 20), expand = FALSE) +
    
    annotation_scale(
        location   = "br",
        width_hint = 0.15,
        line_width = 0.3,
        text_cex   = 0.7
    ) +
    
    # ── Axis Formatting (20-degree intervals) ─────────────────
    scale_x_continuous(
        breaks = seq(20, 160, by = 20),
        labels = function(x) paste0(abs(x), "\u00b0E")
    ) +
    scale_y_continuous(
        breaks = seq(-40, 20, by = 20),
        labels = function(y) ifelse(y < 0,
                                    paste0(abs(y), "\u00b0S"),
                                    ifelse(y == 0, "0\u00b0", paste0(y, "\u00b0N")))
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        axis.title       = element_blank(),
        axis.text        = element_text(size = 9, color = "black"),
        text             = element_text(color = "black"),
        panel.border     = element_rect(color = "gray40", fill = NA,
                                        linewidth = 0.5),
        legend.position  = "none" # Legend removed
    )

# ── Export focal map ───────────────────────────────────────────
# The 12x6 dimensions perfectly balance the land shapes in this broad view
ggsave("RA_IndoPacific_zoomed_out.png",
       plot   = focal_map,
       width  = 12,
       height = 6,
       dpi    = 600)

ggsave("RA_IndoPacific_zoomed_out.pdf",
       plot   = focal_map,
       width  = 12,
       height = 6)

message("Done — Broad Indo-Pacific map saved successfully without legend")

