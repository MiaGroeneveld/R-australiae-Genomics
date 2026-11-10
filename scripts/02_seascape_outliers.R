# ============================================================
# Seascape Genomics Pipeline — outlier detection and four-dataset compilation
# Species: Rhynchobatus australiae
# ============================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(vcfR)       # VCF loading and manipulation
library(adegenet)   # genind object and population genetics
library(OutFLANK)   # differentiation-based outlier detection
library(pcadapt)    # structure-based outlier detection
library(qvalue)     # FDR control via Storey q-values
library(ggplot2)    # visualisation
library(readxl)     # Excel file loading (not needed if all CSV)

# ── Parameters — define once, use throughout ──────────────────────────────────
# All key analytical values defined here
# Change values here only — not throughout the script

VCF_FILE      <- "RApop_set1.snps.vcf"  # set1: MAC>=3, n=71, 27,562 SNPs
POPMAP_FILE   <- "RApopmap.txt"          # two-column: ID pop
COORDS_FILE   <- "RA_coordinates.csv"    # comma-delimited: ID, lat, lon
N_IND         <- 71       # total individuals
N_SNP         <- 27562    # total SNPs
K_PCADAPT     <- 4        # optimal K from scree plot + GIF inspection
K_LFMM        <- 4        # consistent with pcadapt K — see Section 9
FDR_PCADAPT   <- 0.05     # Storey q-value threshold for pcadapt outliers
FDR_LFMM      <- 0.01     # FDR threshold for LFMM outliers
RDA_SD        <- 3        # SD threshold for RDA candidate loci (Forester 2018)
N_BOOT_FST    <- 1000     # bootstrap replicates for pairwise FST
N_BOOT_DIV    <- 500      # bootstrap replicates for diversity statistics
N_PERM_IBD    <- 10000    # permutations for IBD/IBE GLS test
N_PERM_RDA    <- 999      # permutations for RDA ANOVA significance test
MISS_THRESH   <- 0.20     # individual missingness threshold (20%)

# Population order
# Mainland Africa (MZ, TZ) → Islands (MD, RI, SY) → Indo-Pacific (MA, AU)
POP_ORDER <- c("MZ","TZ","MD","RI","SY","MA","AU")
N_POPS    <- 7

# Population display names for figures
POP_DISPLAY <- c(
  "MZ" = "Mozambique",
  "TZ" = "Tanzania",
  "MD" = "Madagascar",
  "RI" = "Réunion Island",
  "SY" = "Seychelles",
  "MA" = "Malaysia",
  "AU" = "Australia"
)

# ── File verification — run before any analysis ───────────────────────────────
required_files <- c(VCF_FILE, POPMAP_FILE, COORDS_FILE)
missing_files  <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
} else {
  cat("All input files found ✓\n")
}

# ============================================================
# SECTION 1: Load raw inputs
# ============================================================

# ── 1.1 Popmap ────────────────────────────────────────────────────────────────
# Two-column tab-delimited file: individual ID and population code
# Population codes: MZ=Mozambique, TZ=Tanzania, MD=Madagascar,
#                   RI=Réunion Island, SY=Seychelles, MA=Malaysia, AU=Australia
popmap <- read.table(POPMAP_FILE,
                     header    = FALSE,
                     col.names = c("ID","pop"))
popmap$pop <- factor(popmap$pop, levels = POP_ORDER)

cat("Popmap loaded:", nrow(popmap), "individuals\n")
cat("Populations:", levels(popmap$pop), "\n")
cat("Population sizes:\n"); print(table(popmap$pop))
# Expected: TZ=8, MD=2, SY=12, AU=2, RI=9, MZ=28, MA=10
# Total = 71 ✓

stopifnot(nrow(popmap) == N_IND)
cat("Popmap n check passed ✓\n")

# ── 1.2 Coordinates ───────────────────────────────────────────────────────────
# Comma-delimited CSV: ID, lat, lon
coords <- read.csv(COORDS_FILE,
                   header           = TRUE,
                   sep              = ",",
                   stringsAsFactors = FALSE)

# Rename columns if needed — expected: ID, lat, lon
colnames(coords) <- c("ID","lat","lon")

# Convert to numeric — handles any decimal comma issues
coords$lat <- as.numeric(gsub(",",".", trimws(coords$lat)))
coords$lon <- as.numeric(gsub(",",".", trimws(coords$lon)))

# Match to popmap order
coords <- coords[match(popmap$ID, coords$ID), ]
coords$pop <- popmap$pop[match(coords$ID, popmap$ID)]

cat("Coordinates loaded:", nrow(coords), "individuals\n")
cat("Missing coordinates:\n")
print(coords[is.na(coords$lat) | is.na(coords$lon), "ID"])

# Verify order matches popmap
stopifnot(all(coords$ID == popmap$ID))
cat("Coordinates-popmap order consistent ✓\n")

# Check Réunion Island coordinates are offshore (not on land)
# Réunion Island approximate land bbox: lat -21.4 to -20.8, lon 55.2 to 55.8
ri_coords <- coords[coords$pop == "RI", c("ID","lat","lon")]
cat("Réunion Island coordinates (verify offshore):\n")
print(ri_coords)

# ── 1.3 VCF and genotype matrix ───────────────────────────────────────────────
cat("\nLoading VCF...\n")
obj.vcfR <- read.vcfR(VCF_FILE, verbose = FALSE)

cat("VCF loaded:", nrow(obj.vcfR@fix), "loci |",
    ncol(obj.vcfR@gt) - 1, "individuals\n")

# Hard stops — VCF must match expected dimensions exactly
stopifnot(ncol(obj.vcfR@gt) - 1 == N_IND)
stopifnot(nrow(obj.vcfR@fix) == N_SNP)
cat("VCF dimension check passed ✓\n")

# Verify VCF individual order matches popmap
stopifnot(all(colnames(obj.vcfR@gt)[-1] == popmap$ID))
cat("VCF-popmap order consistent ✓\n")

# Extract raw genotypes (loci x individuals matrix)
geno       <- extract.gt(obj.vcfR)
position   <- getPOS(obj.vcfR)
chromosome <- getCHROM(obj.vcfR)

# Locus names: CHR_POS format
# VCF POS is already 1-based — no offset adjustment needed
locus_names_full <- paste(chromosome, position, sep = "_")

# Convert GT strings to dosage (0/1/2)
# 0/0 or 0|0 = homozygous reference = 0 copies of alt allele
# 0/1 etc    = heterozygous          = 1 copy of alt allele
# 1/1 or 1|1 = homozygous alternate  = 2 copies of alt allele
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0","0|0")] <- 0
G[geno %in% c("0/1","1/0","1|0","0|1")] <- 1
G[geno %in% c("1/1","1|1")] <- 2

cat("Missing data rate before imputation:",
    round(mean(is.na(G)) * 100, 2), "%\n")
# Expected: ~6.32% — consistent with 80% per-population genotyping threshold
# Values >20% per individual flagged as potential quality issue

# Transpose: individuals as rows, loci as columns
# Required for RDA, LFMM, sNMF input format
G <- t(G)
rownames(G) <- colnames(obj.vcfR@gt)[-1]  # individual IDs
colnames(G) <- locus_names_full             # locus names CHR_POS

# Verify dimensions: 71 individuals x 27,562 loci
cat("G dimensions:", dim(G), "\n")
stopifnot(nrow(G) == N_IND)
stopifnot(ncol(G) == N_SNP)

# ── 1.4 Mean imputation of missing genotypes ──────────────────────────────────
# Most common genotype per locus
# used for RDA, LFMM, sNMF
# Raw G matrix retained for diversity stats and FST calculation
G_imp <- apply(G, 2, function(x) {
  miss <- is.na(x)
  if (any(miss)) {
    x[miss] <- as.numeric(
      names(sort(table(x[!miss]), decreasing = TRUE)[1]))
  }
  x
})
rownames(G_imp) <- rownames(G)
colnames(G_imp) <- locus_names_full

stopifnot(sum(is.na(G_imp)) == 0)
cat("Imputation complete — missing data remaining:",
    sum(is.na(G_imp)), "✓\n")

# Verify dimensions unchanged after imputation
stopifnot(nrow(G_imp) == N_IND)
stopifnot(ncol(G_imp) == N_SNP)
cat("G_imp dimensions:", dim(G_imp), "\n")

# ── 1.5 Build genind for population genetics analyses ─────────────────────────
gind_full <- vcfR2genind(obj.vcfR)

# Verify and assign locus names — CHR_POS format
# vcfR2genind uses CHR:POS format internally — convert to CHR_POS
locNames(gind_full) <- locus_names_full

# Verify locus names match G_imp
stopifnot(all(locNames(gind_full) == colnames(G_imp)))
cat("Genind-G_imp locus name consistency ✓\n")

# Assign population factor in correct order
pop_factor       <- popmap$pop[match(indNames(gind_full), popmap$ID)]
pop(gind_full)   <- factor(pop_factor, levels = POP_ORDER)

cat("Genind built:", nInd(gind_full), "individuals |",
    nLoc(gind_full), "loci\n")
cat("Population sizes:\n"); print(table(pop(gind_full)))
# Expected: MZ=28, TZ=8, MD=2, RI=9, SY=12, MA=10, AU=2
# Doesn't follow this order; reorder before producing results

# ── 1.6 Final consistency check ───────────────────────────────────────────────
# All individual orders must be identical across all objects
stopifnot(all(rownames(G_imp) == popmap$ID))
stopifnot(all(coords$ID == popmap$ID))
stopifnot(all(indNames(gind_full) == popmap$ID))
cat("All individual orders consistent across objects ✓\n")

# ── 1.7 Save section 1 objects ────────────────────────────────────────────────
save(popmap, coords, obj.vcfR, geno, G, G_imp,
     locus_names_full, position, chromosome,
     gind_full, pop_factor,
     file = "RA_set1_section1_objects.RData")
cat("Section 1 complete — objects saved to RA_set1_section1_objects.RData\n")
cat("Objects saved: popmap, coords, obj.vcfR, geno, G, G_imp,\n")
cat("               locus_names_full, position, chromosome,\n")
cat("               gind_full, pop_factor\n")

# ============================================================
# SECTION 2: OutFLANK
# Fits chi-squared distribution to trimmed FST distribution
# NOTE: zero outliers expected — strong neutral population
# structure inflates FST null distribution beyond chi-sq fit
# ============================================================

# ── Using objects from Section 1 ──────────────────────────────────────────────
# Required: G, locus_names_full, popmap
# If running standalone, uncomment:
# library(vcfR)
# obj.vcfR   <- read.vcfR("RApop_set1.snps.vcf", verbose=FALSE)
# geno       <- extract.gt(obj.vcfR)
# position   <- getPOS(obj.vcfR)
# chromosome <- getCHROM(obj.vcfR)
# G <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
# G[geno %in% c("0/0","0|0")] <- 0
# G[geno %in% c("0/1","1/0","1|0","0|1")] <- 1
# G[geno %in% c("1/1","1|1")] <- 2
# G <- t(G)
# locus_names_full <- paste(getCHROM(obj.vcfR), getPOS(obj.vcfR), sep="_")
# popmap <- read.table("RApopmap.txt", header=FALSE)
# colnames(popmap) <- c("ID","pop")

library(OutFLANK)

G_of <- G
G_of[is.na(G_of)] <- 9

# MakeDiploidFSTMat: SNPmat = loci x individuals (transposed)
# popNames must be character vector — not factor
my_fst <- MakeDiploidFSTMat(
  SNPmat     = t(G_of),
  locusNames = locus_names_full,
  popNames   = as.character(popmap$pop)
)

cat("Mean FST:", round(mean(my_fst$FST, na.rm = TRUE), 4), "\n")

# Run OutFLANK — first argument positional
# NumberOfSamples = 7 populations for RA
out <- OutFLANK(
  my_fst,
  LeftTrimFraction  = 0.05,
  RightTrimFraction = 0.05,
  Hmin              = 0.1,
  NumberOfSamples   = N_POPS,
  qthreshold        = 0.05
)

# FST distribution plot — check chi-sq fit
# Poor fit confirms justification for exclusion
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_outflank_fst_distribution.pdf",
                               width = 7, height = 6)
  else png("RA_set1_outflank_fst_distribution.png",
           width = 7, height = 6, units = "in", res = 300)
  OutFLANKResultsPlotter(out,
                          withOutliers = TRUE,
                          NoCorr       = TRUE,
                          Hmin         = 0.1,
                          binwidth     = 0.005,
                          titletext    = "")
  dev.off()
}
cat("OutFLANK FST distribution plot saved\n")

# Extract outliers
P1 <- pOutlierFinderChiSqNoCorr(
  my_fst,
  Fstbar     = out$FSTNoCorrbar,
  dfInferred = out$dfInferred,
  qthreshold = 0.05,
  Hmin       = 0.1
)

outliers_outflank <- P1$LocusName[which(P1$OutlierFlag == TRUE)]
cat("OutFLANK outliers (q<0.05):", length(outliers_outflank), "\n")
# Expected: 0 — excluded from four-dataset framework

# Save
save(G_of, my_fst, out, P1, outliers_outflank,
     file = "RA_set1_outflank_results.RData")
cat("Section 2 complete — OutFLANK results saved\n")


# ============================================================
# SECTION 3: pcadapt
# Structure-based individual-level outlier detection
# K=4 selected: lowest GIF=1.286, corroborated by scree plot
# Threshold: Storey q-value < 0.05 (Storey & Tibshirani 2003)
# Result: 1 270 outliers (4.6% of 27 562 SNPs)
# ============================================================

# ── Using objects from Section 1 ──────────────────────────────────────────────
# Required: locus_names_full, popmap
# pcadapt reads VCF directly — most reliable input format
# Individual order in VCF confirmed = popmap order (Section 1) ✓
# If running standalone, uncomment:
# popmap <- read.table("RApopmap.txt", header=FALSE)
# colnames(popmap) <- c("ID","pop")
# locus_names_full <- paste(chromosome, position, sep="_")

library(pcadapt)
library(qvalue)

ra_pcadapt <- read.pcadapt("RApop_set1.snps.vcf", type = "vcf")

# ── GIF across K=1-20 ────────────────────────────────────────────────────────
# Select K where GIF is lowest and closest to 1
# Corroborated by visual inspection of scree plot
gif_values <- data.frame(K = 1:20, GIF = NA)
for (k in 1:20) {
  x_temp <- pcadapt(input = ra_pcadapt, K = k)
  gif_values$GIF[k] <- x_temp$gif
}
cat("GIF values K=1-20:\n")
print(gif_values, digits = 6)
cat("Optimal K (lowest GIF):",
    gif_values$K[which.min(gif_values$GIF)],
    "GIF =", round(min(gif_values$GIF), 5), "\n")

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_pcadapt_gif.pdf",
                               width = 8, height = 5)
  else png("RA_set1_pcadapt_gif.png",
           width = 8, height = 5, units = "in", res = 300)
  plot(gif_values$K, gif_values$GIF,
       type = "b", pch = 19,
       xlab = "K", ylab = "GIF",
       main = "GIF across K — R. australiae set1",
       las  = 1)
  abline(h = 1, col = "grey50", lty = 2)
  abline(v = K_PCADAPT, col = "red", lty = 2)
  legend("topright", legend = paste0("Selected K=", K_PCADAPT),
         col = "red", lty = 2, bty = "n", cex = 0.9)
  dev.off()
}

# ── Scree plot ────────────────────────────────────────────────────────────────
x_K20 <- pcadapt(input = ra_pcadapt, K = 20)

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_pcadapt_screeplot.pdf",
                               width = 8, height = 5)
  else png("RA_set1_pcadapt_screeplot.png",
           width = 8, height = 5, units = "in", res = 300)
  plot(x_K20, option = "screeplot",
       main = "R. australiae set1 — pcadapt scree plot")
  dev.off()
}

# ── Final run at K=4 ─────────────────────────────────────────────────────────
x <- pcadapt(input = ra_pcadapt, K = K_PCADAPT)
cat("pcadapt GIF at K=", K_PCADAPT, ":", round(x$gif, 5), "\n")
# Expected: 1.28573
# NOTE: report this in methods — NOT 2.033 which is LFMM GIF

# ── p-value histogram ─────────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_pcadapt_pval_hist.pdf",
                               width = 7, height = 5)
  else png("RA_set1_pcadapt_pval_hist.png",
           width = 7, height = 5, units = "in", res = 300)
  hist(x$pvalues,
       xlab   = "p-values",
       main   = "pcadapt p-value distribution — R. australiae set1",
       breaks = 50, col = "orange")
  dev.off()
}

# ── Q-Q plot ──────────────────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_pcadapt_qqplot.pdf",
                               width = 6, height = 6)
  else png("RA_set1_pcadapt_qqplot.png",
           width = 6, height = 6, units = "in", res = 300)
  plot(x, option = "qqplot",
       main = "R. australiae set1 — pcadapt Q-Q (K=4)")
  dev.off()
}

# ── Outlier detection — Storey q-value < 0.05 ────────────────────────────────
qval                   <- qvalue(x$pvalues)$qvalues
outliers_pcadapt       <- which(qval < FDR_PCADAPT)
outliers_pcadapt_names <- locus_names_full[outliers_pcadapt]

cat("pcadapt outliers (q <", FDR_PCADAPT, "):",
    length(outliers_pcadapt_names), "\n")
cat("Percentage:",
    round(length(outliers_pcadapt_names)/N_SNP*100, 1), "%\n")
# Expected: 1 270 outliers (4.6%)

# ── Manhattan plot with FDR threshold line ────────────────────────────────────
# fdr_line = -log10 of maximum p-value among q<0.05 outliers
fdr_threshold <- max(x$pvalues[qval < FDR_PCADAPT], na.rm = TRUE)
fdr_line      <- -log10(fdr_threshold)
cat("FDR threshold (-log10 p):", round(fdr_line, 3), "\n")

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_pcadapt_manhattan.pdf",
                               width = 10, height = 5)
  else png("RA_set1_pcadapt_manhattan.png",
           width = 10, height = 5, units = "in", res = 300)
  plot(seq_along(x$pvalues),
       -log10(x$pvalues),
       pch      = 19, cex = 0.3,
       col      = adjustcolor("grey40", alpha.f = 0.6),
       xlab     = "SNP index",
       ylab     = expression(-log[10](italic(p)~value)),
       main     = "R. australiae set1 — pcadapt Manhattan (K=4)",
       cex.lab  = 1.2,
       cex.axis = 1.0,
       las      = 1,
       bty      = "l")
  abline(h   = fdr_line,
         col = "#87c1f7",
         lty = 2,
         lwd = 1.5)
  legend("topright",
         legend = "FDR (q < 0.05)",
         col    = "#87c1f7",
         lty    = 2,
         lwd    = 1.5,
         bty    = "n",
         cex    = 0.9)
  dev.off()
}
cat("Manhattan plot saved\n")

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(data.frame(locus = outliers_pcadapt_names),
          "RA_set1_outliers_pcadapt.csv",
          row.names = FALSE)

save(x, gif_values, qval, outliers_pcadapt, outliers_pcadapt_names,
     file = "RA_set1_pcadapt_results.RData")
cat("Section 3 complete — pcadapt results saved\n")
cat("pcadapt outliers:", length(outliers_pcadapt_names),
    "(", round(length(outliers_pcadapt_names)/N_SNP*100, 1), "%)\n")

# ============================================================
# SECTION 4: Environmental Variable Extraction
# Bio-ORACLE v3.0 (Assis et al. 2024) — biooracler package
# MARSPEC (Sbrocco & Barber 2013) — sdmpredictors package
# 11 starting variables
# Chlorophyll excluded: r=1.00 with phytoplankton (phyc)
# MARSPEC bathymetry used instead of Bio-ORACLE v3 bathymetry:
#   Bio-ORACLE v3 cell-mean depth averaging produced anomalous
#   values (-909m) for shallow coastal Réunion Island samples
#   due to steep volcanic submarine slopes
# All Bio-ORACLE variables use surface layers (depthsurf)
# Extraction at individual sampling coordinates
#
# NOTE: extraction already completed - load CSV if available
# If rerunning extraction, ensure corrected Réunion Island
# coordinates are used (RA_coordinates.csv — already corrected)
# ============================================================

library(sdmpredictors)
library(biooracler)
library(raster)

# ── Option A: Load existing CSV (recommended) ─────────────────────────────────
# Use this if extraction has already been completed
env_data_11 <- read.csv("RA_environmental_data_11vars.csv",
                         stringsAsFactors = FALSE)
env_data_11 <- env_data_11[match(popmap$ID, env_data_11$ID), ]
env_data_11$pop <- as.character(popmap$pop)

cat("Environmental data loaded:", nrow(env_data_11), "individuals\n")
cat("Variables:", colnames(env_data_11), "\n")
cat("Missing values:\n"); print(colSums(is.na(env_data_11)))
stopifnot(all(env_data_11$ID == popmap$ID))
cat("Env data-popmap order consistent ✓\n")

# ── Option B: Re-extract from scratch (uncomment if needed) ──────────────────
# Requires internet connection and biooracler package
# Use corrected RA_coordinates.csv (Réunion Island corrections applied)
#
# Bio-ORACLE v3.0 point extraction function
# extract_biooracle_points <- function(dataset_id, variables, coords) {
#   results <- data.frame(ID = coords$ID)
#   for (var in variables) {
#     vals <- numeric(nrow(coords))
#     for (i in 1:nrow(coords)) {
#       tryCatch({
#         tiny <- list(
#           latitude  = c(coords$lat[i] - 0.05, coords$lat[i] + 0.05),
#           longitude = c(coords$lon[i] - 0.05, coords$lon[i] + 0.05))
#         d <- download_layers(dataset_id  = dataset_id,
#                              variables   = var,
#                              constraints = tiny,
#                              fmt         = "csv")
#         vals[i] <- mean(d[[var]], na.rm = TRUE)
#         cat("Done:", coords$ID[i], var, "=", round(vals[i], 3), "\n")
#       }, error = function(e) {
#         vals[i] <<- NA
#         cat("NA:", coords$ID[i], var, "\n")
#       })
#     }
#     results[[var]] <- vals
#   }
#   return(results)
# }
#
# # MARSPEC via sdmpredictors
# # bathymetry = MS_bathy_5m; SSS mean = MS_biogeo08_sss_mean_5m;
# # SSS range = MS_biogeo11_sss_range_5m
# options(sdmpredictors_datadir = ".")
# marspec_rasters <- load_layers(c("MS_biogeo08_sss_mean_5m",
#                                   "MS_biogeo11_sss_range_5m",
#                                   "MS_bathy_5m"))
# marspec_values  <- raster::extract(marspec_rasters,
#                                     coords[, c("lon","lat")])
# marspec_df      <- data.frame(ID = coords$ID, marspec_values)
# cat("MARSPEC extracted\n")
# print(colSums(is.na(marspec_df)))
#
# # Bio-ORACLE variables
# # Variable selection rationale:
# # ph_mean     — ocean pH; physiological stress, acid-base regulation
# # thetao_mean — sea surface temperature mean; metabolic rate, thermal tolerance
# # thetao_range— temperature seasonality; thermal variability stress
# # o2_ltmin    — dissolved oxygen minimum; hypoxia tolerance
# # phyc_mean   — phytoplankton carbon mean; prey availability proxy
# # phyc_range  — phytoplankton carbon seasonality; productivity variability
# # sws_mean    — current velocity mean; dispersal potential, energy expenditure
# # sws_range   — current velocity range; flow variability
#
# thetao_df <- extract_biooracle_points(
#   "thetao_baseline_2000_2019_depthsurf",
#   c("thetao_mean","thetao_range"), coords)
# o2_df     <- extract_biooracle_points(
#   "o2_baseline_2000_2018_depthsurf",
#   c("o2_ltmin"), coords)
# phyc_df   <- extract_biooracle_points(
#   "phyc_baseline_2000_2020_depthsurf",
#   c("phyc_mean","phyc_range"), coords)
# sws_df    <- extract_biooracle_points(
#   "sws_baseline_2000_2019_depthsurf",
#   c("sws_mean","sws_range"), coords)
# ph_df     <- extract_biooracle_points(
#   "ph_baseline_2000_2018_depthsurf",
#   c("ph_mean"), coords)
#
# # Combine all variables
# env_data_11 <- Reduce(function(x,y) merge(x,y, by="ID"),
#                        list(marspec_df, thetao_df, o2_df,
#                             phyc_df, sws_df, ph_df))
# env_data_11$pop <- popmap$pop[match(env_data_11$ID, popmap$ID)]
#
# cat("Dimensions:", nrow(env_data_11), "x", ncol(env_data_11), "\n")
# cat("Missing values:\n"); print(colSums(is.na(env_data_11)))
#
# write.csv(env_data_11, "RA_environmental_data_11vars.csv",
#           row.names = FALSE)
# cat("Environmental data saved\n")

# ── Environmental correlation heatmap ─────────────────────────────
# Pearson correlation matrix — all 11 variables
# Justifies variable selection by showing collinearity structure
library(corrplot)

env_vars_11 <- c("ph_mean","thetao_mean","thetao_range",
                  "o2_ltmin","phyc_mean","phyc_range",
                  "sws_mean","sws_range",
                  "MS_bathy_5m",
                  "MS_biogeo08_sss_mean_5m",
                  "MS_biogeo11_sss_range_5m")

# Clean display names for figure
clean_names <- c(
  "ph_mean"                  = "pH (mean)",
  "thetao_mean"              = "Temp (mean)",
  "thetao_range"             = "Temp (range)",
  "o2_ltmin"                 = "O2 (min)",
  "phyc_mean"                = "Phyto (mean)",
  "phyc_range"               = "Phyto (range)",
  "sws_mean"                 = "Current (mean)",
  "sws_range"                = "Current (range)",
  "MS_bathy_5m"              = "Bathymetry",
  "MS_biogeo08_sss_mean_5m"  = "SSS (mean)",
  "MS_biogeo11_sss_range_5m" = "SSS (range)"
)

env_mat        <- as.matrix(env_data_11[, env_vars_11])
colnames(env_mat) <- clean_names[env_vars_11]
cor_mat        <- cor(env_mat, method = "pearson",
                       use = "complete.obs")

cat("Correlation matrix computed\n")
cat("High correlations (|r| > 0.7):\n")
cor_upper <- cor_mat
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA
high_cor  <- which(abs(cor_upper) > 0.7, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  for (i in 1:nrow(high_cor)) {
    cat(rownames(cor_upper)[high_cor[i,1]], "—",
        colnames(cor_upper)[high_cor[i,2]], ":",
        round(cor_upper[high_cor[i,1], high_cor[i,2]], 3), "\n")
  }
} else cat("No pairs with |r| > 0.7\n")

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_env_correlation_heatmap.pdf",
                               width = 8, height = 7)
  else png("RA_set1_env_correlation_heatmap.png",
           width = 8, height = 7, units = "in", res = 300)
  corrplot(cor_mat,
           method   = "color",
           type     = "upper",
           tl.col   = "black",
           tl.cex   = 0.8,
           tl.srt   = 45,
           addCoef.col = "black",
           number.cex  = 0.65,
           col      = colorRampPalette(c("#2166AC","white","#B2182B"))(200),
           diag     = FALSE,
           mar      = c(0, 0, 1, 0))
  dev.off()
}
cat("Environmental correlation heatmap saved\n")

# Save
save(env_data_11, cor_mat,
     file = "RA_set1_env_data.RData")
cat("Section 4 complete\n")

# ============================================================
# SECTION 5: Marine Distances and dbMEMs
# marineDistances (Assis et al. — https://github.com/jorgeassis/marineDistances)
# Shortest-path algorithm with infinite resistance of landmasses
# Resolution = 0.5 degrees; runtime ~8 hours (ran overnight, kept laptop open)
# Requires GSHHS_h_L1.shp and associated files in working directory
# Download: https://www.soest.hawaii.edu/pwessel/gshhg/
# PCoA correction applied to obtain Euclidean distances
# dbMEMs computed from Euclidean distances (adespatial v0.3.29)
# Only positive eigenvalues retained (n=2) — broad-scale
# spatial gradients representing positive autocorrelation
# dbMEMs are identical for set1 and set2 — derived from
# coordinates, not SNP data
# NOTE: marine distances already calculated — load if available
# ============================================================

library(adespatial)

# ── Option A: Load existing spatial data (recommended) ───────────────────────
load("RA_spatial_data.RData")
# Objects loaded: marine_dist, marine_euclid, dbmems

cat("Marine distances loaded:", dim(marine_dist), "\n")
cat("dbMEMs:", ncol(dbmems), "variables\n")
cat("Positive eigenvalues retained:", ncol(dbmems), "\n")

# Verify alignment with coords
stopifnot(all(rownames(marine_dist) == coords$ID))
stopifnot(all(rownames(as.data.frame(dbmems)) == coords$ID))
cat("Marine distances-coords alignment ✓\n")

# ── Option B: Re-run marineDistances (uncomment if needed) ───────────────────
# Requires GSHHS_h_L1.shp shapefile in working directory
#
# write.table(
#   data.frame(Name = coords$ID, Lon = coords$lon, Lat = coords$lat),
#   file      = "RA_individuals.txt",
#   sep       = " ",
#   row.names = FALSE,
#   col.names = FALSE,
#   quote     = FALSE
# )
#
# source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")
# contour(
#   global.polygon = "GSHHS_h_L1.shp",
#   file           = "RA_individuals.txt",
#   file.sep       = " ",
#   file.dec       = ".",
#   file.strucutre = 1,
#   file.header    = FALSE,
#   resolution     = 0.5,
#   buffer         = c(5, 5, 5, 5),
#   export.file    = TRUE
# )
#
# marine_dist <- as.matrix(
#   read.table("Contour _ Pairwise Marine Distances.txt",
#              header = TRUE, row.names = 1,
#              sep = " ", dec = ".")
# )
# cat("Marine distances computed:", dim(marine_dist), "\n")
#
# # Sanity checks — approximate expected distances
# cat("TZ-MZ (expect ~2000-3000 km):",
#     round(marine_dist["TZ2","MZ1"], 0), "\n")
# cat("MZ-AU (expect ~8000-10000 km):",
#     round(marine_dist["MZ1","AU1"], 0), "\n")
#
# # PCoA → Euclidean distances
# # Corrects for non-Euclidean nature of marine distances
# marine_pcoa   <- cmdscale(as.dist(marine_dist),
#                            k   = nrow(marine_dist) - 1,
#                            eig = TRUE)
# n_pos         <- sum(marine_pcoa$eig > 0)
# cat("Positive eigenvalues:", n_pos, "\n")
# marine_euclid <- marine_pcoa$points[, 1:n_pos]
# rownames(marine_euclid) <- coords$ID
# marine_euclid <- dist(marine_euclid)
# ***
# # dbMEMs from Euclidean distances
# dbmems <- dbmem(marine_euclid)
# cat("dbMEMs computed:", ncol(dbmems), "variables\n")
# cat("MEM eigenvalues (first 5):\n")
# print(attributes(dbmems)$values[1:5])
#
# save(marine_dist, marine_euclid, dbmems,
#      file = "RA_spatial_data.RData")
# cat("Spatial data saved\n")

# ── Prepare dbMEM CSV for HPC ─────────────────────────────────────────────────
# dbMEMs must be uploaded to HPC for ordistep and pRDA
dbmem_df <- data.frame(ID   = coords$ID,
                        MEM1 = as.numeric(dbmems[, 1]),
                        MEM2 = as.numeric(dbmems[, 2]))

write.csv(dbmem_df, "RA_dbmems.csv", row.names = FALSE)
cat("dbMEM CSV saved for HPC upload\n")

# Save
save(marine_dist, marine_euclid, dbmems, dbmem_df,
     file = "RA_spatial_data.RData")
cat("Section 5 complete\n")

# ============================================================
# SECTION 6: Prepare ordistep input for HPC
# Env-only ordistep run on HPC — 11 starting variables
# Bidirectional stepwise selection (Blanchet et al. 2008)
# 500 permutations — vegan ordistep function
# MEMs not ordistep-ped — orthogonal by construction
# (VIF = 1 confirmed), both retained
# After ordistep: VIF filtering locally (Section 7)
# Then pRDA prepared and submitted to HPC (Section 8)
#
# NOTE: ordistep already run on HPC — load results below
# HPC scripts: set1_RDA_ordistep11.R and set1_RDA_ordistep11.sh
# Results:    RA_set1_ordistep_results_11.rds
# ============================================================

# ── Prepare and save input files for HPC ─────────────────────────────────────
# These files are uploaded to HPC before running ordistep
# G_imp.rds and locus_names_full.rds already saved in Section 1

# Environmental data — 11 variables
env_vars_11 <- c("ph_mean","thetao_mean","thetao_range",
                  "o2_ltmin","phyc_mean","phyc_range",
                  "sws_mean","sws_range",
                  "MS_bathy_5m",
                  "MS_biogeo08_sss_mean_5m",
                  "MS_biogeo11_sss_range_5m")

# Must match popmap order — same order as G_imp rows
env_11 <- env_data_11[match(popmap$ID, env_data_11$ID),
                       c("ID", env_vars_11)]

# Hard verification
stopifnot(all(env_11$ID == popmap$ID))
stopifnot(all(env_11$ID == rownames(G_imp)))
stopifnot(sum(is.na(env_11)) == 0)
cat("Env-popmap alignment ✓\n")
cat("Env-G_imp alignment ✓\n")
cat("Missing values: 0 ✓\n")
cat("Dimensions:", dim(env_11), "\n")

write.csv(env_11, "RA_environmental_data_11vars.csv",
          row.names = FALSE)
write.csv(coords,  "RA_coordinates.csv",
          row.names = FALSE)

# Save G_imp as RDS for HPC
saveRDS(G_imp,           "G_imp.rds")
saveRDS(locus_names_full, "locus_names_full.rds")

cat("HPC input files saved:\n")
cat("  G_imp.rds\n")
cat("  locus_names_full.rds\n")
cat("  RA_environmental_data_11vars.csv\n")
cat("  RA_coordinates.csv\n")
cat("  RA_dbmems.csv\n")
cat("Upload these to the HPC working directory for the RDA step\n")

# ── HPC ordistep script ───────────────────────────────────────────────────────
# set1_RDA_ordistep11.R and set1_RDA_ordistep11.sh

# ── Load ordistep results ─────────────────────────────────────────────────────
# Download from HPC before running this section
ordistep_res <- readRDS("RA_set1_ordistep_results_11.rds")
selected_env <- ordistep_res$selected_env

cat("Ordistep selected variables:\n")
print(selected_env)
# Expected: "MS_biogeo11_sss_range_5m" "ph_mean" "thetao_mean" "MS_biogeo08_sss_mean_5m" "thetao_range" "phyc_range" "sws_mean

cat("Section 6 complete\n")

# ============================================================
# SECTION 7: VIF filtering — local
# Step 1: env-only VIF — threshold <10
#         remove variables iteratively in order of highest VIF
# Step 2: SSS range excluded despite ordistep selection
#         r=0.93 with MEM1, VIF>15 in combined env+spatial set
# Step 3: combined env + MEMs VIF — threshold <15
#         MEMs are orthogonal (VIF=1) by mathematical construction
# Final set1 variables: ph_mean, thetao_range, phyc_range, sws_mean
# ============================================================

library(usdm)

# ── Step 1: env-only VIF on ordistep-selected variables ──────────────────────
env_selected <- env_data_11[, selected_env]
rownames(env_selected) <- env_data_11$ID

cat("=== Env-only VIF (threshold < 10) ===\n")
vif_env <- vif(as.data.frame(env_selected))
print(vif_env)

# Iterative removal if any VIF > 10
# vifstep removes highest VIF variable until all < threshold
vif_step_env <- vifstep(as.data.frame(env_selected), th = 10)
cat("Variables retained after env-only VIF < 10:\n")
print(vif_step_env@results$Variables)
env_final <- selected_env[selected_env %in% vif_step_env@results$Variables]
cat("Final env variables after VIF step 1:", env_final, "\n")

# MS_biogeo11_sss_range_5m ph_mean thetao_range phyc_range sws_mean

# ── Step 2: combined env + MEMs VIF ──────────────────────────────────────────
env_mem_combined <- data.frame(
  env_data_11[match(popmap$ID, env_data_11$ID), env_final],
  MEM1 = as.numeric(dbmems[, 1]),
  MEM2 = as.numeric(dbmems[, 2])
)

cat("=== Combined env + spatial VIF (threshold < 15) ===\n")
vif_combined <- vif(as.data.frame(env_mem_combined))
print(vif_combined)

# Check SSS range - MEM1 correlation
r_sss_mem1 <- cor(
  env_data_11$MS_biogeo11_sss_range_5m[match(popmap$ID,
                                               env_data_11$ID)],
  as.numeric(dbmems[, 1]),
  use = "complete.obs")
cat("SSS range correlation with MEM1:", round(r_sss_mem1, 3), "\n")
# Expected: r = 0.93

# NOTE: vifstep() would automatically remove MEM1 (highest VIF)
# We override this decision — MEMs are retained by design as
# spatial conditioning variables in the partial RDA framework
# SSS range is removed instead: r=0.93 with MEM1 confirms it
# captures the same broad-scale geographic-salinity gradient
# rather than an independent environmental signal

env_final_4 <- env_final[env_final != "MS_biogeo11_sss_range_5m"]
cat("Final env variables after manual SSS range removal:",
    env_final_4, "\n")
# Expected: ph_mean, thetao_range, phyc_range, sws_mean

# Re-run VIF on final 4 env + MEMs to confirm all < 15
env_mem_final <- data.frame(
  env_data_11[match(popmap$ID, env_data_11$ID), env_final_4],
  MEM1 = as.numeric(dbmems[, 1]),
  MEM2 = as.numeric(dbmems[, 2])
)

vif_final <- vif(as.data.frame(env_mem_final))
cat("\n=== Final VIF — 4 env + 2 MEMs ===\n")
print(vif_final)
cat("VIF range:", round(min(vif_final$VIF), 3),
    "-", round(max(vif_final$VIF), 3), "\n")

# ── Save pRDA input — FINAL 4 variables only ─────────────────────────────────
# SSS range excluded before pRDA — captures spatial not environmental signal
# pRDA uses: ph_mean, thetao_range, phyc_range, sws_mean

env_final_df <- env_data_11[match(popmap$ID, env_data_11$ID),
                              c("ID", env_final_4)]
stopifnot(all(env_final_df$ID == popmap$ID))
stopifnot(sum(is.na(env_final_df)) == 0)
stopifnot(ncol(env_final_df) == 5)  # ID + 4 variables

write.csv(env_final_df, "RA_env_final_set1_11vars.csv",
          row.names = FALSE)
cat("pRDA input saved: RA_env_final_set1_11vars.csv\n")
cat("Contains 4 final variables:", env_final_4, "\n")
cat("Upload to HPC before submitting pRDA\n")

# ── Confirm pRDA input is correct ────────────────────────────────────────────
env_check <- read.csv("RA_env_final_set1_11vars.csv")
cat("Columns in pRDA input:", colnames(env_check), "\n")
cat("Number of env variables:", ncol(env_check) - 1, "\n")
stopifnot(ncol(env_check) - 1 == 4)
cat("pRDA input verified — 4 variables ✓\n")

# ── Save VIF results ──────────────────────────────────────────────────────────
save(vif_combined, vif_final, env_final_4, r_sss_mem1,
     file = "RA_set1_VIF_results.RData")
cat("Section 7 complete\n")
cat("Final variables for pRDA:", env_final_4, "\n")
cat("VIF range:", round(min(vif_final$VIF), 3),
    "-", round(max(vif_final$VIF), 3), "\n")

# NOTE: filename says 11vars but contains FINAL 4 variables only

# ============================================================
# SECTION 8: pRDA — HPC
# Partial RDA conditioning env variables on dbMEMs (space)
# Isolates pure environmental genomic signal from neutral
# spatial population structure
# 999-permutation ANOVA for model and axis significance
# Variance partitioning: env alone, space alone, shared, residual
#
# NOTE: pRDA already run on HPC — load results below
# HPC scripts: set1_partialRDA.R and set1_partialRDA.sh
# Results:    RA_set1_partialRDA_results_11vars.rds
# ============================================================

# ── HPC pRDA script ─────────────────────────
# set1_partialRDA.R and set1_partialRDA.sh

# ── Load pRDA results ─────────────────────────────────────────────────────────
# Download from HPC before running this section
# Note: local file confirmed — RA_set1_partialRDA_results_11.rds
library(vegan)

pRDA_results <- readRDS("RA_set1_partialRDA_results_11.rds")

rda_partial        <- pRDA_results$rda_partial
anova_partial      <- pRDA_results$anova_partial
anova_partial_axes <- pRDA_results$anova_partial_axes
vp                 <- pRDA_results$vp
env_only_final     <- pRDA_results$env_only_final
mem_df             <- pRDA_results$mem_df

cat("pRDA results loaded\n")
cat("Partial RDA summary:\n")
cat("Adj R² (pure env):",
    round(RsquareAdj(rda_partial)$adj.r.squared, 4), "\n")
cat("Model significance:\n"); print(anova_partial)
cat("Axis significance:\n");  print(anova_partial_axes)
cat("Variance partitioning:\n"); print(vp)

cat("=== Significance tests from HPC pRDA ===\n")
cat("[a] Pure env (partial RDA):\n")
print(pRDA_results$anova_partial)

cat("\n[a+c] Env-only:\n")
print(pRDA_results$anova_env)

cat("\n[b+c] Space-only:\n")
print(pRDA_results$anova_sp)

cat("\n[a+b+c] Full model:\n")
print(pRDA_results$anova_full)

# Report key results
cat("\n── Key results for methods/results text ──\n")
cat("Pure env [a] adj R²:",
    round(vp$part$indfract[1,3], 4), "\n")
cat("Pure space [b] adj R²:",
    round(vp$part$indfract[2,3], 4), "\n")
cat("Shared [c] adj R²:",
    round(vp$part$indfract[3,3], 4), "\n")
cat("Residuals [d]:",
    round(vp$part$indfract[4,3], 4), "\n")

save(rda_partial, anova_partial, anova_partial_axes,
     vp, env_only_final, mem_df,
     file = "RA_set1_pRDA_results_local.RData")
cat("Section 8 complete\n")

# ============================================================
# SECTION 9: RDA outlier detection
# Candidate loci under environmental selection identified
# using 3 SD threshold on each significant RDA axis
# Applied to RDA1, RDA2, RDA3 only (RDA4 p=0.359 excluded)
# Result: 846 outliers (3.1% of 27 562 SNPs)
# ============================================================

library(vegan)

# ── Extract locus scores on constrained axes ──────────────────────────────────
# Species scores in RDA = locus loadings on each axis
# scaling=1: distances among sites are meaningful
rda_scores <- scores(rda_partial,
                     display  = "species",
                     scaling  = 1,
                     choices  = 1:3)  # RDA1-3 only (significant)

cat("Locus scores dimensions:", dim(rda_scores), "\n")
cat("Axes:", colnames(rda_scores), "\n")

# ── 3SD threshold per axis ────────────────────────────────────────────────────
# Loci exceeding ±3 SD from mean on any significant axis
# are considered candidates under environmental selection

outliers_per_axis <- list()
for (ax in colnames(rda_scores)) {
  scores_ax  <- rda_scores[, ax]
  mean_ax    <- mean(scores_ax)
  sd_ax      <- sd(scores_ax)
  thresh_pos <-  mean_ax + RDA_SD * sd_ax
  thresh_neg <-  mean_ax - RDA_SD * sd_ax

  out_idx <- which(scores_ax > thresh_pos | scores_ax < thresh_neg)
  outliers_per_axis[[ax]] <- locus_names_full[out_idx]

  cat(ax, "— threshold: ±", round(RDA_SD * sd_ax, 4),
      "| outliers:", length(out_idx), "\n")
}

# ── Combine outliers across all significant axes ──────────────────────────────
# Union — locus flagged on ANY significant axis is a candidate
outliers_rda       <- unique(unlist(outliers_per_axis))
outliers_rda_idx   <- which(locus_names_full %in% outliers_rda)

cat("\nRDA outliers per axis:\n")
for (ax in names(outliers_per_axis)) {
  cat(" ", ax, ":", length(outliers_per_axis[[ax]]), "\n")
}
cat("Total unique RDA outliers:", length(outliers_rda), "\n")
cat("Percentage of SNPs:",
    round(length(outliers_rda)/N_SNP*100, 1), "%\n")
# Expected: 846 outliers (3.1%)

# ============================================================
# SNP cluster plots — RDA1 vs RDA2
# Unlabelled (no arrows) — version for publication figures
# Unlabelled with arrows — Canva input (add labels in Canva)
# Labelled with arrows — for inspection
# Candidate loci: orange (#F5793A) — exceed ±3SD on RDA1-3
# Reference loci: blue (#85C0F9)
# Environmental arrows scaled to 40% of plot range
# NOTE: unlabelled-with-arrows version is Canva input
# ============================================================

library(vegan)

# ── Axis variance percentages ─────────────────────────────────────────────────
eig_vals <- summary(eigenvals(rda_partial,
                               model = "constrained"))
ax1_var  <- round(eig_vals["Proportion Explained", "RDA1"] * 100, 1)
ax2_var  <- round(eig_vals["Proportion Explained", "RDA2"] * 100, 1)
cat("RDA1:", ax1_var, "% | RDA2:", ax2_var, "%\n")
# Expected: RDA1=47.2%, RDA2=26.7%

# ── Locus scores ──────────────────────────────────────────────────────────────
rda_all  <- scores(rda_partial, display = "species",
                   scaling = 1, choices = 1:3)
snp_rda1 <- rda_all[, "RDA1"]
snp_rda2 <- rda_all[, "RDA2"]

cat("Locus scores dimensions:", dim(rda_all), "\n")
# Expected: 27562 x 3

# ── Candidate/reference classification ───────────────────────────────────────
is_candidate <- locus_names_full %in% outliers_rda
snp_col <- ifelse(is_candidate,
                  adjustcolor("#F5793A", alpha.f = 0.7),
                  adjustcolor("#85C0F9", alpha.f = 0.3))
snp_cex <- ifelse(is_candidate, 0.5, 0.3)

cat("Candidate loci:", sum(is_candidate), "\n")
cat("Reference loci:", sum(!is_candidate), "\n")

# ── Plot limits ───────────────────────────────────────────────────────────────
x_snp <- max(abs(snp_rda1)) * 1.2
y_snp <- max(abs(snp_rda2)) * 1.2

# ── Environmental arrows ──────────────────────────────────────────────────────
bp_scores <- scores(rda_partial, display = "bp", scaling = 1)
cat("Environmental variable scores:\n")
print(round(bp_scores, 4))

# Scale arrows to 40% of plot range
# Ensures arrows visible but do not dominate the SNP cloud
arrow_scale <- 0.6 * min(x_snp, y_snp) /
               max(sqrt(bp_scores[,1]^2 + bp_scores[,2]^2))
cat("Arrow scale factor:", round(arrow_scale, 4), "\n")

env_labels <- c(
  "ph_mean"      = "pH",
  "thetao_range" = "Temp range",
  "phyc_range"   = "Phyto range",
  "sws_mean"     = "Current mean"
)

# ── Legend ────────────────────────────────────────────────────────────────────
add_legend <- function() {
  legend("topright",
         legend = c(paste0("Candidate loci (n=", sum(is_candidate), ")"),
                    paste0("Reference loci (n=", sum(!is_candidate), ")")),
         col    = c(adjustcolor("#F5793A", alpha.f = 0.7),
                    adjustcolor("#85C0F9", alpha.f = 0.3)),
         pch    = 19,
         pt.cex = c(0.8, 0.6),
         bty    = "n",
         cex    = 0.85)
}

# ── Base plot function ────────────────────────────────────────────────────────
draw_base <- function() {
  plot(snp_rda1, snp_rda2,
       pch     = 19,
       col     = snp_col,
       cex     = snp_cex,
       xlim    = c(-x_snp, x_snp),
       ylim    = c(-y_snp, y_snp),
       xlab    = paste0("RDA1 (", ax1_var, "%)"),
       ylab    = paste0("RDA2 (", ax2_var, "%)"),
       cex.lab = 1.2,
       cex.axis = 1.0,
       las     = 1,
       bty     = "o")
  abline(h = 0, v = 0, lty = 2, col = "grey70", lwd = 0.5)
}

# ── Arrow drawing function ────────────────────────────────────────────────────
draw_arrows <- function(labelled = TRUE) {
  for (i in 1:nrow(bp_scores)) {
    arrows(0, 0,
           bp_scores[i, 1] * arrow_scale,
           bp_scores[i, 2] * arrow_scale,
           length = 0.12, lwd = 2.0,
           col = "black")
    if (labelled) {
      text(bp_scores[i, 1] * arrow_scale * 1.2,
           bp_scores[i, 2] * arrow_scale * 1.2,
           labels = env_labels[rownames(bp_scores)[i]],
           cex    = 0.9,
           font   = 3,
           col    = "black")
    }
  }
}

# ── 1. Unlabelled — no arrows — publication figure ─────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_SNP_cluster_unlabelled.pdf",
                               width = 7, height = 7)
  else png("RA_set1_SNP_cluster_unlabelled.png",
           width = 7, height = 7, units = "in", res = 300)
  draw_base()
  add_legend()
  dev.off()
}
cat("Unlabelled SNP cluster plot saved\n")

# ── 2. Unlabelled arrows — Canva input (add labels in Canva) ─────────────────
# NOTE: labels added manually in Canva — do not use standalone
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_SNP_cluster_arrows_unlabelled.pdf",
                               width = 7, height = 7)
  else png("RA_set1_SNP_cluster_arrows_unlabelled.png",
           width = 7, height = 7, units = "in", res = 300)
  draw_base()
  draw_arrows(labelled = FALSE)
  add_legend()
  dev.off()
}
cat("Unlabelled arrows SNP cluster plot saved — Canva input\n")

# ── 3. Labelled arrows — for inspection ─────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_SNP_cluster_labelled.pdf",
                               width = 7, height = 7)
  else png("RA_set1_SNP_cluster_labelled.png",
           width = 7, height = 7, units = "in", res = 300)
  draw_base()
  draw_arrows(labelled = TRUE)
  add_legend()
  dev.off()
}
cat("Labelled SNP cluster plot saved\n")

cat("All SNP cluster plots complete\n")
cat("Files produced:\n")
cat("  RA_set1_SNP_cluster_unlabelled.pdf/png — publication figure\n")
cat("  RA_set1_SNP_cluster_arrows_unlabelled.pdf/png — Canva input\n")
cat("  RA_set1_SNP_cluster_labelled.pdf/png — inspection\n")

# ── Candidate-environment correlations ────────────────────────────────────────
# Pearson correlation between RDA candidate locus genotype
# frequencies and each retained environmental variable
# n=846 RDA environmental candidates only
# Significance: p<0.05 (two-tailed Pearson, df=69)
# Multiple testing: BH FDR correction across all 846×4=3,384 tests
# Strong effect: |r|>0.5 (Cohen 1988 large effect, r²>0.25)
# Variables ranked by % candidates showing significant associations

R_STRONG <- 0.5  # large effect size threshold (Cohen 1988)

env_vars_final <- c("ph_mean","thetao_range","phyc_range","sws_mean")
env_display    <- c(
  "ph_mean"      = "pH",
  "thetao_range" = "Temperature range",
  "phyc_range"   = "Phytoplankton range",
  "sws_mean"     = "Current velocity mean"
)

# Verify alignment
stopifnot(all(colnames(G_imp) == locus_names_full))
cat("G_imp column alignment ✓\n")

G_candidates <- G_imp[, outliers_rda]
cat("Candidate genotype matrix:", dim(G_candidates), "\n")

env_mat <- as.matrix(
  env_data_11[match(popmap$ID, env_data_11$ID), env_vars_final])
rownames(env_mat) <- popmap$ID
stopifnot(all(rownames(env_mat) == rownames(G_imp)))
cat("Env matrix aligned ✓\n")

# ── Calculate Pearson r AND p-value for each candidate x env var ──────────────
n_cand <- length(outliers_rda)
n_env  <- length(env_vars_final)

cor_r <- matrix(NA, nrow = n_cand, ncol = n_env,
                 dimnames = list(outliers_rda, env_vars_final))
cor_p <- matrix(NA, nrow = n_cand, ncol = n_env,
                 dimnames = list(outliers_rda, env_vars_final))

cat("Calculating correlations...\n")
for (j in seq_along(env_vars_final)) {
  for (i in seq_along(outliers_rda)) {
    ct           <- cor.test(G_candidates[, i],
                              env_mat[, j],
                              method = "pearson")
    cor_r[i, j] <- ct$estimate
    cor_p[i, j] <- ct$p.value
  }
}
cat("Correlations complete\n")

# ── BH FDR correction across all tests ───────────────────────────────────────
# 846 candidates × 4 variables = 3,384 tests
# BH correction applied across all tests simultaneously
p_vec     <- as.vector(cor_p)           # 3,384 p-values
p_adj_vec <- p.adjust(p_vec, method = "BH")
cor_p_adj <- matrix(p_adj_vec,
                     nrow = n_cand, ncol = n_env,
                     dimnames = list(outliers_rda, env_vars_final))

cat("BH FDR correction applied across",
    n_cand * n_env, "tests\n")

# ── Per-locus summary ─────────────────────────────────────────────────────────
max_r    <- apply(abs(cor_r), 1, max)
best_env <- env_vars_final[apply(abs(cor_r), 1, which.max)]

cor_summary <- data.frame(
  locus    = outliers_rda,
  best_env = best_env,
  max_r    = round(max_r, 3),
  stringsAsFactors = FALSE
)

# Sanity check — best_env should correspond to highest |r|
cat("Sanity check — best_env assignments:\n")
check <- sapply(seq_along(outliers_rda), function(i) {
  which.max(abs(cor_r[i,])) ==
    which(env_vars_final == cor_summary$best_env[i])
})
cat("Correct assignments:", sum(check), "/", length(check), "\n")
stopifnot(all(check))
cat("All best_env assignments correct ✓\n")

# ── Summary table per environmental variable ──────────────────────────────────
# n significant = FDR-corrected p<0.05
# % significant = % of ALL 846 candidates (not just assigned)
# mean |r| = across all 846 loci (not just significant)
# max |r| = maximum for any single locus
# n strong = |r|>0.5 across all 846 loci

cor_table <- do.call(rbind, lapply(env_vars_final, function(v) {
  r_vals   <- abs(cor_r[, v])
  p_adj    <- cor_p_adj[, v]
  n_sig    <- sum(p_adj < 0.05, na.rm = TRUE)
  n_strong <- sum(r_vals >= R_STRONG, na.rm = TRUE)
  data.frame(
    Variable        = env_display[v],
    n_significant   = n_sig,
    pct_significant = round(n_sig / n_cand * 100, 1),
    mean_abs_r      = round(mean(r_vals, na.rm = TRUE), 3),
    max_abs_r       = round(max(r_vals,  na.rm = TRUE), 3),
    n_strong        = n_strong,
    stringsAsFactors = FALSE
  )
}))

# Rank by % significant (descending)
cor_table <- cor_table[order(cor_table$pct_significant,
                               decreasing = TRUE), ]
rownames(cor_table) <- NULL

cat("\n── Candidate-environment correlation table ──\n")
cat("Ranked by % significant (FDR-corrected p<0.05)\n")
print(cor_table)

cat("\nCandidates assigned per variable (best |r|):\n")
print(table(cor_summary$best_env))
cat("Mean |r| across all candidates:",
    round(mean(max_r), 3), "\n")
cat("Max |r| across all candidates:",
    round(max(max_r), 3), "\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
write.csv(cor_summary,
          "RA_set1_candidate_correlations_per_locus.csv",
          row.names = FALSE)
write.csv(cor_table,
          "RA_set1_candidate_correlations_summary.csv",
          row.names = FALSE)

# ── Save Section 9 objects ────────────────────────────────────────────────────
save(rda_scores, outliers_per_axis, outliers_rda,
     outliers_rda_idx, is_candidate,
     snp_rda1, snp_rda2,
     cor_r, cor_p, cor_p_adj, cor_results,
     cor_summary, cor_table,
     R_STRONG,
     file = "RA_set1_RDA_outlier_results.RData")
cat("Section 9 complete\n")
cat("RDA outliers:", length(outliers_rda),
    "(", round(length(outliers_rda)/N_SNP*100, 1), "%)\n")

# ============================================================
# SECTION 10: LFMM
# Latent Factor Mixed Models — univariate GEA
# lfmm2 least-squares estimation (Caye et al. 2019)
# LEA v3.16.0 (Frichot & François 2015)
# Input: 11 environmental variables → PCA → broken stick
# K=4 latent factors — consistent with pcadapt optimal K
# GIF=2.033 — justified by Q-Q plot inspection:
#   tail-only departure indicates true outliers not inflation
#   genome-wide inflation would show linear Q-Q departure
# FDR threshold: q<0.01 (stricter than pcadapt q<0.05)
#   rationale: LFMM univariate tests produce more candidates
#   per locus — stricter threshold controls false positives
# Result: 1 762 outliers (6.4% of 27,562 SNPs)
# ============================================================

library(LEA)
library(qvalue)
library(vegan)  # for bstick()

# ── Option A: Load existing results (recommended) ─────────────────────────────
  load("RA_set1_LFMM_results_11vars.RData")
  cat("LFMM results loaded from saved file\n")
  cat("Objects loaded:", ls(), "\n")
  cat("LFMM outliers:", length(outliers_lfmm_names), "\n")
  cat("GIF:", round(gif, 3), "\n")
  # Expected: 1 762 outliers, GIF=2.033

# ── Option B: Rerun LFMM (if results file missing/corrupted) ─────────────────

# ── 10.1 Prepare environmental PCA input ─────────────────────────────────────
# All 11 starting variables used for PCA
# PCA handles collinearity — no VIF filtering needed here
# This is methodologically distinct from RDA which uses VIF+ordistep
env_vars_11 <- c("ph_mean","thetao_mean","thetao_range",
                  "o2_ltmin","phyc_mean","phyc_range",
                  "sws_mean","sws_range",
                  "MS_bathy_5m",
                  "MS_biogeo08_sss_mean_5m",
                  "MS_biogeo11_sss_range_5m")

env_all_lfmm <- env_data_11[match(popmap$ID, env_data_11$ID),
                              env_vars_11]
rownames(env_all_lfmm) <- popmap$ID

# Verify alignment
stopifnot(all(rownames(env_all_lfmm) == rownames(G_imp)))
stopifnot(sum(is.na(env_all_lfmm)) == 0)
cat("LFMM env data aligned ✓\n")
cat("Variables:", ncol(env_all_lfmm), "\n")

# ── 10.2 Environmental PCA ────────────────────────────────────────────────────
env_pca <- prcomp(env_all_lfmm, scale = TRUE, center = TRUE)

# Broken stick criterion — retain PCs with eigenvalues
# exceeding those expected from random distribution
eigenvalues <- env_pca$sdev^2 / sum(env_pca$sdev^2)
bstick_vals <- bstick(length(env_pca$sdev))
n_pcs       <- sum(eigenvalues > bstick_vals)

cat("PCs retained (broken stick):", n_pcs, "\n")
cat("Variance explained by retained PCs:",
    round(sum(eigenvalues[1:n_pcs]) * 100, 1), "%\n")
# Expected: 3 PCs explaining 81.4%

# Broken stick plot
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_LFMM_brokenstick.pdf",
                               width = 6, height = 5)
  else png("RA_set1_LFMM_brokenstick.png",
           width = 6, height = 5, units = "in", res = 300)
  plot(eigenvalues, type = "b", pch = 19,
       xlab     = "Principal component",
       ylab     = "Proportion of variance",
       ylim     = c(0, max(c(eigenvalues, bstick_vals)) + 0.05),
       cex.lab  = 1.2, cex.axis = 1.0,
       las = 1, bty = "l")
  lines(bstick_vals, col = "red", lty = 2)
  abline(v = n_pcs + 0.5, col = "grey50", lty = 3)
  legend("topright",
         legend = c("Observed eigenvalue",
                    "Broken stick criterion",
                    paste0("Retained: ", n_pcs, " PCs")),
         col    = c("black","red","grey50"),
         lty    = c(1, 2, 3),
         pch    = c(19, NA, NA),
         bty    = "n", cex = 0.9)
  dev.off()
}
cat("Broken stick plot saved\n")

# Extract retained PCs
env_pcs <- env_pca$x[, 1:n_pcs, drop = FALSE]
cat("env_pcs dimensions:", dim(env_pcs), "\n")
# Expected: 71 x 3

# ── 10.3 Run LFMM2 ───────────────────────────────────────────────────────────
# K=4 latent factors — consistent with pcadapt optimal K
# lfmm2 uses least-squares estimation (Caye et al. 2019)
# More computationally efficient than ridge regression
# without loss of power or precision
cat("Running LFMM2 (K=4)...\n")
cat("Start time:", format(Sys.time()), "\n")

mod_lfmm <- lfmm2(
  input = G_imp,
  env   = env_pcs,
  K     = K_LFMM  # K=4 from parameters section
)

# ── 10.4 P-values and GIF ────────────────────────────────────────────────────
pv <- lfmm2.test(
  object = mod_lfmm,
  input  = G_imp,
  env    = env_pcs,
  full   = TRUE  # combined test across all env PCs
)

gif <- pv$gif
cat("GIF:", round(gif, 3), "\n")
# Expected: 2.033
# Justification: Q-Q plot shows tail-only departure
# indicating true biological outliers not genome-wide inflation
# See Q-Q plot below for confirmation

cat("End time:", format(Sys.time()), "\n")

# ── 10.5 Q-Q plot — GIF justification ────────────────────────────────────────
# Critical diagnostic — must show tail-only departure
# Linear departure throughout = inflation requiring correction
# Tail-only departure = true outliers, GIF inflation acceptable
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_LFMM_qqplot.pdf",
                               width = 6, height = 6)
  else png("RA_set1_LFMM_qqplot.png",
           width = 6, height = 6, units = "in", res = 300)
  pvals_sorted <- sort(pv$pvalues)
  n            <- length(pvals_sorted)
  expected     <- (1:n) / n
  plot(-log10(expected), -log10(pvals_sorted),
       pch     = 19, cex = 0.3,
       col     = adjustcolor("grey40", alpha.f = 0.5),
       xlab    = expression("Expected" ~ -log[10](italic(p))),
       ylab    = expression("Observed" ~ -log[10](italic(p))),
       cex.lab = 1.2, cex.axis = 1.0,
       las = 1, bty = "l")
  abline(0, 1, col = "red", lty = 2, lwd = 1.5)
  mtext(paste0("GIF = ", round(gif, 3)),
        side = 3, adj = 1, cex = 0.9, font = 3)
  dev.off()
}
cat("Q-Q plot saved — inspect for tail-only departure\n")

# ── 10.6 P-value histogram ────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_LFMM_pval_hist.pdf",
                               width = 7, height = 5)
  else png("RA_set1_LFMM_pval_hist.png",
           width = 7, height = 5, units = "in", res = 300)
  hist(pv$pvalues,
       breaks   = 50,
       col      = adjustcolor("#85C0F9", alpha.f = 0.8),
       border   = "white",
       xlab     = expression(italic(p)~value),
       ylab     = "Frequency",
       main     = "",
       cex.lab  = 1.2, cex.axis = 1.0,
       las = 1, bty = "l")
  mtext(paste0("GIF = ", round(gif, 3)),
        side = 3, adj = 1, cex = 0.9, font = 3)
  dev.off()
}
cat("P-value histogram saved\n")

# ── 10.7 Outlier detection — FDR q<0.01 ──────────────────────────────────────
# Stricter threshold than pcadapt (q<0.05)
# Rationale: univariate LFMM tests each locus independently
# producing more candidates — stricter FDR controls false positives
qv              <- qvalue(pv$pvalues)$qvalues
outliers_lfmm   <- which(qv < FDR_LFMM)  # FDR_LFMM = 0.01
outliers_lfmm_names <- locus_names_full[outliers_lfmm]

cat("LFMM outliers (q <", FDR_LFMM, "):",
    length(outliers_lfmm_names), "\n")
cat("Percentage:",
    round(length(outliers_lfmm_names)/N_SNP*100, 1), "%\n")
# Expected: 1,762 outliers (6.4%)

# FDR threshold for Manhattan plot
fdr_threshold_lfmm <- max(pv$pvalues[qv < FDR_LFMM],
                           na.rm = TRUE)
fdr_line_lfmm      <- -log10(fdr_threshold_lfmm)
cat("FDR threshold (-log10 p):", round(fdr_line_lfmm, 3), "\n")

# ── 10.8 Manhattan plot ───────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_LFMM_manhattan.pdf",
                               width = 10, height = 5)
  else png("RA_set1_LFMM_manhattan.png",
           width = 10, height = 5, units = "in", res = 300)
  plot(seq_along(pv$pvalues),
       -log10(pv$pvalues),
       pch      = 19,
       cex      = 0.3,
       col      = adjustcolor("grey40", alpha.f = 0.6),
       xlab     = "SNP index",
       ylab     = expression(-log[10](italic(p)~value)),
       cex.lab  = 1.2, cex.axis = 1.0,
       las = 1, bty = "l")
  abline(h   = fdr_line_lfmm,
         col = "#87c1f7",
         lty = 2, lwd = 1.5)
  legend("topright",
         legend = paste0("FDR (q < ", FDR_LFMM, ")"),
         col    = "#87c1f7",
         lty    = 2, lwd = 1.5,
         bty    = "n", cex = 0.9)
  dev.off()
}
cat("Manhattan plot saved\n")

# ── 10.9 Save ─────────────────────────────────────────────────────────────────
save(mod_lfmm, pv, gif, qv, outliers_lfmm, outliers_lfmm_names,
     env_pca, env_pcs, n_pcs, eigenvalues, bstick_vals,
     env_all_lfmm,
     file = "RA_set1_LFMM_results_11vars.RData")
cat("LFMM results saved\n")

cat("Section 10 complete\n")
cat("LFMM outliers:", length(outliers_lfmm_names),
    "(", round(length(outliers_lfmm_names)/N_SNP*100,1), "%)\n")
cat("GIF:", round(gif, 3), "\n")

# ============================================================
# SECTION 11: Four-dataset compilation
# Outlier framework using three methods:
# RDA (env GEA), pcadapt (structure-based), LFMM (univariate GEA)
# OutFLANK excluded — 0 outliers due to strong population
# structure inflating FST null distribution (Section 2)
# High-confidence: detected by ≥2 methods
# ============================================================

# ── Verify all three outlier sets are available ───────────────────────────────
stopifnot(exists("outliers_rda"))
stopifnot(exists("outliers_pcadapt_names"))
stopifnot(exists("outliers_lfmm_names"))

cat("RDA outliers:", length(outliers_rda), "\n")
cat("pcadapt outliers:", length(outliers_pcadapt_names), "\n")
cat("LFMM outliers:", length(outliers_lfmm_names), "\n")

# ── Four datasets ─────────────────────────────────────────────────────────────
# 1. All putative outliers — union of all three methods
outliers_any <- unique(c(outliers_rda,
                          outliers_pcadapt_names,
                          outliers_lfmm_names))

# 2. Environmental candidates — RDA only
# These loci respond specifically to environmental gradients
# after conditioning on spatial effects
outliers_enviro <- outliers_rda

# 3. High-confidence — detected by ≥2 methods
# Cross-validation across complementary statistical frameworks
# minimises false positives (Forester et al. 2018)
outliers_highconf <- unique(c(
  intersect(outliers_rda,          outliers_pcadapt_names),
  intersect(outliers_rda,          outliers_lfmm_names),
  intersect(outliers_pcadapt_names, outliers_lfmm_names)
))

# 4. Neutral — not detected by any method
outliers_neutral <- locus_names_full[
  !locus_names_full %in% outliers_any]

# ── Verify counts match confirmed values ──────────────────────────────────────
cat("\n── Four-dataset summary ──\n")
cat("Neutral:          ", length(outliers_neutral),
    "(", round(length(outliers_neutral)/N_SNP*100,1), "%)\n")
cat("All putative:     ", length(outliers_any),
    "(", round(length(outliers_any)/N_SNP*100,1), "%)\n")
cat("Environmental:    ", length(outliers_enviro),
    "(", round(length(outliers_enviro)/N_SNP*100,1), "%)\n")
cat("High confidence:  ", length(outliers_highconf),
    "(", round(length(outliers_highconf)/N_SNP*100,1), "%)\n")

# Hard checks — must match confirmed values
stopifnot(length(outliers_neutral)  == 24415)
stopifnot(length(outliers_any)      == 3147)
stopifnot(length(outliers_enviro)   == 846)
stopifnot(length(outliers_highconf) == 720)
cat("All SNP counts verified ✓\n")

# ── Method overlap summary ────────────────────────────────────────────────────
cat("\n── Method overlap ──\n")
cat("RDA ∩ pcadapt:",
    length(intersect(outliers_rda, outliers_pcadapt_names)), "\n")
cat("RDA ∩ LFMM:",
    length(intersect(outliers_rda, outliers_lfmm_names)), "\n")
cat("pcadapt ∩ LFMM:",
    length(intersect(outliers_pcadapt_names, outliers_lfmm_names)), "\n")
cat("All three:",
    length(intersect(intersect(outliers_rda,
                                outliers_pcadapt_names),
                     outliers_lfmm_names)), "\n")

# ── Subset genotype matrices ──────────────────────────────────────────────────
# Required for all downstream population genomics analyses
G_neutral  <- G_imp[, outliers_neutral]
G_any      <- G_imp[, outliers_any]
G_enviro   <- G_imp[, outliers_enviro]
G_highconf <- G_imp[, outliers_highconf]

cat("\nGenotype matrix dimensions:\n")
cat("G_neutral: ",  dim(G_neutral),  "\n")
cat("G_any:     ",  dim(G_any),      "\n")
cat("G_enviro:  ",  dim(G_enviro),   "\n")
cat("G_highconf:",  dim(G_highconf), "\n")

# ── Save ──────────────────────────────────────────────────────────────────────
save(outliers_neutral, outliers_any, outliers_enviro, outliers_highconf,
     G_neutral, G_any, G_enviro, G_highconf,
     outliers_rda, outliers_pcadapt_names, outliers_lfmm_names,
     file = "RA_set1_four_datasets_v2.RData")
cat("Section 11 complete — four datasets saved\n")

# ============================================================
# SECTION 12: Venn diagram
# Visualises overlap between three outlier detection methods
# OutFLANK excluded — documented in Section 2 and figure caption
# ggVennDiagram v1.5.7 (Gao et al. 2024)
# ============================================================

library(ggVennDiagram)
library(ggplot2)


venn_list <- list(
  RDA     = outliers_rda,
  pcadapt = outliers_pcadapt_names,
  LFMM    = outliers_lfmm_names
)

cat("Venn input sizes:\n")
cat("RDA:", length(venn_list$RDA), "\n")
cat("pcadapt:", length(venn_list$pcadapt), "\n")
cat("LFMM:", length(venn_list$LFMM), "\n")

p_venn <- ggVennDiagram(
  venn_list,
  label_alpha    = 0,
  set_size       = 4.5,
  label_size     = 3.5
) +
  scale_fill_gradientn(
    colours = c("white","#C6DBEF","#6BAED6")) +
  scale_color_manual(
    values = c("#333333","#333333","#333333")) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill  = "white",
                                    colour = NA))

for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_venn_labelled.pdf",
                               width = 6, height = 6)
  else png("RA_set1_venn_labelled.png",
           width = 6, height = 6,
           units = "in", res = 300)
  print(p_venn)
  dev.off()
}
cat("Venn diagram saved\n")
cat("Section 12 complete\n")

# ============================================================
# SECTION 13: RDA triplots
# RDA1 vs RDA2 (Panel A) and RDA1 vs RDA3 (Panel B)
# Individual points coloured by population
# Environmental arrows labelled
# Consistent with RA population colour scheme throughout
# ============================================================

library(vegan)

# ── Population colours — confirm consistent with all other figures ────────────
pop_cols_RA <- c(
  "MZ" = "#0072B2",
  "TZ" = "#E69F00",
  "MD" = "#56B4E9",
  "RI" = "#D55E00",
  "SY" = "#009E73",
  "MA" = "#FF69B4",
  "AU" = "#8B008B"
)

pop_labels_RA <- c(
  "MZ" = "Mozambique",
  "TZ" = "Tanzania",
  "MD" = "Madagascar",
  "RI" = "Réunion Island",
  "SY" = "Seychelles",
  "MA" = "Malaysia",
  "AU" = "Australia"
)

# ── Scores ────────────────────────────────────────────────────────────────────
site_scores <- scores(rda_partial, display = "sites",
                       scaling = 2, choices = 1:3)
bp_scores   <- scores(rda_partial, display = "bp",
                       scaling = 2)
sp_scores   <- scores(rda_partial, display = "species",
                       scaling = 2, choices = 1:3)

# Axis variance percentages
eig_vals <- summary(eigenvals(rda_partial,
                               model = "constrained"))
ax1_var  <- round(eig_vals["Proportion Explained","RDA1"]*100, 1)
ax2_var  <- round(eig_vals["Proportion Explained","RDA2"]*100, 1)
ax3_var  <- round(eig_vals["Proportion Explained","RDA3"]*100, 1)
cat("RDA1:", ax1_var, "% | RDA2:", ax2_var,
    "% | RDA3:", ax3_var, "%\n")

# Arrow scale — 60% of plot range
env_labels <- c(
  "ph_mean"      = "pH",
  "thetao_range" = "Temp range",
  "phyc_range"   = "Phyto range",
  "sws_mean"     = "Current mean"
)

# ── Triplot function ──────────────────────────────────────────────────────────
make_triplot_RA <- function(ax_x, ax_y, ax_x_var, ax_y_var,
                             filename, labelled = TRUE) {

  # Extract scores for the specific axes requested
  site_sc <- scores(rda_partial, display = "sites",
                    scaling = 2,
                    choices = c(which(c("RDA1","RDA2","RDA3") == ax_x),
                                which(c("RDA1","RDA2","RDA3") == ax_y)))
  bp_sc   <- scores(rda_partial, display = "bp",
                    scaling = 2,
                    choices = c(which(c("RDA1","RDA2","RDA3") == ax_x),
                                which(c("RDA1","RDA2","RDA3") == ax_y)))
  sp_sc   <- scores(rda_partial, display = "species",
                    scaling = 2,
                    choices = c(which(c("RDA1","RDA2","RDA3") == ax_x),
                                which(c("RDA1","RDA2","RDA3") == ax_y)))

  ind_x <- site_sc[, 1]
  ind_y <- site_sc[, 2]
  bp_x  <- bp_sc[, 1]
  bp_y  <- bp_sc[, 2]
  sp_x  <- sp_sc[, 1]
  sp_y  <- sp_sc[, 2]

  x_lim <- max(abs(c(sp_x, ind_x,
                      bp_x * 1.3))) * 1.2
  y_lim <- max(abs(c(sp_y, ind_y,
                      bp_y * 1.3))) * 1.2

  arrow_scale <- 0.6 * min(x_lim, y_lim) /
                 max(sqrt(bp_x^2 + bp_y^2))

  for (ext in c("pdf","png")) {
    if (ext == "pdf") cairo_pdf(paste0(filename, ".pdf"),
                                 width = 7, height = 7)
    else png(paste0(filename, ".png"),
             width = 7, height = 7,
             units = "in", res = 300)

    # 1. Grey SNP cloud — bottom layer
    plot(sp_x, sp_y,
         pch     = 19, cex = 0.3,
         col     = adjustcolor("grey70", alpha.f = 0.4),
         xlim    = c(-x_lim, x_lim),
         ylim    = c(-y_lim, y_lim),
         xlab    = paste0(ax_x, " (",
                           ax_x_var,
                           "% of constrained var.)"),
         ylab    = paste0(ax_y, " (",
                           ax_y_var,
                           "% of constrained var.)"),
         cex.lab  = 1.1, cex.axis = 1.0,
         las = 1, bty = "o")
    abline(h = 0, v = 0,
           lty = 2, col = "grey70", lwd = 0.5)

    # 2. Individual points — middle layer
    for (pop in POP_ORDER) {
      idx <- which(pop_factor == pop)
      points(ind_x[idx], ind_y[idx],
             pch = 21, cex = 1.8,
             bg  = pop_cols_RA[pop],
             col = "black", lwd = 0.8)
    }

    # 3. Environmental arrows — top layer
    # Drawn last so they appear above individual points
    for (i in 1:length(bp_x)) {
      arrows(0, 0,
             bp_x[i] * arrow_scale,
             bp_y[i] * arrow_scale,
             length = 0.12, lwd = 2.0,
             col = "black")
      if (labelled) {
        text(bp_x[i] * arrow_scale * 1.2,
             bp_y[i] * arrow_scale * 1.2,
             labels = env_labels[rownames(bp_sc)[i]],
             cex = 0.85, font = 3,
             col = "black")
      }
    }

    # 4. Legend
    if (labelled) {
      legend("topright",
             legend = pop_labels_RA[POP_ORDER],
             pt.bg  = pop_cols_RA[POP_ORDER],
             col    = "black", pch = 21,
             pt.cex = 1.5, pt.lwd = 0.8,
             bty = "n", cex = 0.85)
    }
    dev.off()
  }
  cat(filename, "saved\n")
}

# ── Panel A: RDA1 vs RDA2 ─────────────────────────────────────────────────────
make_triplot_RA("RDA1","RDA2", ax1_var, ax2_var,
                "RA_set1_triplot_RDA1_RDA2_labelled",
                labelled = TRUE)
make_triplot_RA("RDA1","RDA2", ax1_var, ax2_var,
                "RA_set1_triplot_RDA1_RDA2_unlabelled",
                labelled = FALSE)

# ── Panel B: RDA1 vs RDA3 ─────────────────────────────────────────────────────
make_triplot_RA("RDA1","RDA3", ax1_var, ax3_var,
                "RA_set1_triplot_RDA1_RDA3_labelled",
                labelled = TRUE)
make_triplot_RA("RDA1","RDA3", ax1_var, ax3_var,
                "RA_set1_triplot_RDA1_RDA3_unlabelled",
                labelled = FALSE)

cat("Section 13 — all triplot panels saved\n")

# ============================================================
# SECTION 14: Top 10 candidate loci
# 10 loci with highest max |r| across all 4 env variables
# Selected from RDA environmental candidates (n=846)
# Points coloured by population — all 7 populations included
# FDR-corrected p-values shown (BH, Section 9)
# ============================================================

# ── Use cor_r already computed in Section 9 ──────────────────────────────────
# cor_r = 846 loci x 4 env variables Pearson r matrix
# cor_p_adj = BH FDR-corrected p-values
# cor_summary = per-locus summary with max_r and best_env

stopifnot(exists("cor_r"))
stopifnot(exists("cor_p_adj"))
stopifnot(exists("cor_summary"))
cat("Correlation objects confirmed ✓\n")

# ── Select top 10 loci by max |r| across all variables ───────────────────────
max_r_per_locus <- apply(abs(cor_r), 1, max)
top10_idx       <- order(max_r_per_locus, decreasing = TRUE)[1:10]
top10_loci      <- rownames(cor_r)[top10_idx]
top10_best_var  <- env_final_4[
  apply(abs(cor_r[top10_idx, ]), 1, which.max)]
top10_r         <- round(max_r_per_locus[top10_idx], 3)
top10_p_adj     <- sapply(seq_along(top10_loci), function(i) {
  cor_p_adj[top10_loci[i], top10_best_var[i]]
})

cat("Top 10 candidate loci:\n")
print(data.frame(
  rank    = 1:10,
  locus   = top10_loci,
  best_var = top10_best_var,
  max_r   = top10_r,
  p_adj   = formatC(top10_p_adj, format = "e", digits = 2)
))

# ── Variable labels ───────────────────────────────────────────────────────────
var_labels <- c(
  "ph_mean"      = "pH",
  "thetao_range" = "Temp range",
  "phyc_range"   = "Phyto range",
  "sws_mean"     = "Current mean"
)

# ── Point colours — all 7 populations ────────────────────────────────────────
pt_cols <- pop_cols_RA[as.character(pop_factor)]

# ── Environmental matrix — final 4 variables ─────────────────────────────────
env_mat_final <- as.matrix(
  env_data_11[match(popmap$ID, env_data_11$ID),
               env_vars_final])
rownames(env_mat_final) <- popmap$ID
stopifnot(all(rownames(env_mat_final) == rownames(G_imp)))

# ── Plot ──────────────────────────────────────────────────────────────────────
for (ext in c("pdf","png")) {
  if (ext == "pdf") cairo_pdf("RA_set1_top10_candidate_loci.pdf",
                               width = 10, height = 14)
  else png("RA_set1_top10_candidate_loci.png",
           width = 10, height = 14,
           units = "in", res = 300)

  par(mfrow = c(5, 2), mar = c(4, 4, 3, 1))

  for (i in 1:10) {
    locus   <- top10_loci[i]
    best_v  <- top10_best_var[i]
    env_var <- env_mat_final[, best_v]
    geno    <- G_imp[, locus]
    r_val   <- round(cor_r[locus, best_v], 3)
    p_val   <- top10_p_adj[i]

    plot(env_var, geno,
         pch     = 21,
         bg      = pt_cols,
         col     = "black",
         cex     = 1.0,
         lwd     = 0.5,
         xlab    = var_labels[best_v],
         ylab    = "Genotype (0/1/2)",
         main    = paste0("Locus ", i,
                           " (", var_labels[best_v], ")",
                           "\nr = ", r_val,
                           ", FDR p = ",
                           formatC(p_val, format = "e",
                                   digits = 2)),
         cex.main = 0.8,
         cex.lab  = 0.9,
         las      = 1,
         bty      = "l")
    abline(lm(geno ~ env_var),
           col = "red", lwd = 1.5)
  }
  dev.off()
}
cat("Top 10 candidate loci plot saved\n")
cat("Section 14 complete\n")

# ── FINAL CHECKPOINT — all objects ──────────────────────────────────
save(
  # Section 1 — raw inputs
  popmap, coords, G, G_imp, locus_names_full,
  pop_factor, gind_full,
  # Section 2-3 — outlier methods
  outliers_outflank, outliers_pcadapt_names, gif_values,
  # Section 4 — environmental data
  env_data_11, env_vars_11,
  # Section 5 — spatial
  marine_dist, dbmems, dbmem_df,
  # Section 6-8 — RDA pipeline
  selected_env, env_final, env_final_4,
  rda_partial, anova_partial, anova_partial_axes, vp,
  # Section 9 — RDA outliers
  outliers_rda, cor_r, cor_p_adj, cor_summary, cor_table,
  # Section 10 — LFMM
  outliers_lfmm_names, gif, mod_lfmm,
  # Section 11 — four datasets
  outliers_neutral, outliers_any, outliers_enviro,
  outliers_highconf,
  G_neutral, G_any, G_enviro, G_highconf,
  # Figure objects
  pop_cols_RA, pop_labels_RA,
  file = "RA_seascape_FINAL_checkpoint.RData"
)
cat("FINAL CHECKPOINT SAVED ✓\n")
