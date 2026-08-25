#!/usr/bin/env Rscript

# Generate LD for the observed admixed genotype at each region.
# EUR and AFR RDS files contain local-ancestry dosage components for the same
# individuals, so the merged genotype is EUR + AFR (not a row-bind of samples).

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bigstatsr))

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix"
eur_dir <- file.path(base_dir, "real_data_eur_100regions")
afr_dir <- file.path(base_dir, "real_data_afr_100regions")
out_dir <- file.path(base_dir, "merged_ld_eur_afr")

args <- commandArgs(trailingOnly = TRUE)
regions <- if (length(args) == 0L) 1:100 else as.integer(args)
if (anyNA(regions) || any(regions < 1L | regions > 100L)) {
  stop("Provide zero or more region IDs in 1:100.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (region in regions) {
  eur_file <- file.path(eur_dir, paste0("region", region, ".rds"))
  afr_file <- file.path(afr_dir, paste0("region", region, ".rds"))
  if (!file.exists(eur_file) || !file.exists(afr_file)) {
    stop("Missing EUR or AFR RDS input for region ", region, ".", call. = FALSE)
  }

  eur_raw <- readRDS(eur_file) # SNP x individual
  afr_raw <- readRDS(afr_file) # SNP x individual

  if (!identical(dim(eur_raw), dim(afr_raw)) ||
      !identical(rownames(eur_raw), rownames(afr_raw)) ||
      !identical(colnames(eur_raw), colnames(afr_raw))) {
    stop("EUR and AFR inputs are not aligned for region ", region, ".", call. = FALSE)
  }

  # Store combined dosages on disk.  This avoids retaining a dense double
  # genotype matrix in RAM for the largest regions.
  n_snp <- nrow(eur_raw)
  n_ind <- ncol(eur_raw)
  backingfile <- file.path(out_dir, paste0(".region_", region, "_geno"))
  geno <- FBM(n_ind, n_snp, type = "double", backingfile = backingfile)
  for (j in seq_len(n_snp)) {
    geno[, j] <- eur_raw[j, ] + afr_raw[j, ]
  }
  rm(eur_raw, afr_raw)
  invisible(gc())

  geno_stats <- big_colstats(geno)
  keep <- is.finite(geno_stats$var) & geno_stats$var > 0
  if (sum(keep) < 2L) {
    stop("Fewer than two variable SNPs in region ", region, ".", call. = FALSE)
  }

  ld <- big_cor(geno, ind.col = which(keep))
  ld <- (ld + t(ld)) / 2
  diag(ld) <- 1

  prefix <- file.path(out_dir, paste0("region_", region))
  fwrite(as.data.table(ld), paste0(prefix, ".ld"), sep = " ", col = FALSE)
  fwrite(data.table(SNP = rownames(readRDS(eur_file))[keep]),
         paste0(prefix, ".snps"), sep = "\t")
  unlink(paste0(backingfile, c(".bk", ".rds")))
  message("Wrote region ", region, ": ", sum(keep), " variable SNPs")
}
