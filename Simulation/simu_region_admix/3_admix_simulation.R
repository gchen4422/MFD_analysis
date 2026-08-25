# Admixed simulation for multi-ancestry fine-mapping benchmarking
# Generates ancestry-specific summary statistics from admixed individuals
# Key difference from homogeneous (shared_50): single phenotype, z-scores are correlated
#
# Naming convention (aligned with shared_50):
#   EUR = ancestry 1 (anc0 in original pipeline) -> zscore_1
#   AFR = ancestry 2 (anc1 in original pipeline) -> zscore_2
#
# Z-score model (admixed, differs from homogeneous):
#   E[z_EUR] = sqrt(N) * (R_EUR_EUR * beta_EUR + R_EUR_AFR * beta_AFR)
#   E[z_AFR] = sqrt(N) * (R_AFR_EUR * beta_EUR + R_AFR_AFR * beta_AFR)
# where R_EUR_AFR = X_EUR^T X_AFR / N is the cross-ancestry LD.

library(data.table)
library(mvtnorm)

# --- Constants ---
h2_locus <- 0.005  # per-locus heritability (CARMAX setting, compensates for small N)

geno_dir_eur <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/real_data_eur_100regions"
geno_dir_afr <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/real_data_afr_100regions"

# --- Helper: standardize genotype columns ---
standardize_geno <- function(X) {
  apply(X, 2, function(x) {
    (x - mean(x)) / sd(x)
  })
}

for (LD_BLOCK in 1:100) {

  # Seed initialization (same scheme as shared_50)
  set.seed(1128)
  seed_LOCIrator <- sample.int(1000, 100)[LD_BLOCK]
  set.seed(seed_LOCIrator)

  # --- Load ancestry-specific genotypes from RDS (SNP x individual) ---
  EUR_raw <- readRDS(file.path(geno_dir_eur, paste0("region", LD_BLOCK, ".rds")))
  AFR_raw <- readRDS(file.path(geno_dir_afr, paste0("region", LD_BLOCK, ".rds")))

  # SNP identifiers (CHR:POS format)
  snp_names <- rownames(EUR_raw)
  snp_parts <- strsplit(snp_names, ":")
  snp_chr <- sapply(snp_parts, `[`, 1)
  snp_pos <- as.integer(sapply(snp_parts, `[`, 2))

  # Transpose to individual x SNP
  EU_geno <- t(EUR_raw)
  BB_geno <- t(AFR_raw)
  colnames(EU_geno) <- snp_names
  colnames(BB_geno) <- snp_names

  # --- Remove zero-SD SNPs in either ancestry ---
  sd_EUR <- apply(EU_geno, 2, sd)
  sd_AFR <- apply(BB_geno, 2, sd)
  keep <- which(sd_EUR > 0 & sd_AFR > 0)

  if (length(keep) < 10) {
    cat("Region", LD_BLOCK, ": too few non-zero-SD SNPs, skipping\n")
    next
  }

  EU_geno <- EU_geno[, keep, drop = FALSE]
  BB_geno <- BB_geno[, keep, drop = FALSE]
  snp_names <- snp_names[keep]
  snp_chr <- snp_chr[keep]
  snp_pos <- snp_pos[keep]

  # Standardize (mean 0, sd 1 per column)
  EU_c <- standardize_geno(EU_geno)
  BB_c <- standardize_geno(BB_geno)

  N <- nrow(EU_c)
  N_SNP <- ncol(EU_c)

  # --- Compute LD matrices (once per region, reused across causal scenarios) ---
  EU_cov <- t(EU_c) %*% EU_c / N       # R_EUR_EUR: within-EUR LD
  BB_cov <- t(BB_c) %*% BB_c / N       # R_AFR_AFR: within-AFR LD
  cross_cov <- t(EU_c) %*% BB_c / N    # R_EUR_AFR: cross-ancestry LD

  for (num_causal in 1:3) {

    wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/causal_num_", num_causal, "/")
    data_dir <- paste0(wrk_dir, "summary_data/")
    result_dir <- paste0(wrk_dir, "result/")
    out_dir <- paste0(wrk_dir, "out/")
    system(paste0("mkdir -p ", data_dir))
    system(paste0("mkdir -p ", result_dir))
    system(paste0("mkdir -p ", out_dir))

    # Number of causal SNPs: 2, 6, 10
    num_causal_SNP <- c(2, 4, 6)[num_causal]

    # Stochastic causal config: 1=EUR-only(25%), 2=AFR-only(25%), 3=shared(50%)
    causal_SNP_config <- sample(c(1, 2, 3), num_causal_SNP, replace = TRUE, prob = c(0.25, 0.25, 0.5))

    # Sample causal SNP indices (uniform prior)
    causal_idx <- sample(seq_len(N_SNP), num_causal_SNP, replace = FALSE)

    # --- Draw raw betas ---
    beta_EUR_all <- rep(0, N_SNP)
    beta_AFR_all <- rep(0, N_SNP)
    signal <- rep(0, N_SNP)

    for (ci in seq_len(num_causal_SNP)) {
      idx <- causal_idx[ci]

      if (causal_SNP_config[ci] == 1) {
        # EUR-only: only EUR ancestry has effect
        beta_EUR_all[idx] <- rnorm(1, 0, 1)
        signal[idx] <- 1

      } else if (causal_SNP_config[ci] == 2) {
        # AFR-only: only AFR ancestry has effect
        beta_AFR_all[idx] <- rnorm(1, 0, 1)
        signal[idx] <- 2

      } else {
        # Shared: correlated betas (rho=0.8)
        beta_shared <- rmvnorm(1, mean = rep(0, 2),
                               sigma = matrix(c(1, 0.8, 0.8, 1), ncol = 2, nrow = 2))
        beta_EUR_all[idx] <- beta_shared[1]
        beta_AFR_all[idx] <- beta_shared[2]
        signal[idx] <- 3
      }
    }

    # --- Scale jointly to hit h2_locus ---
    # Single phenotype: y = X_EUR * beta_EUR + X_AFR * beta_AFR + e
    # Total genetic variance from this locus must equal h2_locus
    genetic_value <- EU_c %*% beta_EUR_all + BB_c %*% beta_AFR_all
    var_g <- drop(var(genetic_value))
    scale_factor <- sqrt(h2_locus / var_g)
    beta_EUR_all <- beta_EUR_all * scale_factor
    beta_AFR_all <- beta_AFR_all * scale_factor

    # --- Generate single admixed phenotype ---
    e <- rnorm(N, 0, sqrt(1 - h2_locus))
    y <- EU_c %*% beta_EUR_all + BB_c %*% beta_AFR_all + e

    # --- Compute z-scores via direct projection ---
    # z_EUR_j = X_EUR[:,j]^T y / sqrt(N)
    # This naturally captures:
    #   E[z_EUR] = sqrt(N) * (R_EUR_EUR * beta_EUR + R_EUR_AFR * beta_AFR)
    #   E[z_AFR] = sqrt(N) * (R_AFR_EUR * beta_EUR + R_AFR_AFR * beta_AFR)
    # The cross-ancestry LD leakage and correlated noise are both captured
    # because the same y is projected onto both ancestry genotype matrices.
    z_EU <- drop(t(EU_c) %*% y) / sqrt(N)
    z_BB <- drop(t(BB_c) %*% y) / sqrt(N)

    # --- Write z-file (same format as shared_50) ---
    z_file <- data.frame(
      CHR = as.numeric(snp_chr[1]),
      POS = snp_pos,
      RSID = snp_names,
      zscore_1 = z_EU,
      zscore_2 = z_BB,
      Signal = signal,
      N_1 = N,
      N_2 = N
    )

    z_file_out <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_1")
    write.table(z_file, z_file_out, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " ")

    # --- Write LD matrices ---
    fwrite(as.data.table(as.matrix(EU_cov)),
           paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_1.LD1"),
           sep = " ", col.names = FALSE)
    fwrite(as.data.table(as.matrix(BB_cov)),
           paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_1.LD2"),
           sep = " ", col.names = FALSE)
    fwrite(as.data.table(as.matrix(cross_cov)),
           paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_1.LD_cross"),
           sep = " ", col.names = FALSE)

    print(paste("Done: region", LD_BLOCK, "causal_num", num_causal))
  }
}
