args <- as.integer(commandArgs(trailingOnly = TRUE))
if (length(args) != 3L || anyNA(args)) {
  stop("Usage: Rscript 4_shared_50_analysis.R <region> <h2> <num_causal>", call. = FALSE)
}
LD_BLOCK <- args[1]
h2_num <- args[2]
num_causal <- args[3]

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(susieR))

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix"
wrk_dir <- file.path(base_dir, paste0("causal_num_", num_causal))
data_dir <- file.path(wrk_dir, "summary_data")
result_dir <- file.path(wrk_dir, "result")
merged_ld_dir <- file.path(base_dir, "merged_ld_eur_afr")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- sprintf("CAUSAL_%d_LOCI_%d_h2_%d", num_causal, LD_BLOCK, h2_num)
z_path <- file.path(data_dir, prefix)
metal_path <- file.path(result_dir, sprintf("Metal_%d_h2_%d_result1.tbl", num_causal, h2_num))
merged_ld_path <- file.path(merged_ld_dir, sprintf("region_%d.ld", LD_BLOCK))
merged_snp_path <- file.path(merged_ld_dir, sprintf("region_%d.snps", LD_BLOCK))

required_paths <- c(z_path, metal_path, merged_ld_path, merged_snp_path)
if (!all(file.exists(required_paths))) {
  stop("Missing required input(s): ", paste(required_paths[!file.exists(required_paths)], collapse = ", "), call. = FALSE)
}

normalise_ld <- function(R, label) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) stop(label, " must be square.", call. = FALSE)
  R <- (R + t(R)) / 2
  d <- diag(R)
  if (any(!is.finite(d)) || any(d <= 0)) stop(label, " has an invalid diagonal.", call. = FALSE)
  R <- sweep(sweep(R, 1L, sqrt(d), "/"), 2L, sqrt(d), "/")
  diag(R) <- 1
  R
}

# Simulation SNP order; this is the order used by LD1 and LD2.
zfile <- fread(z_path)
if (!all(c("RSID", "N_1", "N_2") %in% names(zfile)) || anyDuplicated(zfile$RSID)) {
  stop("Simulation file must contain unique RSID, N_1, and N_2 columns.", call. = FALSE)
}
eur_path <- paste0(z_path, "_eur")
if (!file.exists(eur_path)) stop("Missing EUR allele file: ", eur_path, call. = FALSE)
eur_alleles <- fread(eur_path, select = c("RSID", "A1", "A2"))
if (anyDuplicated(eur_alleles$RSID)) stop("Duplicate RSIDs in EUR allele file.", call. = FALSE)
setnames(eur_alleles, "RSID", "SNP")
sim_snps <- data.table(SNP = zfile$RSID)[eur_alleles, on = "SNP"]
if (anyNA(sim_snps$A1) || anyNA(sim_snps$A2)) stop("Missing allele annotation for simulation SNPs.", call. = FALSE)
sim_snps[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]

# METAL output is combined across regions. Restrict and align by SNP ID.
metal <- fread(metal_path)
required_metal <- c("MarkerName", "Allele1", "Allele2", "Effect", "P-value")
if (!all(required_metal %in% names(metal))) stop("METAL output is missing required columns.", call. = FALSE)
metal <- metal[MarkerName %chin% sim_snps$SNP]
if (anyDuplicated(metal$MarkerName)) stop("Duplicate region SNPs in METAL output.", call. = FALSE)

# The pooled-admixed LD order is authoritative. It omits monomorphic total-dosage SNPs.
merged_snps <- fread(merged_snp_path)$SNP
merged_cov <- as.matrix(fread(merged_ld_path))
if (nrow(merged_cov) != ncol(merged_cov) || length(merged_snps) != nrow(merged_cov) || anyDuplicated(merged_snps)) {
  stop("Merged LD matrix and SNP-order file are inconsistent.", call. = FALSE)
}
analysis_snps <- merged_snps[merged_snps %chin% sim_snps$SNP & merged_snps %chin% metal$MarkerName]
if (length(analysis_snps) < 2L) stop("Fewer than two analysis SNPs overlap.", call. = FALSE)

sim_idx <- match(analysis_snps, sim_snps$SNP)
metal_idx <- match(analysis_snps, metal$MarkerName)
merged_idx <- match(analysis_snps, merged_snps)
target <- sim_snps[sim_idx]
metal <- metal[metal_idx]
metal[, `:=`(Allele1 = toupper(Allele1), Allele2 = toupper(Allele2))]
flip <- fifelse(metal$Allele1 == target$A1, 1,
               fifelse(metal$Allele1 == target$A2, -1, NA_real_))
if (anyNA(flip)) stop("Unresolvable METAL allele alignment.", call. = FALSE)
p <- pmax(as.numeric(metal[["P-value"]]), .Machine$double.xmin)
z_metal <- stats::qnorm(p / 2, lower.tail = FALSE) * sign(as.numeric(metal$Effect)) * flip
if (any(!is.finite(z_metal))) stop("Non-finite METAL Z score.", call. = FALSE)

# Persist the allele check so every METAL-to-LD alignment is inspectable.
allele_check <- data.table(
  SNP = analysis_snps,
  target_A1 = target$A1,
  target_A2 = target$A2,
  metal_Allele1 = metal$Allele1,
  metal_Allele2 = metal$Allele2,
  alignment = fifelse(flip == 1, "direct", "swapped"),
  flip = flip,
  metal_effect = as.numeric(metal$Effect),
  metal_pvalue = as.numeric(metal[["P-value"]]),
  Z_aligned = z_metal
)
allele_check_path <- file.path(
  result_dir,
  sprintf("Metal_SuSiE_CAUSAL_%d_LOCI_%d_h2_%d_allele_alignment.tsv",
          num_causal, LD_BLOCK, h2_num)
)
fwrite(allele_check, allele_check_path, sep = "\t")

# Ancestry-specific LD stays in simulation order; subset it to analysis_snps.
eur_cov <- as.matrix(fread(paste0(z_path, ".LD1")))
afr_cov <- as.matrix(fread(paste0(z_path, ".LD2")))
if (any(dim(eur_cov) != nrow(zfile)) || any(dim(afr_cov) != nrow(zfile))) {
  stop("Ancestry-specific LD dimensions do not match the simulation SNP list.", call. = FALSE)
}
sim_ld_idx <- match(analysis_snps, zfile$RSID)
eur_cov <- normalise_ld(eur_cov[sim_ld_idx, sim_ld_idx, drop = FALSE], "EUR LD")
afr_cov <- normalise_ld(afr_cov[sim_ld_idx, sim_ld_idx, drop = FALSE], "AFR LD")
merged_cov <- normalise_ld(merged_cov[merged_idx, merged_idx, drop = FALSE], "Merged admixed LD")

n_eur <- stats::median(zfile$N_1)
n_afr <- stats::median(zfile$N_2)
if (!is.finite(n_eur) || !is.finite(n_afr) || n_eur <= 0 || n_afr <= 0) stop("Invalid sample size.", call. = FALSE)
weighted_cov <- normalise_ld((n_eur * eur_cov + n_afr * afr_cov) / (n_eur + n_afr), "Weighted LD")

start_time <- Sys.time()
susie_weighted <- susie_rss(z = z_metal, R = weighted_cov, n = n_eur, check_prior = FALSE)
susie_weighted.time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
start_time <- Sys.time()
susie_merged <- susie_rss(z = z_metal, R = merged_cov, n = n_eur, check_prior = FALSE)
susie_merged.time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

output_path <- file.path(result_dir, sprintf("Metal_SuSiE_CAUSAL_%d_LOCI_%d_h2_%d_updated.RData", num_causal, LD_BLOCK, h2_num))
save(susie_weighted, susie_merged, susie_weighted.time, susie_merged.time, analysis_snps, file = output_path)
message("Saved ", output_path, " using ", length(analysis_snps), " SNPs.")
message("Saved METAL allele-alignment audit: ", allele_check_path)
