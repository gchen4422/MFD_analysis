library(data.table)


LD_BLOCK = c(1:100)
h2_num = c(1,2)
num_causal = c(1,2,3)
# assumes: num_causal, h2_num, LD_BLOCK (e.g., LD_BLOCK <- 1:100)

safe_read <- function(p) {
  if (!file.exists(p)) return(NULL)
  fread(p, showProgress = FALSE)
}

for (num in num_causal){
  for (h2 in h2_num){

    # input and output paths
    wrk_dir    <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_", num, "/")
    data_dir   <- paste0(wrk_dir, "summary_data/")
    result_dir <- paste0(wrk_dir, "result/")
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

    # ---------- EUR ----------
    eur_parts <- vector("list", length(LD_BLOCK))
    k <- 1L
    for (loci in LD_BLOCK){
      f_eur <- paste0(data_dir, "CAUSAL_", num, "_LOCI_", loci, "_h2_", h2, "_eur")
      dt <- safe_read(f_eur)
      if (!is.null(dt)) {
        dt[, `:=`(loci = as.integer(loci), num_causal = num, h2 = h2)]
        eur_parts[[k]] <- dt
      }
      k <- k + 1L
    }
    eur_parts <- eur_parts[!vapply(eur_parts, is.null, logical(1))]
    if (length(eur_parts) > 0L) {
      eur_combined <- rbindlist(eur_parts, use.names = TRUE, fill = TRUE)
      out_eur <- file.path(data_dir, sprintf("combined_%s_h2_%s_eur.tsv", num, as.character(h2)))
      fwrite(eur_combined, out_eur, sep = "\t")
      message(sprintf("Wrote %s (n=%s rows)", out_eur, nrow(eur_combined)))
    } else {
      warning(sprintf("No EUR inputs for num=%s, h2=%s", num, h2))
    }

    # ---------- AFR ----------
    afr_parts <- vector("list", length(LD_BLOCK))
    k <- 1L
    for (loci in LD_BLOCK){
      f_afr <- paste0(data_dir, "CAUSAL_", num, "_LOCI_", loci, "_h2_", h2, "_afr")
      dt <- safe_read(f_afr)
      if (!is.null(dt)) {
        dt[, `:=`(loci = as.integer(loci), num_causal = num, h2 = h2)]
        afr_parts[[k]] <- dt
      }
      k <- k + 1L
    }
    afr_parts <- afr_parts[!vapply(afr_parts, is.null, logical(1))]
    if (length(afr_parts) > 0L) {
      afr_combined <- rbindlist(afr_parts, use.names = TRUE, fill = TRUE)
      out_afr <- file.path(data_dir, sprintf("combined_%s_h2_%s_afr.tsv", num, as.character(h2)))
      fwrite(afr_combined, out_afr, sep = "\t")
      message(sprintf("Wrote %s (n=%s rows)", out_afr, nrow(afr_combined)))
    } else {
      warning(sprintf("No AFR inputs for num=%s, h2=%s", num, h2))
    }
  }
}




write_metal_per_combo <- function(num_causal, h2_num, file_ext = ".tsv") {
  for (num in num_causal) {
    for (h2 in h2_num) {
      # input and output paths
      wrk_dir    <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_", num, "/")
      data_dir   <- paste0(wrk_dir, "summary_data/")  # retained for directory consistency
      result_dir <- paste0(wrk_dir, "result/")

      # combined inputs (EUR/AFR) for this (num, h2)
      eur_file <- file.path(data_dir, sprintf("combined_%s_h2_%s_eur%s", num, as.character(h2), file_ext))
      afr_file <- file.path(data_dir, sprintf("combined_%s_h2_%s_afr%s", num, as.character(h2), file_ext))

      # where to write the METAL script
      metal_txt <- file.path(data_dir, sprintf("combined_%s_h2_%s_metal.txt", num, as.character(h2)))

      out_prefix <- file.path(result_dir, sprintf("Metal_%s_h2_%s_result", num, as.character(h2)))

      lines <- c(
        "# classical approach, uses effect size estimates and standard errors",
        "SCHEME STDERR",
        "",
        "# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===",
        "MARKER RSID",
        "ALLELE A1 A2",
        "EFFECT Beta",
        "PVALUE P",
        "STDERR Se",
        sprintf("PROCESS %s", eur_file),
        "",
        "# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===",
        sprintf("PROCESS %s", afr_file),
        "",
        sprintf("OUTFILE %s .tbl", out_prefix),   # METAL output prefix
        "ANALYZE HETEROGENEITY"                   # Run heterogeneity analysis
      )

      writeLines(lines, metal_txt)
      message("Wrote METAL script: ", metal_txt)
    }
  }
}


write_metal_per_combo(num_causal = c(1,2,3), h2_num = c(1,2), file_ext = ".tsv")


# Path to METAL binary
metal_bin <- "~/gwas_software/METAL/generic-metal/executables/metal"

# The six METAL scripts used by the workflow
metal_scripts <- c(
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_1/summary_data/combined_1_h2_1_metal.txt",
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_1/summary_data/combined_1_h2_2_metal.txt",
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_2/summary_data/combined_2_h2_1_metal.txt",
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_2/summary_data/combined_2_h2_2_metal.txt",
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_3/summary_data/combined_3_h2_1_metal.txt",
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_3/summary_data/combined_3_h2_2_metal.txt"
)

# Run each script; write stdout/stderr to a .log next to the script
for (f in metal_scripts) {
  if (!file.exists(f)) {
    message("Missing: ", f)
    next
  }
  message("Running METAL: ", f)
  command = paste0(metal_bin, " ", f)
  system(command)
}
