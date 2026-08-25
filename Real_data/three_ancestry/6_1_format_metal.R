library(data.table)


traits_list <- c("BMI", "DBP", "SBP")


write_metal_per_combo <- function(traits_list, file_ext = ".tsv") {
  for (trait in traits_list) {

    # input and output paths
    wrk_dir   <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/")
    data_dir  <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/")
    result_dir <- paste0(wrk_dir, "result/")

    # combined inputs (EUR/AFR/EAS)
    eur_file <- paste0(data_dir, trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv")
    afr_file <- paste0(data_dir, trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv")
    eas_file <- paste0(data_dir, trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap.csv")

    # where to write the METAL script
    metal_txt  <- file.path(data_dir, sprintf("%s_metal.txt", trait))
    out_prefix <- file.path(result_dir, sprintf("Metal_%s_result", trait))

    lines <- c(
      "# Define the separator for CSV files",
      "SEPARATOR COMMA",
      "",
      "# Classical approach, uses effect size estimates and standard errors",
      "SCHEME STDERR",
      "",
      "# === COLUMN DESCRIPTIONS (same format for all input files) ===",
      "MARKER CHR_POS",
      "ALLELE ALT REF",
      "EFFECT BETA",
      "PVALUE PVAL",
      "STDERR SE",
      "",
      "# Track sample size",
      "CUSTOMVARIABLE TotalN",
      "LABEL TotalN AS N",
      "",
      "# === PROCESS INPUT FILES ===",
      sprintf("PROCESS %s", eur_file),
      sprintf("PROCESS %s", afr_file),
      sprintf("PROCESS %s", eas_file),
      "",
      sprintf("OUTFILE %s .tbl", out_prefix),
      "ANALYZE HETEROGENEITY"
    )

    writeLines(lines, metal_txt)
    message("Wrote METAL script: ", metal_txt)
  }
}


write_metal_per_combo(traits_list, file_ext = ".tsv")


# Path to METAL binary
metal_bin <- "~/gwas_software/METAL/generic-metal/executables/metal"

# The six METAL scripts used by the workflow
metal_scripts <- c(
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/BMI_metal.txt",
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/DBP_metal.txt",
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/SBP_metal.txt"
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
