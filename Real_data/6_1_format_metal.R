library(data.table)


traits_list <- c("BMI", "DBP", "SBP")


write_metal_per_combo <- function(traits_list, file_ext = ".tsv") {
  for (trait in traits_list) {
    
      # input and output paths
      wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_", trait, "/")
      data_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/")
      result_dir <- paste0(wrk_dir, "result/")

      # combined inputs (EUR/AFR) for this (num, h2)
      eur_file <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv")
      afr_file <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv")
      
      # where to write the METAL script
      metal_txt <- file.path(data_dir, sprintf("%s_metal.txt", as.character(trait)))
      
      out_prefix <- file.path(result_dir, sprintf("Metal_%s_result", as.character(trait)))
      
      lines <- c(
        "# Define the separator for CSV files",
        "SEPARATOR COMMA",
        
        "# classical approach, uses effect size estimates and standard errors",
        "SCHEME STDERR",
        
        "# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===",
        "MARKER CHR_POS",
        "ALLELE ALT REF",
        "EFFECT BETA",
        "PVALUE PVAL",
        "STDERR SE",
        
        # --- ADD THESE LINES ---
        "CUSTOMVARIABLE TotalN",   # Define a new variable for the output
        "LABEL TotalN AS N",       # Map the input file's 'N' column to it
        # -----------------------
        
        sprintf("PROCESS %s", eur_file),
        
        "",
        "# === THE SECOND INPUT FILE HAS THE SAME FORMAT ===",
        sprintf("PROCESS %s", afr_file),
        
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
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/BMI_metal.txt",
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/DBP_metal.txt",
  "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/SBP_metal.txt"
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

