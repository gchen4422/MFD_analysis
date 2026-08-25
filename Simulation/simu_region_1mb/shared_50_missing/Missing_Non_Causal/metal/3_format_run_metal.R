library(data.table)
library(dplyr)

safe_read <- function(p) {
  if (!file.exists(p)) return(NULL)
  fread(p, showProgress = FALSE)
}


# ranges used in this analysis
num_causal_vec      <- 1:3          # e.g. 1,2
h2_num_vec          <- 1:2          # e.g. 1,2
LD_BLOCK_range      <- 1:100        # or c(5,12,25,33,54,56,63,82,87)
causal_index_vec    <- 1:2          # 1 = Both, 2 = One
External_index_vec  <- 1:2          # 1 = "", 2 = "External_"

for (causal_index in causal_index_vec) {
  for (External_index in External_index_vec) {
    
    causal_index_name   <- c("Both", "One")[causal_index]
    External_index_name <- c("", "External_")[External_index]
    
    for (num in num_causal_vec) {
      for (h2 in h2_num_vec) {
        
        # work/result dirs do NOT depend on LD_BLOCK
        wrk_dir <- paste0(
          "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/",
          "shared_50_missing/Missing_Non_Causal/",
          External_index_name, causal_index_name,
          "/causal_num_", num, "/"
        )
        data_dir   <- paste0(wrk_dir, "summary_data/")
        result_dir <- paste0(wrk_dir, "result/")
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        
        # accumulators across ALL LD_BLOCK for this (num, h2)
        eur_parts <- list()
        afr_parts <- list()
        k_eur <- 1L
        k_afr <- 1L
        
        for (LD_BLOCK in LD_BLOCK_range) {
          
          # read zfile for this locus
          zfile_path <- paste0(
            data_dir,
            "CAUSAL_", num, "_LOCI_", LD_BLOCK, "_h2_", h2
          )
          zfile <- safe_read(zfile_path)
          
          summary_stat_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_",LD_BLOCK,".bim")) %>% filter(V2 %in% zfile$RSID)
          
          zfile = zfile %>% mutate(A1 = summary_stat_bim$V5,A2 = summary_stat_bim$V6)
          

          if (is.null(zfile) || nrow(zfile) == 0L) {
            message("[zfile] missing or empty: ", zfile_path)
            next
          }
          
          zfile <- as.data.frame(zfile)
          
          # ----- missingness filtering -----
          if (causal_index == 1) {
            # drop SNPs missing in EU or BB
            EU_zfile <- zfile %>% filter(EU_missing == 0)
            BB_zfile <- zfile %>% filter(BB_missing == 0)
            missing_EU_index <- which(zfile$EU_missing != 0)
            missing_BB_index <- which(zfile$BB_missing != 0)
            
          } else if (causal_index == 2) {
            # drop SNPs missing in EU only
            EU_zfile <- zfile %>% filter(EU_missing == 0)
            BB_zfile <- zfile
            missing_EU_index <- which(zfile$EU_missing != 0)
            missing_BB_index <- integer(0)
            
          } else {
            message("[causal_index] invalid value: ", causal_index)
            next
          }
          
          # ----- add meta columns and store for combination -----
          if (nrow(EU_zfile) > 0L) {
            EU_zfile$Beta       <- EU_zfile$zscore_1/sqrt(EU_zfile$N_1)
            EU_zfile$Se       <- 1/sqrt(EU_zfile$N_1)
            EU_zfile$loci       <- as.integer(LD_BLOCK)
            EU_zfile$num_causal <- num
            EU_zfile$h2         <- h2
            eur_parts[[k_eur]]  <- EU_zfile
            k_eur <- k_eur + 1L
          }
          
          if (nrow(BB_zfile) > 0L) {
            BB_zfile$Beta       <- BB_zfile$zscore_2/sqrt(BB_zfile$N_2)
            BB_zfile$Se       <- 1/sqrt(BB_zfile$N_2)
            BB_zfile$loci       <- as.integer(LD_BLOCK)
            BB_zfile$num_causal <- num
            BB_zfile$h2         <- h2
            afr_parts[[k_afr]]  <- BB_zfile
            k_afr <- k_afr + 1L
          }
          
          # (per-locus analysis can be performed here if needed)
          
        } # end LD_BLOCK loop
        
        # ---------- EUR combined across loci ----------
        eur_parts <- eur_parts[!vapply(eur_parts, is.null, logical(1))]
        if (length(eur_parts) > 0L) {
          eur_combined <- rbindlist(eur_parts, use.names = TRUE, fill = TRUE)
          out_eur <- file.path(
            data_dir,
            sprintf("combined_%s_h2_%s_eur.tsv", num, as.character(h2))
          )
          fwrite(eur_combined, out_eur, sep = "\t")
          message(sprintf("Wrote %s (n=%s rows)", out_eur, nrow(eur_combined)))
        } else {
          warning(sprintf("No EUR inputs for num=%s, h2=%s", num, h2))
        }
        
        # ---------- AFR combined across loci ----------
        afr_parts <- afr_parts[!vapply(afr_parts, is.null, logical(1))]
        if (length(afr_parts) > 0L) {
          afr_combined <- rbindlist(afr_parts, use.names = TRUE, fill = TRUE)
          out_afr <- file.path(
            data_dir,
            sprintf("combined_%s_h2_%s_afr.tsv", num, as.character(h2))
          )
          fwrite(afr_combined, out_afr, sep = "\t")
          message(sprintf("Wrote %s (n=%s rows)", out_afr, nrow(afr_combined)))
        } else {
          warning(sprintf("No AFR inputs for num=%s, h2=%s", num, h2))
        }
        
      } # h2
    }   # num
  }     # External_index
}       # causal_index
    
    
    






write_metal_per_combo <- function(num_causal_vec,
                                  h2_num_vec,
                                  causal_index_vec,
                                  External_index_vec,
                                  file_ext = ".tsv") {
  
  for (causal_index in causal_index_vec) {
    for (External_index in External_index_vec) {
      
      causal_index_name   <- c("Both", "One")[causal_index]
      External_index_name <- c("", "External_")[External_index]
      
      for (num in num_causal_vec) {
        for (h2 in h2_num_vec) {
          
          # paths for the combined EUR and AFR inputs
          wrk_dir <- paste0(
            "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/",
            "shared_50_missing/Missing_Non_Causal/",
            External_index_name, causal_index_name,
            "/causal_num_", num, "/"
          )
          
          data_dir   <- file.path(wrk_dir, "summary_data")
          result_dir <- file.path(wrk_dir, "result")
          dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
          
          # combined inputs (EUR/AFR) for this (causal_index, External_index, num, h2)
          eur_file <- file.path(
            data_dir,
            sprintf("combined_%s_h2_%s_eur%s", num, as.character(h2), file_ext)
          )
          afr_file <- file.path(
            data_dir,
            sprintf("combined_%s_h2_%s_afr%s", num, as.character(h2), file_ext)
          )
          
          # skip if inputs are missing
          if (!file.exists(eur_file) || !file.exists(afr_file)) {
            warning("Missing EUR/AFR combined file for: ",
                    "causal_index=", causal_index_name,
                    ", External_index=", External_index_name,
                    ", num=", num, ", h2=", h2)
            next
          }
          
          # where to write the METAL script (per combo, same name within its folder)
          metal_txt <- file.path(
            data_dir,
            sprintf("combined_%s_h2_%s_metal.txt", num, as.character(h2))
          )
          
          out_prefix <- file.path(
            result_dir,
            sprintf("Metal_%s_h2_%s_result", num, as.character(h2))
          )
          
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
            sprintf("OUTFILE %s .tbl", out_prefix),
            "ANALYZE HETEROGENEITY"
          )
          
          writeLines(lines, metal_txt)
          message("Wrote METAL script: ", metal_txt)
        }
      }
    }
  }
}


write_metal_per_combo(num_causal_vec,h2_num_vec,causal_index_vec, External_index_vec, file_ext = ".tsv")


# Path to METAL binary
metal_bin <- "~/gwas_software/METAL/generic-metal/executables/metal"

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal"

num_causal_vec     <- c(1,2,3)
h2_num_vec         <- c(1, 2)
causal_index_names <- c("Both", "One")
external_names     <- c("", "External_")

combos <- expand.grid(
  causal  = causal_index_names,
  external = external_names,
  num     = num_causal_vec,
  h2      = h2_num_vec,
  stringsAsFactors = FALSE
)

# pattern: "Both", "External_Both", "One", "External_One"
combos$pattern <- ifelse(
  combos$external == "",
  combos$causal,
  paste0(combos$external, combos$causal)
)

metal_scripts <- with(
  combos,
  sprintf(
    "%s/%s/causal_num_%d//summary_data/combined_%d_h2_%d_metal.txt",
    base_dir, pattern, num, num, h2
  )
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

