library(data.table)
library(dplyr)
library(tidyr)

LD_BLOCK_vec  <- 1:100
h2_vec        <- c(1, 2)
num_causal_vec <- c(1, 2, 3)

for (num_causal in num_causal_vec) {
  
  # Directories that depend only on num_causal
  wrk_dir    <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/causal_num_", num_causal, "/")
  data_dir   <- paste0(wrk_dir, "summary_data/")
  result_dir <- paste0(wrk_dir, "result/")
  paintor_dir <- data_dir
  
  # make sure dirs exist
  system(paste0("mkdir -p ", wrk_dir))
  system(paste0("mkdir -p ", data_dir))
  system(paste0("mkdir -p ", result_dir))
  system(paste0("mkdir -p ", data_dir, "metal_slalom/"))
  
  for (h2_num in h2_vec) {
    
    # METAL result file does NOT depend on LD_BLOCK in the naming convention
    metal_results <- fread(
      paste0(result_dir, "Metal_", num_causal, "_h2_", h2_num, "_result1.tbl")
    )
    
    for (LD_BLOCK in LD_BLOCK_vec) {
      
      ##############################################
      #
      #        Format METAL Input for this locus
      #
      ##############################################
      
      zfile <- read.table(
        paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num),
        header = TRUE
      )
      
      metal_results_loci <- metal_results %>%
        filter(MarkerName %in% zfile$RSID) %>%
        separate(MarkerName, into = c("CHR","POS"),
                 sep = ":", remove = FALSE, convert = TRUE) %>%
        arrange(CHR, POS) %>%
        mutate(Z = qnorm(1 - `P-value` / 2) * sign(Effect))
      
      summary_stat_slalom <- metal_results_loci %>%
        dplyr::rename(
          chromosome = CHR,
          position   = POS,
          rsid       = MarkerName,
          allele1    = Allele1,
          allele2    = Allele2,
          beta       = Effect,
          se         = StdErr,
          p          = `P-value`
        ) %>%
        mutate(
          allele1 = toupper(allele1),
          allele2 = toupper(allele2),
          n_nfe   = 300000,
          n_afr   = 300000
        )
      
      out_file <- paste0(
        data_dir, "metal_slalom/CAUSAL_",
        num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_metal.snp"
      )
      
      write.table(
        summary_stat_slalom,
        file      = out_file,
        col.names = TRUE,
        row.names = FALSE,
        quote     = FALSE,
        sep       = "\t"
      )
    }
  }
}
