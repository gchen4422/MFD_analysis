library(data.table)
library(dplyr)

eur_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_eur"
afr_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_afr"

loci_num <- c(1:100)


for (loci in loci_num) {
  command_eur = paste0("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_eur/loci_1kg_",loci," -freq -out /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_eur/alt_freqs_loci_1kg_",loci)
  system(command_eur)
  command_afr = paste0("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_afr/loci_1kg_",loci," -freq -out /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/risk_loci_ld_afr/alt_freqs_loci_1kg_",loci)
  system(command_afr)
}







# helper to pick the ALT freq column even if named slightly differently
get_alt_freq <- function(dt) {
  cand <- intersect(c("ALT_FREQS","ALT_FRQ","ALT_FREQ"), names(dt))
  if (length(cand) == 0) stop("No ALT frequency column found.")
  dt[[cand[1]]]
}

for (loci in loci_num) {
  eu_file <- file.path(eur_dir, sprintf("alt_freqs_loci_1kg_%d.afreq", loci))
  bb_file <- file.path(afr_dir, sprintf("alt_freqs_loci_1kg_%d.afreq", loci))
  
  if (!file.exists(eu_file)) { message("Missing: ", eu_file, " — skipping"); next }
  if (!file.exists(bb_file)) { message("Missing: ", bb_file, " — skipping"); next }
  
  maf_EU <- fread(eu_file)
  maf_BB <- fread(bb_file)
  
  maf_EU_vec <- get_alt_freq(maf_EU)
  maf_BB_vec <- get_alt_freq(maf_BB)
  
  maf_EU_vec_output <- file.path(eur_dir, sprintf("alt_freqs_loci_1kg_%d.txt", loci))
  maf_BB_vec_output <- file.path(afr_dir, sprintf("alt_freqs_loci_1kg_%d.txt", loci))
  
  fwrite(data.table(maf_EU_vec), maf_EU_vec_output, col.names = FALSE, sep = "\t")
  fwrite(data.table(maf_BB_vec), maf_BB_vec_output, col.names = FALSE, sep = "\t")
}
