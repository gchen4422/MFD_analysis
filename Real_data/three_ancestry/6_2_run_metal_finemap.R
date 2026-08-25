library(dplyr)
library(data.table)
library(susieR)
library(tidyr)


traits_list <- c("BMI", "DBP", "SBP")




for (trait in traits_list) {
  # Directories
  wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/")
  data_dir <- paste0(wrk_dir, "summary_data/")
  result_dir <- paste0(wrk_dir, "result/")
  out_dir <- paste0(wrk_dir, "out/")

  # Read Candidate Regions
  candidate_region_file <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_risk_loci")
  candidate_region <- fread(candidate_region_file)


  metal_result = fread(paste0(result_dir,"Metal_",trait,"_result1.tbl"))
  metal_result = metal_result %>%  separate(MarkerName, into = c("CHR","POS"), sep = "_", remove = FALSE, convert = TRUE)


  for(current_block_id in 1:nrow(candidate_region)){

    chr = candidate_region$seqnames[current_block_id]
    bp_start = candidate_region$start[current_block_id]
    bp_end = candidate_region$end[current_block_id]

    metal_results_loci_raw <- metal_result %>%
      filter(CHR == chr) %>%
      filter(POS >= bp_start) %>%
      filter(POS <= bp_end) %>%

      mutate(
        Allele1 = toupper(Allele1),
        Allele2 = toupper(Allele2)
      )  %>% arrange(CHR, POS) %>%
      mutate(Z = qnorm(log(`P-value`/2), lower.tail = FALSE, log.p = TRUE) * sign(Effect))

    load(paste0(result_dir,"MF_result_LOCI_",current_block_id,".RData"))

    mfd_result_all = mfd_result_all %>% mutate(CHR_POS = paste(CHR,POS,sep = "_"))
    metal_results_loci = metal_results_loci_raw %>% filter(MarkerName %in% mfd_result_all$CHR_POS)

    susie_metal<-susie_rss(metal_results_loci$Z,R = diag(length(metal_results_loci$Z)), n = median(metal_results_loci$TotalN), L = 1)


    if(!all(mfd_result_all$CHR_POS == metal_results_loci$MarkerName)) {
      warning("Mismatch in Credible Set variants after merge! Check SNP IDs error.")
    }

    susie_metal_cs<-rep(0,length(susie_metal$pip))
    susie_metal_cs[unlist(susie_metal$sets$cs)]<-1

    mfd_result_all = mfd_result_all %>% mutate(SuSiE_merged_PIP_Either = susie_metal$pip,SuSiE_metal_cs = susie_metal_cs) %>% select(-CHR_POS)
    mf_output_name<-paste0(result_dir,"MF_result_LOCI_merged_",current_block_id,".RData")
    save(mfd_result_all,file =mf_output_name )

  }
}
