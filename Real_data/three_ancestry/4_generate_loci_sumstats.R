library(data.table)
library(susieR)
library(snpStats)
library(dplyr)

# args <- as.numeric(commandArgs(TRUE))
# i <- args[1]
# h2_num <- args[2]
# num_causal <- args[3]


preprocess_data <- function(data_matrix) {
  # Replace missing values with column means and scale columns
  data_matrix <- apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x)
    return(x)
  })

  return(data_matrix)
}




# List all files
all_files <- list.files(path = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/formated_summstats/", pattern = "_format_noMHC.csv.gz", full.names = TRUE)

# Identify eas and EUR files
eas_files <- all_files[grep("_EAS_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]
afr_files <- all_files[grep("_afr_", all_files)]


for (eas_file in eas_files) {

  trait <- sub("^[^_]+_([^_]+)_.*", "\\1", basename(eas_file))
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_3anc_loci_not_copied.txt"))
  colnames(candidate_region)[2] = "CHR"
  colnames(candidate_region)[3] = "MinBP"
  colnames(candidate_region)[4] = "MaxBP"
  candidate_region$loci_id <- as.numeric(gsub("loci_", "", candidate_region$`3anc_loci_id`))
  #eas_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF_ukb_tpmi/1kg_eas_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #eur_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #afr_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #chr = 0

  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_eur_",trait,"/")
  LD_dir_eur = geno_dir_eur

  geno_dir_eas <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_eas_",trait,"/")
  LD_dir_eas = geno_dir_eas

  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_afr_",trait,"/")
  LD_dir_afr = geno_dir_afr

  #annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_annot_",trait,"/")

  eas_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"))

  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))


  #for(i in 1:20){
  for(j in 1:nrow(candidate_region)){
    i <- candidate_region$loci_id[j]
    #for(i in 1:5){
    #region_chr = candidate_region$CHR[i]
    #region_start = candidate_region$MinBP[i]
    #region_end = candidate_region$MaxBP[i]

    #sumstats_eas <- eas_sumstats_sig %>%
    #  filter(
    #    CHR == region_chr,
    #    POS >= region_start,
    #    POS <= region_end
    #  )

    #sumstats_eur <- eur_sumstats_sig %>%
    #  filter(
    #    CHR == region_chr,
    #    POS >= region_start,
    #    POS <= region_end
    #  )



    snplist_ori_eas = fread(paste0(geno_dir_eas,"loci_",i),header = F)%>% pull(V1)
    snp_to_exclude_path_eas <- paste0(geno_dir_eas, "loci_", i, "_snp_to_exclude")

    # Check if the file is empty
    if (file.info(snp_to_exclude_path_eas)$size > 0) {
      # Read the file if it's not empty
      snp_to_exclude_eas <- fread(snp_to_exclude_path_eas, header = FALSE) %>% pull(V1)
      snplist_filtered_eas <- setdiff(snplist_ori_eas, snp_to_exclude_eas)
    } else {
      # Handle the case where the file is empty
      snplist_filtered_eas  = snplist_ori_eas
    }


    ##Check for LD mismatch, remove mismatched SNPs
    idx_to_remove = c()
    while(!(identical(idx_to_remove, integer(0)))){



      sumstats_eas = eas_sumstats_sig %>% filter(SNP %in% snplist_filtered_eas)



      eas_plink <- read.plink(paste0(geno_dir_eas,"loci_",i,".bed"))
      eas_plink_geno <- as(eas_plink$genotypes, "numeric")
      #unique_value_cols_eas <- apply(eas_plink_geno, 2, function(col) length(unique(col)) == 1)
      #unique_value_col_indices_eas <- which(unique_value_cols_eas)
      eas_mat = as.matrix(eas_plink_geno)
      eas_mat <- eas_mat[, colnames(eas_mat) %in% snplist_filtered_eas]
      eas_mat = preprocess_data(eas_mat)
      cov_matrix_eas <- cov2cor(crossprod(eas_mat))
      EA_cov <- (cov_matrix_eas + t(cov_matrix_eas)) / 2

      EA_diagnostic <- kriging_rss(sumstats_eas$Z, EA_cov, n = median(sumstats_eas$N))
      idx_to_remove_eas = EA_diagnostic$plot$plot_env$idx
      idx_to_remove = idx_to_remove_eas

      if (!(identical(idx_to_remove, integer(0)))){
        snplist_filtered_eas = snplist_filtered_eas[-idx_to_remove]
      }

      #gc()

    }

    snplist_ori_eur = fread(paste0(geno_dir_eur,"loci_",i),header = F)%>% pull(V1)
    snp_to_exclude_path_eur <- paste0(geno_dir_eur, "loci_", i, "_snp_to_exclude")

    # Check if the file is empty
    if (file.info(snp_to_exclude_path_eur)$size > 0) {
      # Read the file if it's not empty
      snp_to_exclude_eur <- fread(snp_to_exclude_path_eur, header = FALSE) %>% pull(V1)
      snplist_filtered_eur <- setdiff(snplist_ori_eur, snp_to_exclude_eur)
    } else {
      # Handle the case where the file is empty
      snplist_filtered_eur  = snplist_ori_eur
    }


    ##Check for LD mismatch, remove mismatched SNPs
    idx_to_remove = c()
    while(!(identical(idx_to_remove, integer(0)))){



      sumstats_eur = eur_sumstats_sig %>% filter(SNP %in% snplist_filtered_eur)



      eur_plink <- read.plink(paste0(geno_dir_eur,"loci_",i,".bed"))
      eur_plink_geno <- as(eur_plink$genotypes, "numeric")
      #unique_value_cols_eur <- apply(eur_plink_geno, 2, function(col) length(unique(col)) == 1)
      #unique_value_col_indices_eur <- which(unique_value_cols_eur)
      eur_mat = as.matrix(eur_plink_geno)
      eur_mat <- eur_mat[, colnames(eur_mat) %in% snplist_filtered_eur]
      eur_mat = preprocess_data(eur_mat)
      cov_matrix_eur <- cov2cor(crossprod(eur_mat))
      EU_cov <- (cov_matrix_eur + t(cov_matrix_eur)) / 2

      EU_diagnostic <- kriging_rss(sumstats_eur$Z, EU_cov, n = median(sumstats_eur$N))
      idx_to_remove_eur = EU_diagnostic$plot$plot_env$idx
      idx_to_remove = idx_to_remove_eur

      if (!(identical(idx_to_remove, integer(0)))){
        snplist_filtered_eur = snplist_filtered_eur[-idx_to_remove]
      }
    }

    snplist_ori_afr = fread(paste0(geno_dir_afr,"loci_",i),header = F)%>% pull(V1)
    snp_to_exclude_path_afr <- paste0(geno_dir_afr, "loci_", i, "_snp_to_exclude")

    # Check if the file is empty
    if (file.info(snp_to_exclude_path_afr)$size > 0) {
      # Read the file if it's not empty
      snp_to_exclude_afr <- fread(snp_to_exclude_path_afr, header = FALSE) %>% pull(V1)
      snplist_filtered_afr <- setdiff(snplist_ori_afr, snp_to_exclude_afr)
    } else {
      # Handle the case where the file is empty
      snplist_filtered_afr  = snplist_ori_afr
    }


    ##Check for LD mismatch, remove mismatched SNPs
    idx_to_remove = c()
    while(!(identical(idx_to_remove, integer(0)))){



      sumstats_afr = afr_sumstats_sig %>% filter(SNP %in% snplist_filtered_afr)



      afr_plink <- read.plink(paste0(geno_dir_afr,"loci_",i,".bed"))
      afr_plink_geno <- as(afr_plink$genotypes, "numeric")
      #unique_value_cols_afr <- apply(afr_plink_geno, 2, function(col) length(unique(col)) == 1)
      #unique_value_col_indices_afr <- which(unique_value_cols_afr)
      afr_mat = as.matrix(afr_plink_geno)
      afr_mat <- afr_mat[, colnames(afr_mat) %in% snplist_filtered_afr]
      afr_mat = preprocess_data(afr_mat)
      cov_matrix_afr <- cov2cor(crossprod(afr_mat))
      BB_cov <- (cov_matrix_afr + t(cov_matrix_afr)) / 2

      BB_diagnostic <- kriging_rss(sumstats_afr$Z, BB_cov, n = median(sumstats_afr$N))
      idx_to_remove_afr = BB_diagnostic$plot$plot_env$idx
      idx_to_remove = idx_to_remove_afr

      if (!(identical(idx_to_remove, integer(0)))){
        snplist_filtered_afr = snplist_filtered_afr[-idx_to_remove]
      }

      #gc()

    }

      gc()



    ### For MESuSiE, retain SNPs shared by the EAS and EUR summary statistics and LD matrices.
    ###Can start the testing for marginal ancestry-specific signals from here




    fwrite(sumstats_eas,file = paste0(data_dir,"eas_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    fwrite(sumstats_eur,file = paste0(data_dir,"eur_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    fwrite(sumstats_afr, file = paste0(data_dir,"afr_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)


    fwrite(EA_cov, paste0(data_dir,"eas_loci_",i,".ld"),sep =" ")
    fwrite(EU_cov, paste0(data_dir,"eur_loci_",i,".ld"),sep =" ")
    fwrite(BB_cov, paste0(data_dir,"afr_loci_",i,".ld"), sep = " ")
  }


  message(paste("Processed files for trait:", trait))

}






for (eas_file in eas_files) {

  trait <- sub("^[^_]+_([^_]+)_.*", "\\1", basename(eas_file))
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  #eas_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF_ukb_tpmi/1kg_eas_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #eur_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #afr_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  #chr = 0


  geno_dir_eas <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_eas_",trait,"/")
  LD_dir_eas = geno_dir_eas



  #annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_annot_",trait,"/")

  eas_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap.csv"))

  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))


  #for(i in 1:20){
  for(i in 1:nrow(candidate_region)){

    if(file.exists(paste0(data_dir,"eas_loci_",i,".ld"))) next
    #for(i in 1:5){
    #region_chr = candidate_region$CHR[i]
    #region_start = candidate_region$MinBP[i]
    #region_end = candidate_region$MaxBP[i]

    #sumstats_eas <- eas_sumstats_sig %>%
    #  filter(
    #    CHR == region_chr,
    #    POS >= region_start,
    #    POS <= region_end
    #  )

    #sumstats_eur <- eur_sumstats_sig %>%
    #  filter(
    #    CHR == region_chr,
    #    POS >= region_start,
    #    POS <= region_end
    #  )



    snplist_ori_eas = fread(paste0(geno_dir_eas,"loci_",i),header = F)%>% pull(V1)
    snp_to_exclude_path_eas <- paste0(geno_dir_eas, "loci_", i, "_snp_to_exclude")

    # Check if the file is empty
    if (file.info(snp_to_exclude_path_eas)$size > 0) {
      # Read the file if it's not empty
      snp_to_exclude_eas <- fread(snp_to_exclude_path_eas, header = FALSE) %>% pull(V1)
      snplist_filtered_eas <- setdiff(snplist_ori_eas, snp_to_exclude_eas)
    } else {
      # Handle the case where the file is empty
      snplist_filtered_eas  = snplist_ori_eas
    }


    ##Check for LD mismatch, remove mismatched SNPs
    idx_to_remove = c()
    while(!(identical(idx_to_remove, integer(0)))){



      sumstats_eas = eas_sumstats_sig %>% filter(SNP %in% snplist_filtered_eas)



      eas_plink <- read.plink(paste0(geno_dir_eas,"loci_",i,".bed"))
      eas_plink_geno <- as(eas_plink$genotypes, "numeric")
      #unique_value_cols_eas <- apply(eas_plink_geno, 2, function(col) length(unique(col)) == 1)
      #unique_value_col_indices_eas <- which(unique_value_cols_eas)
      eas_mat = as.matrix(eas_plink_geno)
      eas_mat <- eas_mat[, colnames(eas_mat) %in% snplist_filtered_eas]
      eas_mat = preprocess_data(eas_mat)
      cov_matrix_eas <- cov2cor(crossprod(eas_mat))
      EA_cov <- (cov_matrix_eas + t(cov_matrix_eas)) / 2

      EA_diagnostic <- kriging_rss(sumstats_eas$Z, EA_cov, n = median(sumstats_eas$N))
      idx_to_remove_eas = EA_diagnostic$plot$plot_env$idx
      idx_to_remove = idx_to_remove_eas

      if (!(identical(idx_to_remove, integer(0)))){
        snplist_filtered_eas = snplist_filtered_eas[-idx_to_remove]
      }

      #gc()

    }





  fwrite(sumstats_eas,file = paste0(data_dir,"eas_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)

  fwrite(EA_cov, paste0(data_dir,"eas_loci_",i,".ld"),sep =" ")

  }
  message(paste("Processed files for trait:", trait))

}


# for (eas_file in eas_files) {
#   trait <- gsub(".*/(.*?)_.*", "\\1", eas_file)
#   candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/summstats_to_finemap/", trait, "_risk_loci"))
#   if(nrow(candidate_region)<10){
#     message(paste("Processed files for trait:", trait))
#   }
# }
#
# Processed files for trait: DEP
# Processed files for trait: SOD



#####################
#Get the Omega matrix for XMAP
#####################
library(XMAP)
ldscore <- data.frame()
for (chr in 1:22) {
  ldscore_chr <- fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/LD_score/LDscore_eas_brit_afr_chr", chr, ".txt"))
  ldscore <- rbind(ldscore, ldscore_chr)
  cat("CHR", chr, "\n")
}

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]


traits = c("BMI","DBP","SBP")
for (trait in traits) {

  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))

  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eas_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap.csv"))

  # pre-process: overlap SNPs
  snps <- Reduce(intersect, list(ldscore$rsid, afr_sumstats_sig$SNP, eur_sumstats_sig$SNP,eas_sumstats_sig$SNP))
  sumstat_AFR_ldsc <- afr_sumstats_sig[match(snps, afr_sumstats_sig$SNP),]
  sumstat_EUR_ldsc <- eur_sumstats_sig[match(snps, eur_sumstats_sig$SNP),]
  sumstat_EAS_ldsc <- eas_sumstats_sig[match(snps, eas_sumstats_sig$SNP),]
  ldscore_sub <- ldscore[match(snps, ldscore$rsid),]
  # pre-process: flip alleles

  z_afr <- sumstat_AFR_ldsc$Z
  z_eur <- sumstat_EUR_ldsc$Z
  z_eas <- sumstat_EAS_ldsc$Z

  idx_flip <- which(sumstat_AFR_ldsc$ALT != ldscore_sub$allele1 & sumstat_AFR_ldsc$ALT != comple(ldscore_sub$allele1))
  z_afr[idx_flip] <- -z_afr[idx_flip]

  idx_flip <- which(sumstat_EUR_ldsc$ALT != ldscore_sub$allele1 & sumstat_EUR_ldsc$ALT != comple(ldscore_sub$allele1))
  z_eur[idx_flip] <- -z_eur[idx_flip]

  idx_flip <- which(sumstat_EAS_ldsc$ALT != ldscore_sub$allele1 & sumstat_EAS_ldsc$ALT != comple(ldscore_sub$allele1))
  z_eas[idx_flip] <- -z_eas[idx_flip]

  idx1 <- which(z_afr^2 < 30 & z_eur^2 < 30 & z_eas^2 < 30)
  ld_afr_w <- 1 / sapply(ldscore_sub$AFR, function(x) max(x, 1))
  ld_eur_w <- 1 / sapply(ldscore_sub$EUR, function(x) max(x, 1))
  ld_eas_w <- 1 / sapply(ldscore_sub$EAS, function(x) max(x, 1))

  # ---- Pair 1: AFR-EUR ----
  fit_afr_eur_s1 <- estimate_gc(
    data.frame(Z = z_afr[idx1], N = sumstat_AFR_ldsc$N[idx1]),
    data.frame(Z = z_eur[idx1], N = sumstat_EUR_ldsc$N[idx1]),
    ldscore_sub$AFR[idx1], ldscore_sub$EUR[idx1], ldscore_sub$AFR_EUR[idx1],
    reg_w1 = ld_afr_w[idx1], reg_w2 = ld_eur_w[idx1],
    reg_wx = sqrt(ld_afr_w[idx1] * ld_eur_w[idx1]),
    constrain_intercept = F)
  fit_afr_eur_s2 <- estimate_gc(
    data.frame(Z = z_afr, N = sumstat_AFR_ldsc$N),
    data.frame(Z = z_eur, N = sumstat_EUR_ldsc$N),
    ldscore_sub$AFR, ldscore_sub$EUR, ldscore_sub$AFR_EUR,
    reg_w1 = ld_afr_w, reg_w2 = ld_eur_w,
    reg_wx = sqrt(ld_afr_w * ld_eur_w),
    constrain_intercept = T,
    fit_afr_eur_s1$tau1$coefs[1], fit_afr_eur_s1$tau2$coefs[1], fit_afr_eur_s1$theta$coefs[1])

  # ---- Pair 2: AFR-EAS ----
  fit_afr_eas_s1 <- estimate_gc(
    data.frame(Z = z_afr[idx1], N = sumstat_AFR_ldsc$N[idx1]),
    data.frame(Z = z_eas[idx1], N = sumstat_EAS_ldsc$N[idx1]),
    ldscore_sub$AFR[idx1], ldscore_sub$EAS[idx1], ldscore_sub$AFR_EAS[idx1],
    reg_w1 = ld_afr_w[idx1], reg_w2 = ld_eas_w[idx1],
    reg_wx = sqrt(ld_afr_w[idx1] * ld_eas_w[idx1]),
    constrain_intercept = F)
  fit_afr_eas_s2 <- estimate_gc(
    data.frame(Z = z_afr, N = sumstat_AFR_ldsc$N),
    data.frame(Z = z_eas, N = sumstat_EAS_ldsc$N),
    ldscore_sub$AFR, ldscore_sub$EAS, ldscore_sub$AFR_EAS,
    reg_w1 = ld_afr_w, reg_w2 = ld_eas_w,
    reg_wx = sqrt(ld_afr_w * ld_eas_w),
    constrain_intercept = T,
    fit_afr_eas_s1$tau1$coefs[1], fit_afr_eas_s1$tau2$coefs[1], fit_afr_eas_s1$theta$coefs[1])

  # ---- Pair 3: EUR-EAS ----
  fit_eur_eas_s1 <- estimate_gc(
    data.frame(Z = z_eur[idx1], N = sumstat_EUR_ldsc$N[idx1]),
    data.frame(Z = z_eas[idx1], N = sumstat_EAS_ldsc$N[idx1]),
    ldscore_sub$EUR[idx1], ldscore_sub$EAS[idx1], ldscore_sub$EUR_EAS[idx1],
    reg_w1 = ld_eur_w[idx1], reg_w2 = ld_eas_w[idx1],
    reg_wx = sqrt(ld_eur_w[idx1] * ld_eas_w[idx1]),
    constrain_intercept = F)
  fit_eur_eas_s2 <- estimate_gc(
    data.frame(Z = z_eur, N = sumstat_EUR_ldsc$N),
    data.frame(Z = z_eas, N = sumstat_EAS_ldsc$N),
    ldscore_sub$EUR, ldscore_sub$EAS, ldscore_sub$EUR_EAS,
    reg_w1 = ld_eur_w, reg_w2 = ld_eas_w,
    reg_wx = sqrt(ld_eur_w * ld_eas_w),
    constrain_intercept = T,
    fit_eur_eas_s1$tau1$coefs[1], fit_eur_eas_s1$tau2$coefs[1], fit_eur_eas_s1$theta$coefs[1])

  # ---- Build 3x3 OmegaHat ----
  # Diagonal: heritability slopes from each pair (use the estimate where that ancestry is pop1)
  # AFR h2 from AFR-EUR pair, EUR h2 from AFR-EUR pair, EAS h2 from EUR-EAS pair
  h2_afr <- fit_afr_eur_s2$tau1$coefs[2]
  h2_eur <- fit_afr_eur_s2$tau2$coefs[2]
  h2_eas <- fit_eur_eas_s2$tau2$coefs[2]

  OmegaHat <- matrix(0, 3, 3)
  diag(OmegaHat) <- c(h2_afr, h2_eur, h2_eas)

  # Off-diagonal: co-heritability slopes
  OmegaHat[1, 2] <- fit_afr_eur_s2$theta$coefs[2]  # AFR-EUR
  OmegaHat[1, 3] <- fit_afr_eas_s2$theta$coefs[2]  # AFR-EAS
  OmegaHat[2, 3] <- fit_eur_eas_s2$theta$coefs[2]  # EUR-EAS
  OmegaHat[lower.tri(OmegaHat)] <- OmegaHat[upper.tri(OmegaHat)]

  # Intercepts (residual variance inflation)
  c1 <- fit_afr_eur_s2$tau1$coefs[1]  # AFR
  c2 <- fit_afr_eur_s2$tau2$coefs[1]  # EUR
  c3 <- fit_eur_eas_s2$tau2$coefs[1]  # EAS

  save(c1, c2, c3, OmegaHat, file = paste0(data_dir, trait, "_xmap_omega.Rdata"))
}
