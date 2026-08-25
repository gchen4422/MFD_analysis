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
all_files <- list.files(path = "../", pattern = "_regenie_results_From5_20250318_combined_format_noMHC.csv.gz", full.names = TRUE)

# Identify AFR and EUR files
afr_files <- all_files[grep("_afr_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]



for (afr_file in afr_files) {
  
  trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  # eur_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/processed/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  # afr_ref = fread(paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/processed/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped.bim"))
  # chr = 0
  
  geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_ld_eur_",trait,"/")
  LD_dir_eur = geno_dir_eur
  
  geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_ld_afr_",trait,"/")
  LD_dir_afr = geno_dir_afr
  
  #annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_annot_",trait,"/")
  
  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))
  
  
  #for(i in 1:20){
  for(i in 21:nrow(candidate_region)){
    #for(i in 1:5){
    #region_chr = candidate_region$CHR[i]
    #region_start = candidate_region$MinBP[i]
    #region_end = candidate_region$MaxBP[i]
    
    #sumstats_afr <- afr_sumstats_sig %>%
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
      
      gc()
      
    }
    
    ### For MESuSiE, retain SNPs shared by the AFR and EUR summary statistics and LD matrices.
    ###Can start the testing for marginal ancestry-specific signals from here
    
    
    
    
    fwrite(sumstats_afr,file = paste0(data_dir,"afr_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    fwrite(sumstats_eur,file = paste0(data_dir,"eur_loci_",i,"_summ"), quote = FALSE, row.names=FALSE, col.names = T)
    
    fwrite(BB_cov, paste0(data_dir,"afr_loci_",i,".ld"),sep =" ")
    fwrite(EU_cov, paste0(data_dir,"eur_loci_",i,".ld"),sep =" ")
    
    #annot_list = fread(paste0(annot_dir,"loci_",i,".annot"))
    #annot_filtered = annot_list %>% filter(SNP %in% snplist_filtered)
    #write.table(annot_filtered,file = paste0(data_dir,"loci_",i,".annot"), quote = FALSE, row.names=FALSE, col.names = T)
  }
  
  message(paste("Processed files for trait:", trait))
  
}


#####################
#Get the Omega matrix for XMAP
#####################

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
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  result_dir<-paste0(wrk_dir,"result/")
  out_dir<-paste0(wrk_dir,"out/")
  system(paste0("mkdir -p ",data_dir))
  system(paste0("mkdir -p ",result_dir))
  system(paste0("mkdir -p ",out_dir))
  
  afr_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  eur_sumstats_sig = fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"))
  # pre-process: overlap SNPs
  snps <- Reduce(intersect, list(ldscore$rsid, afr_sumstats_sig$SNP, eur_sumstats_sig$SNP))
  sumstat_AFR_ldsc <- afr_sumstats_sig[match(snps, afr_sumstats_sig$SNP),]
  sumstat_EUR_ldsc <- eur_sumstats_sig[match(snps, eur_sumstats_sig$SNP),]
  ldscore_sub <- ldscore[match(snps, ldscore$rsid),]
  # pre-process: flip alleles
  
  z_afr <- sumstat_AFR_ldsc$Z
  z_eur <- sumstat_EUR_ldsc$Z
  
  idx_flip <- which(sumstat_AFR_ldsc$ALT != ldscore_sub$allele1 & sumstat_AFR_ldsc$ALT != comple(ldscore_sub$allele1))
  z_afr[idx_flip] <- -z_afr[idx_flip]
  
  idx_flip <- which(sumstat_EUR_ldsc$ALT != ldscore_sub$allele1 & sumstat_EUR_ldsc$ALT != comple(ldscore_sub$allele1))
  z_eur[idx_flip] <- -z_eur[idx_flip]
  
  
  idx1 <- which(z_afr^2 < 30 & z_eur^2 < 30)
  ld_afr_w <- 1 / sapply(ldscore_sub$AFR, function(x) max(x, 1))
  ld_eur_w <- 1 / sapply(ldscore_sub$EUR, function(x) max(x, 1))
  
  # bi-variate LDSC: AFR-EUR
  # stage 1 of LDSC: estimate intercepts
  fit_step1 <- estimate_gc(data.frame(Z = z_afr[idx1], N = sumstat_AFR_ldsc$N[idx1]), data.frame(Z = z_eur[idx1], N = sumstat_EUR_ldsc$N[idx1]),
                           ldscore_sub$AFR[idx1], ldscore_sub$EUR[idx1], ldscore_sub$AFR_EUR[idx1],
                           reg_w1 = ld_afr_w[idx1], reg_w2 = ld_eur_w[idx1], reg_wx = sqrt(ld_afr_w[idx1] * ld_eur_w[idx1]),
                           constrain_intercept = F)
  # stage 2 of LDSC: fix intercepts and estimate slopes
  fit_step2 <- estimate_gc(data.frame(Z = z_afr, N = sumstat_AFR_ldsc$N), data.frame(Z = z_eur, N = sumstat_EUR_ldsc$N),
                           ldscore_sub$AFR, ldscore_sub$EUR, ldscore_sub$AFR_EUR,
                           reg_w1 = ld_afr_w, reg_w2 = ld_eur_w, reg_wx = sqrt(ld_afr_w * ld_eur_w),
                           constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])
  
  # Assign LDSC estimates to covariance of polygenic effects
  OmegaHat <- diag(c(fit_step2$tau1$coefs[2], fit_step2$tau2$coefs[2])) # AFR EUR
  OmegaHat[1, 2] <- fit_step2$theta$coefs[2] # co-heritability
  OmegaHat[lower.tri(OmegaHat)] <- OmegaHat[upper.tri(OmegaHat)]
  
  c1 <- fit_step2$tau1$coefs[1] # AFR
  c2 <- fit_step2$tau2$coefs[1] # EUR
  
  save(c1,c2,OmegaHat,file = paste0(data_dir,trait,"_xmap_omega.Rdata"))  
}

# for (afr_file in afr_files) {
#   trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
#   candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/summstats_to_finemap/", trait, "_risk_loci"))
#   if(nrow(candidate_region)<10){
#     message(paste("Processed files for trait:", trait))
#   }
# }
# 
# Processed files for trait: DEP
# Processed files for trait: SOD



#LD1 = fread("risk_loci_ld_eur_DBP/loci_1.ld")



#if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

#LD_plot <- if (max(LD1, na.rm = TRUE) <= 1 && min(LD1, na.rm = TRUE) >= -1) LD1^2 else LD1

#pheatmap::pheatmap(
#  LD_plot,
#  color = colorRampPalette(c("white", "red"))(200),
#  cluster_rows = FALSE, cluster_cols = FALSE,           # no clustering
#  show_rownames = FALSE, show_colnames = FALSE,         # hide SNP names
#  treeheight_row = 0, treeheight_col = 0,               # remove dendrogram space
#  legend = TRUE,
#  display_numbers = FALSE,
#  border_color = NA
#)
