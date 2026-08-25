library(dplyr)
library(data.table)
library(MESuSiE)
library(susieR)

args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]
causal_index = args[4]
External_index = args[5]
causal_index_name = c("Both","One")[causal_index]
External_index_name = c("","External_")[External_index]

library(data.table)
library(dplyr)
num_BB = 300000
num_EU = 300000

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
paintor_dir<-data_dir


zfile<-read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
EU_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
rownames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID
rownames(BB_cov)<-zfile$RSID


if(External_index ==1 ){
  num_snp_missing = round(nrow(zfile)*0.3)
}else if(External_index ==2){
  num_snp_missing = round(nrow(zfile)*0.3+nrow(zfile)*(1-0.3)*0.3)
}


non_causal_SNP = zfile%>%filter(Signal==0)%>%pull(RSID)
EU_missing_SNP = c(sample(non_causal_SNP,num_snp_missing-1),zfile%>%filter(Signal!=0)%>%pull(RSID))
BB_missing_SNP = EU_missing_SNP
zfile<-zfile%>%mutate(EU_missing = ifelse(RSID%in%EU_missing_SNP,1,0),BB_missing = ifelse(RSID%in%BB_missing_SNP,1,0))


# --- Define the Decision Function ---
decide_finemapping_method <- function(sum_eur, sum_afr, ld_eur, ld_afr, p_thresh = 5e-8, r2_thresh = 0.6) {
  
  # 1. Pre-processing
  sum_eur[, SNP := as.character(SNP)]
  sum_afr[, SNP := as.character(SNP)]
  
  check_and_convert_cor <- function(mat) {
    if (any(abs(diag(mat) - 1) > 0.1)) cov2cor(mat) else mat
  }
  R_eur <- check_and_convert_cor(ld_eur)
  R_afr <- check_and_convert_cor(ld_afr)
  
  # 2. Identify Variant Sets
  shared_snps <- intersect(sum_eur$SNP, sum_afr$SNP)
  asv_eur <- setdiff(sum_eur$SNP, sum_afr$SNP)
  asv_afr <- setdiff(sum_afr$SNP, sum_eur$SNP)
  
  # 3. Check Condition 1: Are AS-Vs Significant?
  sig_asv_eur <- sum_eur[SNP %in% asv_eur & PVAL < p_thresh, SNP]
  sig_asv_afr <- sum_afr[SNP %in% asv_afr & PVAL < p_thresh, SNP]
  has_sig_asv <- length(sig_asv_eur) > 0 || length(sig_asv_afr) > 0
  
  if (!has_sig_asv) {
    return(list(method = "MESuSiE", reason_code = 1,reason = "AS-V not significant"))
  }
  
  # 4. Prepare for LD Checks
  sig_shared_eur <- sum_eur[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  sig_shared_afr <- sum_afr[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  all_sig_shared <- unique(c(sig_shared_eur, sig_shared_afr))
  
  if (length(all_sig_shared) == 0) { # This case is that no significant shared signals even present
    return(list(method = "SuSiE post-hoc",reason_code = 2,reason = "AS-V significant, but not in high LD with shared signals"))
  }
  
  # 5. Check High LD Logic
  high_ld_found <- FALSE
  high_ld_partners <- character()
  
  # EUR Check
  if (length(sig_asv_eur) > 0) {
    valid_asv <- intersect(sig_asv_eur, rownames(R_eur))
    valid_shared <- intersect(all_sig_shared, colnames(R_eur))
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_eur[valid_asv, valid_shared, drop=FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(high_ld_partners, colnames(ld_sub)[apply(ld_sub^2, 2, max) > r2_thresh])
      }
    }
  }
  
  # AFR Check
  if (length(sig_asv_afr) > 0) {
    valid_asv <- intersect(sig_asv_afr, rownames(R_afr))
    valid_shared <- intersect(all_sig_shared, colnames(R_afr))
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_afr[valid_asv, valid_shared, drop=FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(high_ld_partners, colnames(ld_sub)[apply(ld_sub^2, 2, max) > r2_thresh])
      }
    }
  }
  
  high_ld_partners <- unique(high_ld_partners)
  
  if (!high_ld_found) {
    return(list(method = "SuSiE post-hoc",reason_code = 2,reason = "AS-V significant, but not in high LD with shared signals"))
  }
  
  # 6. Check Condition 3: Nature of Shared SNPs
  p_eur_check <- sum_eur[SNP %in% high_ld_partners, PVAL]
  p_afr_check <- sum_afr[SNP %in% high_ld_partners, PVAL]
  is_sig_both <- (p_eur_check < p_thresh) & (p_afr_check < p_thresh)
  
  if (any(is_sig_both)) {
    return(list(method = "MESuSiE", reason_code = 3,reason = "AS-V significant and in high LD with a Shared SNP significant in ALL ancestries"))
  } else {
    return(list(method = "SuSiE post-hoc",reason_code = 4,reason = "AS-V significant and in high LD with a Shared SNP significant in ONLY specific ancestry"))
  }
}




# DETERMINE BLOCK ID:
current_block_id <- LD_BLOCK

# Use tryCatch to skip regions where files might be missing without stopping the whole script
skip_to_next <- FALSE

cat(paste0("Processing Block: ", current_block_id, "... "))


# --- Run Decision Logic ---

if(causal_index==1){
  susie_EU_zfile<-zfile%>%filter(EU_missing==0)
  susie_BB_zfile<-zfile%>%filter(BB_missing==0)
  missing_EU_index = which(zfile$EU_missing!=0)
  missing_BB_index = which(zfile$BB_missing!=0)
  susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
  susie_BB_cov<-BB_cov[-missing_BB_index,-missing_BB_index]
}else if(causal_index==2){
  susie_EU_zfile<-zfile%>%filter(EU_missing==0)
  susie_BB_zfile<-zfile
  missing_EU_index = which(zfile$EU_missing!=0)
  susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
  susie_BB_cov<-BB_cov
}


summary_stat_1 = data.frame("SNP" = susie_EU_zfile$RSID,"CHR" = susie_EU_zfile$CHR,"POS" = susie_EU_zfile$POS,"Signal" = susie_EU_zfile$Signal,"Beta"=susie_EU_zfile$zscore_1/sqrt(num_EU),"Se"=1/sqrt(num_EU), "Z" =susie_EU_zfile$zscore_1 ,  "N" =num_EU ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
summary_stat_2 = data.frame("SNP" = susie_BB_zfile$RSID,"CHR" = susie_BB_zfile$CHR,"POS" = susie_BB_zfile$POS,"Signal" = susie_BB_zfile$Signal,"Beta"=susie_BB_zfile$zscore_2/sqrt(num_BB),"Se"=1/sqrt(num_BB), "Z" =susie_BB_zfile$zscore_2 ,  "N" =num_BB ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))

setDT(summary_stat_1)
setDT(summary_stat_2)

decision <- decide_finemapping_method(summary_stat_1, summary_stat_2, susie_EU_cov, susie_BB_cov)


if(decision$method == "SuSiE post-hoc"){
  
  susie_EU<-susie_rss(summary_stat_1$Z,susie_EU_cov,check_prior=F)
  susie_BB<-susie_rss(summary_stat_2$Z,susie_BB_cov,check_prior=F)
  
  
  susie_BB_cs<-rep(0,length(susie_BB$pip))
  susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
  for(j in seq_along(susie_BB$sets$cs)) {
    susie_BB_cs[susie_BB$sets$cs[[j]]] <- j
  }
  susie_EU_cs<-rep(0,length(susie_EU$pip))
  susie_EU_cs[unlist(susie_EU$sets$cs)]<-1
  for(k in seq_along(susie_EU$sets$cs)) {
    susie_EU_cs[susie_EU$sets$cs[[k]]] <- k
  }
  
  SuSiE_EU_res_dataframe<-data.frame(SNP = summary_stat_1$SNP,SuSiE_PIP_EU = susie_EU$pip,susie_EU_cs = susie_EU_cs )
  SuSiE_BB_res_dataframe<-data.frame(SNP = summary_stat_2$SNP,SuSiE_PIP_BB = susie_BB$pip,susie_BB_cs = susie_BB_cs )
  
  mf_result_list <- summary_stat_1 %>% 
    select(SNP, CHR, POS) %>% 
    full_join(summary_stat_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
    mutate(
      CHR = coalesce(CHR.x, CHR.y),
      POS = coalesce(POS.x, POS.y)
    ) %>%
    select(SNP, CHR, POS) %>% 
    left_join(SuSiE_EU_res_dataframe, by = "SNP")%>%
    left_join(SuSiE_BB_res_dataframe, by = "SNP")%>% replace(is.na(.), 0)
  
  mf_result_list <- mf_result_list%>%mutate(SuSiE_cs = ifelse(susie_EU_cs+susie_BB_cs==0,0,1),SuSiE_PIP_Either = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),SuSiE_PIP_Shared = pmin(SuSiE_PIP_EU,SuSiE_PIP_BB)) %>% 
    select(SNP,CHR, POS,SuSiE_PIP_Either,SuSiE_PIP_Shared,SuSiE_PIP_EU,SuSiE_PIP_BB,SuSiE_cs)%>% arrange(CHR,POS) 
  
  
}else{
  
  if(causal_index==1){
    missing_index = which(zfile$EU_missing!=0|zfile$BB_missing!=0)
    zfile_subset<-zfile%>%filter(EU_missing==0,BB_missing==0)
    EU_cov_subset<-EU_cov[-missing_index,-missing_index]
    BB_cov_subset<-BB_cov[-missing_index,-missing_index]
  }else if(causal_index==2){
    missing_index = which(zfile$EU_missing!=0)
    zfile_subset<-zfile%>%filter(EU_missing==0)
    EU_cov_subset<-EU_cov[-missing_index,-missing_index]
    BB_cov_subset<-BB_cov[-missing_index,-missing_index]
  }
  
  summary_stat_1_common = data.frame("SNP" = zfile_subset$RSID,"CHR" = zfile_subset$CHR,"POS" = zfile_subset$POS,"Signal" = zfile_subset$Signal,"Beta"=zfile_subset$zscore_1/sqrt(num_EU),"Se"=1/sqrt(num_EU), "Z" =zfile_subset$zscore_1 ,  "N" =num_EU ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
  summary_stat_2_common = data.frame("SNP" = zfile_subset$RSID,"CHR" = zfile_subset$CHR,"POS" = zfile_subset$POS,"Signal" = zfile_subset$Signal,"Beta"=zfile_subset$zscore_2/sqrt(num_BB),"Se"=1/sqrt(num_BB), "Z" =zfile_subset$zscore_2 ,  "N" =num_BB ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
  summary_stat_sd_list = list("EUR" = summary_stat_1_common,"AFR"=summary_stat_2_common )  
  R_mat_list=list("EUR" = EU_cov_subset,"AFR" = BB_cov_subset)
  
  
  MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10)
  
  
  MESuSiE_cs<-rep(0,length(MESuSiE_res$pip))
  for(m in seq_along(MESuSiE_res$cs$cs)) {
    MESuSiE_cs[MESuSiE_res$cs$cs[[m]]] <- m
  }
  MESuSiE_res_dataframe<-data.frame(SNP = summary_stat_1_common$SNP,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )
  mf_result_list <- summary_stat_1 %>% 
    select(SNP, CHR, POS) %>% 
    full_join(summary_stat_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
    mutate(
      CHR = coalesce(CHR.x, CHR.y),
      POS = coalesce(POS.x, POS.y)
    ) %>%
    select(SNP, CHR, POS) %>% 
    left_join(MESuSiE_res_dataframe, by = "SNP") %>% 
    select(SNP,CHR, POS,MESuSiE_PIP_Either, MESuSiE_PIP_Shared, MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_cs) %>% arrange(CHR,POS) %>% replace(is.na(.), 0) 
  
}






susie_EU$sets$cs
susie_BB$sets$cs

shared_snps <- intersect(summary_stat_eur$SNP, summary_stat_afr$SNP)
summary_stat_eur_common = summary_stat_eur %>% filter(SNP %in% shared_snps)
summary_stat_afr_common = summary_stat_afr %>% filter(SNP %in% shared_snps)
EU_cov_common = EU_cov[shared_snps,shared_snps]
BB_cov_common = BB_cov[shared_snps,shared_snps]

summary_stat_sd_list = list("EUR" = summary_stat_eur_common,"AFR"=summary_stat_afr_common ) 
R_mat_list=list("EUR" = EU_cov_common,"AFR" = BB_cov_common)
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10)


MESuSiE_cs<-rep(0,length(MESuSiE_res$pip))
for(m in seq_along(MESuSiE_res$cs$cs)) {
  MESuSiE_cs[MESuSiE_res$cs$cs[[m]]] <- m
}
MESuSiE_res_dataframe<-data.frame(SNP = summary_stat_1_common$SNP,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )
mf_result_list_test <- summary_stat_1 %>% 
  select(SNP, CHR, POS) %>% 
  full_join(summary_stat_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
  mutate(
    CHR = coalesce(CHR.x, CHR.y),
    POS = coalesce(POS.x, POS.y)
  ) %>%
  select(SNP, CHR, POS) %>% 
  left_join(MESuSiE_res_dataframe, by = "SNP") %>% 
  select(SNP,CHR, POS,MESuSiE_PIP_Either, MESuSiE_PIP_Shared, MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_cs) %>% arrange(CHR,POS) %>% replace(is.na(.), 0) 


zfile %>% filter(RSID %in% mf_result_list$SNP[mf_result_list$SuSiE_cs!=0]) %>% pull(Signal) %>% table()
zfile %>% filter(RSID %in% mf_result_list_test$SNP[mf_result_list_test$MESuSiE_cs!=0]) %>% pull(Signal) %>% table()



for(causal_index in 2){
  for(External_index in 1:2){
    for (num_causal in 1){
      for(h2_num in 1:2){
        final_results_list <- list()
        
        for(LD_BLOCK in 1:100){
          causal_index_name = c("Both","One")[causal_index]
          External_index_name = c("","External_")[External_index]

  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  system(paste0("mkdir -p ",data_dir))
  result_dir<-paste0(wrk_dir,"result/")
  system(paste0("mkdir -p ",result_dir))
  paintor_dir<-data_dir
  
  
  zfile<-read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
  EU_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
  colnames(EU_cov)<-zfile$RSID
  rownames(EU_cov)<-zfile$RSID
  BB_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
  colnames(BB_cov)<-zfile$RSID
  rownames(BB_cov)<-zfile$RSID
  
  
  if(External_index ==1 ){
    num_snp_missing = round(nrow(zfile)*0.3)
  }else if(External_index ==2){
    num_snp_missing = round(nrow(zfile)*0.3+nrow(zfile)*(1-0.3)*0.3)
  }
  
  
  non_causal_SNP = zfile%>%filter(Signal==0)%>%pull(RSID)
  EU_missing_SNP = c(sample(non_causal_SNP,num_snp_missing-1),zfile%>%filter(Signal!=0)%>%pull(RSID))
  BB_missing_SNP = EU_missing_SNP
  zfile<-zfile%>%mutate(EU_missing = ifelse(RSID%in%EU_missing_SNP,1,0),BB_missing = ifelse(RSID%in%BB_missing_SNP,1,0))
  
  
  
  # DETERMINE BLOCK ID:
  current_block_id <- LD_BLOCK
  
  # Use tryCatch to skip regions where files might be missing without stopping the whole script
  skip_to_next <- FALSE
  
  cat(paste0("Processing Block: ", current_block_id, "... "))
  
  
  tryCatch({
    
    # --- Run Decision Logic ---
    
    if(causal_index==1){
      susie_EU_zfile<-zfile%>%filter(EU_missing==0)
      susie_BB_zfile<-zfile%>%filter(BB_missing==0)
      missing_EU_index = which(zfile$EU_missing!=0)
      missing_BB_index = which(zfile$BB_missing!=0)
      susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
      susie_BB_cov<-BB_cov[-missing_BB_index,-missing_BB_index]
    }else if(causal_index==2){
      susie_EU_zfile<-zfile%>%filter(EU_missing==0)
      susie_BB_zfile<-zfile
      missing_EU_index = which(zfile$EU_missing!=0)
      susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
      susie_BB_cov<-BB_cov
    }
    
    
    summary_stat_1 = data.frame("SNP" = susie_EU_zfile$RSID,"CHR" = susie_EU_zfile$CHR,"POS" = susie_EU_zfile$POS,"Signal" = susie_EU_zfile$Signal,"Beta"=susie_EU_zfile$zscore_1/sqrt(num_EU),"Se"=1/sqrt(num_EU), "Z" =susie_EU_zfile$zscore_1 ,  "N" =num_EU ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
    summary_stat_2 = data.frame("SNP" = susie_BB_zfile$RSID,"CHR" = susie_BB_zfile$CHR,"POS" = susie_BB_zfile$POS,"Signal" = susie_BB_zfile$Signal,"Beta"=susie_BB_zfile$zscore_2/sqrt(num_BB),"Se"=1/sqrt(num_BB), "Z" =susie_BB_zfile$zscore_2 ,  "N" =num_BB ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
    
    setDT(summary_stat_1)
    setDT(summary_stat_2)
    
    decision <- decide_finemapping_method(sum_eur = summary_stat_1, sum_afr = summary_stat_2, ld_eur = susie_EU_cov, ld_afr = susie_BB_cov,r2_thresh = 0.6)
    
    
    # Store Result
    final_results_list[[LD_BLOCK]] <- data.table(
      Locus_ID = current_block_id,
      Method = decision$method,
      Reason = decision$reason,
      Reason_code = decision$reason_code,
      Status = "Success"
    )
    
    cat(paste0("Decision: ", decision$method, "\n"))
    
    
    
  }, error = function(e) {
    # Error Handler
    cat(paste0("ERROR: ", e$message, "\n"))
    final_results_list[[LD_BLOCK]] <<- data.table(
      Locus_ID = current_block_id,
      Method = NA,
      Reason = paste("Error:", e$message),
      Status = "Failed"
    )
  })
        }
        final_results_df <- rbindlist(final_results_list)
        output_file <- paste0(result_dir, "Finemapping_Method_Decisions_h2_",h2_num, ".csv")
        fwrite(final_results_df, output_file)
        print(paste0("All done. Results saved to: ", output_file))
        
}}}}



# --- Compile and Save Results ---
final_results_df <- rbindlist(final_results_list)

# Generate timestamp for filename
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- paste0(result_dir, "Finemapping_Method_Decisions_", trait, "_", timestamp, ".csv")

fwrite(final_results_df, output_file)
print(paste0("All done. Results saved to: ", output_file))


