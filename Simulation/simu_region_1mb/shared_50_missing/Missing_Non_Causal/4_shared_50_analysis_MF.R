library(dplyr)
library(data.table)
library(MESuSiE)
library(susieR)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/MFD_utility.R")


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

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
paintor_dir<-data_dir

ld_EU_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
ld_BB_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")

zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
EU_cov<-as.matrix(fread(ld_EU_name))
colnames(EU_cov)<-zfile$RSID
rownames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(ld_BB_name))
colnames(BB_cov)<-zfile$RSID
rownames(BB_cov)<-zfile$RSID




  
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
  
  
  #susie_BB_cs<-rep(0,length(susie_BB$pip))
  #susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
  #for(j in seq_along(susie_BB$sets$cs)) {
  #  susie_BB_cs[susie_BB$sets$cs[[j]]] <- j
  #}
  susie_BB_cs <- get_cs_index_vector(
    cs_res = susie_BB$sets, 
    n_snps = length(susie_BB$pip), # Or nrow(your_dataframe)
    renumber = TRUE
  )
  #susie_EU_cs<-rep(0,length(susie_EU$pip))
  #susie_EU_cs[unlist(susie_EU$sets$cs)]<-1
  #for(k in seq_along(susie_EU$sets$cs)) {
  #  susie_EU_cs[susie_EU$sets$cs[[k]]] <- k
  #}
  
  susie_EU_cs <- get_cs_index_vector(
    cs_res = susie_EU$sets, 
    n_snps = length(susie_EU$pip), # Or nrow(your_dataframe)
    renumber = TRUE
  )
  
      
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
    select(SNP,CHR, POS,SuSiE_PIP_Either,SuSiE_PIP_Shared,SuSiE_PIP_EU,SuSiE_PIP_BB,SuSiE_cs,susie_EU_cs,susie_BB_cs)%>% arrange(CHR,POS) 
    
      
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
      
      
  #MESuSiE_cs<-rep(0,length(MESuSiE_res$pip))
  #for(m in seq_along(MESuSiE_res$cs$cs)) {
  #  MESuSiE_cs[MESuSiE_res$cs$cs[[m]]] <- m
  #}
  
  MESuSiE_cs = get_cs_index_vector(
    cs_res = MESuSiE_res$cs, 
    n_snps = length(MESuSiE_res$pip), # Or nrow(your_dataframe)
    renumber = TRUE
  )
  
  #wb_cs_res <- meSuSie_get_cs_specific(
  #  res = MESuSiE_res, 
  #  Xcorr = R_mat_list[1], # Or whichever LD matrix is appropriate
  #  target_idx = 1,          # 3 usually = Shared
  #  cor_threshold = 0.5
  #)
  
  MESuSiE_cs_WB_vec <- rep(0, length(MESuSiE_res$pip)) 
  
  if (any(MESuSiE_res$cs$cs_category == "EUR")) {
    
    for (j in seq_along(MESuSiE_res$cs$cs[MESuSiE_res$cs$cs_category == "EUR"])) {
      
      new_effect_id <- j 
      snp_indices <- MESuSiE_res$cs$cs[MESuSiE_res$cs$cs_category == "EUR"][[j]]
      MESuSiE_cs_WB_vec[snp_indices] <- new_effect_id
    }
  }
  
  
  #bb_cs_res <- meSuSie_get_cs_specific(
  #  res = MESuSiE_res, 
  #  Xcorr = R_mat_list[2], # Or whichever LD matrix is appropriate
  #  target_idx = 1,          # 3 usually = Shared
  #  cor_threshold = 0.5
  #)
  
  MESuSiE_cs_BB_vec <- rep(0, length(MESuSiE_res$pip)) 
  
  if (any(MESuSiE_res$cs$cs_category == "AFR")) {
    
    for (j in seq_along(MESuSiE_res$cs$cs[MESuSiE_res$cs$cs_category == "AFR"])) {
      
      new_effect_id <- j 
      snp_indices <- MESuSiE_res$cs$cs[MESuSiE_res$cs$cs_category == "AFR"][[j]]
      MESuSiE_cs_BB_vec[snp_indices] <- new_effect_id
    }
  }
  
  
  MESuSiE_res_dataframe<-data.frame(SNP = summary_stat_1_common$SNP,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs,MESuSiE_WB_cs =MESuSiE_cs_WB_vec, MESuSiE_BB_cs =MESuSiE_cs_BB_vec )
  mf_result_list <- summary_stat_1 %>% 
    select(SNP, CHR, POS) %>% 
    full_join(summary_stat_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
    mutate(
          CHR = coalesce(CHR.x, CHR.y),
          POS = coalesce(POS.x, POS.y)
    ) %>%
    select(SNP, CHR, POS) %>% 
    left_join(MESuSiE_res_dataframe, by = "SNP") %>% 
    select(SNP,CHR, POS,MESuSiE_PIP_Either, MESuSiE_PIP_Shared, MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_cs,MESuSiE_WB_cs,MESuSiE_BB_cs) %>% arrange(CHR,POS) %>% replace(is.na(.), 0) 
      
}

colnames(mf_result_list) = c("SNP","CHR","POS","PIP_Either","PIP_Shared","PIP_Ancestry_1","PIP_Ancestry_2","CS","CS_Ancestry_1","CS_Ancestry_2")   




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
MESuSiE_res_dataframe<-data.frame(SNP = summary_stat_eur_common$SNP,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )
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
    for (num_causal in 1:3){
      for(h2_num in 1:2){
        final_results_list <- list()
        
          for(LD_BLOCK in 1:100){
          causal_index_name = c("Both","One")[causal_index]
          External_index_name = c("","External_")[External_index]
  
  wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
  system(paste0("mkdir -p ",wrk_dir))
  data_dir<-paste0(wrk_dir,"summary_data/")
  system(paste0("mkdir -p ",data_dir))
  result_dir<-paste0(wrk_dir,"result/")
  system(paste0("mkdir -p ",result_dir))
  paintor_dir<-data_dir
  
  ld_EU_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
  ld_BB_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
  
  zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
  EU_cov<-as.matrix(fread(ld_EU_name))
  colnames(EU_cov)<-zfile$RSID
  rownames(EU_cov)<-zfile$RSID
  BB_cov<-as.matrix(fread(ld_BB_name))
  colnames(BB_cov)<-zfile$RSID
  rownames(BB_cov)<-zfile$RSID
  
  
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
    
    decision <- decide_finemapping_method(summary_stat_1, summary_stat_2, susie_EU_cov, susie_BB_cov)
    

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
