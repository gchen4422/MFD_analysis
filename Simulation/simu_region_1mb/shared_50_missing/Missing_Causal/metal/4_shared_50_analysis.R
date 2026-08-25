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
library(tidyr)
library(MESuSiE)





wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
paintor_dir<-data_dir


##############################################
#
#		Format METAL Input
#
##############################################
zfile<-read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
eur_freq = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_",LD_BLOCK,".bim"),header=F) %>% select(V2,V5,V6) %>% dplyr::rename(RSID = V2,A1=V5,A2=V6)
zfile = zfile %>% left_join(eur_freq, by = "RSID")
#zfile_susiex_input = fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"))
#zfile$EU_missing = ifelse(zfile$RSID %in%zfile_susiex_input$RSID,0,1)
#zfile$BB_missing = zfile$EU_missing
SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")
load(SuSiE_res_name)
zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>% select(RSID,EU_missing,BB_missing)
zfile = zfile %>% left_join(zfile_MESuSiE_SuSiE_Paintor, by = "RSID")

#EU_cov<-as.matrix(fread(paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
#colnames(EU_cov)<-zfile$RSID
#BB_cov<-as.matrix(fread(paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
#colnames(BB_cov)<-zfile$RSID
merged_bim_path <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_merge/loci_", LD_BLOCK, ".bim")
merged_bim <- fread(merged_bim_path, header=FALSE)
if (!all(merged_bim$V5 == zfile$A1 & 
         merged_bim$V6 == zfile$A2)) {
  stop("CHECK 2 FAILED error: Allele mismatch! Merged LD alleles do not match zfile.")
}
cat("Check 2 Passed: LD Matrix is perfectly aligned to zfile.\n")

merged_cov = as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_merge/loci_",LD_BLOCK,".ld")))  ## In-sample LD matrix after merge bfile of eur and afr
colnames(merged_cov)<-zfile$RSID


#n_eur = median(zfile$N_1)
#n_afr = median(zfile$N_2)
#weighted_cov <- (n_eur * EU_cov + n_afr * BB_cov) / (n_eur + n_afr)
#colnames(weighted_cov)<-zfile$RSID


if(causal_index==1){
  missing_index = which(zfile$EU_missing!=0|zfile$BB_missing!=0)
  zfile_subset<-zfile%>%filter(EU_missing==0,BB_missing==0)
  merged_cov_subset<-merged_cov[-missing_index,-missing_index]
}else if(causal_index==2){
  #missing_index = which(zfile$EU_missing!=0)
  zfile_subset<-zfile
  merged_cov_subset<-merged_cov
}




metal_results = fread(paste0(result_dir,"Metal_",num_causal,"_h2_",h2_num,"_result1.tbl"))
#metal_results_loci <- metal_results %>%
#  filter(MarkerName %in% zfile_subset$RSID) %>%
#  separate(MarkerName, into = c("CHR","POS"), sep = ":", remove = FALSE, convert = TRUE) %>% arrange(CHR, POS) %>% 
#  mutate(Z = qnorm(log(`P-value`/2), lower.tail = FALSE, log.p = TRUE) * sign(Effect))

metal_results_loci_raw <- metal_results %>%
  filter(MarkerName %in% zfile_subset$RSID) %>%
  mutate(
    Allele1 = toupper(Allele1),
    Allele2 = toupper(Allele2)
  ) %>%
  separate(MarkerName, into = c("CHR","POS"), sep = ":", remove = FALSE, convert = TRUE) %>% arrange(CHR, POS) %>% 
  mutate(Z = qnorm(log(`P-value`/2), lower.tail = FALSE, log.p = TRUE) * sign(Effect)) 

metal_results_loci_aligned <- metal_results_loci_raw %>%
  mutate(
    # Define Target (LD) and Current (METAL) references for comparison
    Target_Ref = zfile_subset$A1,
    Target_Alt = zfile_subset$A2,
    Current_Ref = Allele1,
    
    # Logic:
    # 1. Direct Match: METAL Ref == LD Ref -> Keep Z (flip_index = 1)
    # 2. Swapped:      METAL Ref == LD Alt -> Flip Z (flip_index = -1)
    # 3. Mismatch:     Neither (e.g. strand error) -> NA
    flip_index = case_when(
      Current_Ref == Target_Ref ~ 1,
      Current_Ref == Target_Alt ~ -1,
      TRUE ~ NA_real_
    ),
    
    # Apply the flip
    Z_final = Z * flip_index,
    
    # Optional: Update the allele columns to match the LD reference
    Allele1 = Target_Ref,
    Allele2 = Target_Alt
  ) %>%
  # Safety: Remove any SNPs that couldn't be matched (flip_index is NA)
  filter(!is.na(flip_index))

if (!all(metal_results_loci_aligned$Allele1 == zfile_subset$A1 & 
         metal_results_loci_aligned$Allele2 == zfile_subset$A2)) {
  stop("CHECK 1 FAILED error: Allele mismatch! METAL results alleles do not match zfile.")
}
cat("Check 1 Passed: METAL results are perfectly aligned to zfile.\n")

if (!all(metal_results_loci_aligned$MarkerName == zfile_subset$RSID)) {
  stop("Check 0 FAILED error: RSID Order Mismatch despite fix!")
}


metal_results_loci = metal_results_loci_aligned %>% select(MarkerName,Z_final) %>% dplyr::rename(Z= Z_final)

##############################################
#
#		Run SuSiE
#
##############################################
library(susieR)
start.time<-Sys.time()
#susie_weighted<-susie_rss(metal_results_loci$Z,weighted_cov,check_prior=F, n = 600000)
susie_merged<-susie_rss(metal_results_loci$Z,merged_cov_subset,check_prior=F, n = 600000)


susie_merged<-susie_rss(metal_results_loci$Z,merged_cov_subset,check_prior=F, n = 600000,L =1)

# Start with full data
#Z_cur   <- metal_results_loci$Z
#cov_cur <- merged_cov_subset

#all_outliers <- integer(0)

#repeat {
#  EU_diagnostic <- kriging_rss(Z_cur, cov_cur)
#  outlier_index <- EU_diagnostic$plot$plot_env$idx
#  
#  # stop if no outliers
#  if (is.null(outlier_index) || length(outlier_index) == 0) break
  
#   # record outliers in original indexing space if needed
#   # (optional: only when tracking mappings)
#   all_outliers <- c(all_outliers, outlier_index)
#   
#   # keep everything except current outliers
#   keep_idx <- setdiff(seq_along(Z_cur), outlier_index)
#   
#   Z_cur   <- Z_cur[keep_idx]
#   cov_cur <- cov_cur[keep_idx, keep_idx, drop = FALSE]
# }
# 
# 
# susie_merged<-susie_rss(Z_cur,cov_cur,check_prior=F, n = 600000)

SuSiE_name<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_updated.RData")


#save(susie_weighted,susie_merged, file = SuSiE_name)

save(susie_merged, file = SuSiE_name)

#time_out_name<-paste0(result_dir,"time_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
#time_vec<-c(MESuSiE.time,susie.time,Paintor.time, xmap.time, susiex.time,carmax.time)
#write.table(time_vec,time_out_name,col.names = F,row.names = F,quote=F)







