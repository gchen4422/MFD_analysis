simulation_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/"
plot_dir<-paste0(simulation_dir,"Figure/")
system(paste0("mkdir -p ",plot_dir))
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/")
system(paste0("mkdir -p ",res_dir))
res_out<-"shared_50_external_ld_updated_xmap_meta_slalom.RData"
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)
library(XMAP)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/utility.R")
num_BB = 300000
num_EU = 300000
find_index_until_threshold <- function(vec, threshold) {
  # Sort the input vector in decreasing order and obtain the sorted indices
  sorted_indices <- order(vec, decreasing = TRUE)
  
  # Calculate the cumulative sum of the sorted vector
  cumsum_vec <- cumsum(vec[sorted_indices])
  
  # Find the index where the cumulative sum reaches or exceeds the threshold
  index <- which(cumsum_vec >= threshold)[1]
  
  # Return the original indices of the elements up to the found index
  sorted_indices[1:index]
}

################################################################
#
#        Result Summarizing
#
#
################################################################
ancestry_all_ROC_data_list<-list()
either_all_ROC_data_list<-list()
shared_all_ROC_data_list<-list()
all_Set_data_list<-list()
pip_res_all_list<-list()
time_all<-c()
data_all<-c()
#susiex_runtime = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/susiex_result/susiex_runtime_log.txt", header = F)
for(causal_num in 1:3){
  for(h2 in 1:2){
    wrk_dir<-paste0(simulation_dir,"causal_num_",causal_num,"/")
    data_dir<-paste0(wrk_dir,"summary_data/slalom/")
    result_dir<-paste0(wrk_dir,"result/slalom/")
    result_dir_meta<-paste0(wrk_dir,"result/")
    
    setwd(result_dir)
    set1 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
    set2 = which(paste0("CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_slalom.mcmc.paintor")%in%list.files())
    replicate_set = Reduce(intersect, list(set1,set2))
    #replicate_set = rep(1:41)
    empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
    
    MESuSiE_ROC <- SuSiE_ROC <- SuSiE_weighted_ROC <- SuSiE_merged_ROC <- Paintor_ROC <- susiex_ROC <- XMAP_ROC <- multisusie_ROC <- carmax_ROC <- 
      MESuSiE_ROC_Either <- SuSiE_ROC_Either <- SuSiE_weighted_ROC_Either <- SuSiE_merged_ROC_Either<- Paintor_ROC_Either <- susiex_ROC_Either <- 
      XMAP_ROC_Either <- multisusie_ROC_Either <- carmax_ROC_Either <- 
      MESuSiE_ROC_EU <- MESuSiE_ROC_BB <- SuSiE_ROC_EU <- SuSiE_ROC_BB <- Paintor_ROC_EU <- Paintor_ROC_BB <- 
      SuSiE_weighted_ROC_EU <- SuSiE_weighted_ROC_BB <- SuSiE_merged_ROC_EU <- SuSiE_merged_ROC_BB <-
      susiex_ROC_EU <- susiex_ROC_BB <- XMAP_ROC_EU <- XMAP_ROC_BB <- 
      multisusie_ROC_EU <- multisusie_ROC_BB <- carmax_ROC_EU <- carmax_ROC_BB <- empty_mat
    
    set_res_all<-c()
    pip_res_all<-c()  
    
    
    
    for(LOCI_num in replicate_set){
      
      SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      load(SuSiE_res_name)
      zfile_ori = fread(paste0(wrk_dir,"summary_data/","CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2))%>%dplyr::rename(SNP = RSID)
      paintor_res<-fread(paste0(result_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_slalom.mcmc.paintor"))%>%dplyr::rename(SNP = RSID)
      SuSiE_res_name_metal<-paste0(result_dir_meta,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_correction_updated.RData")
      load(SuSiE_res_name_metal)
      #XMAP_res_name<-paste0(result_dir,"XMAP_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      #load(XMAP_res_name)
      
      
      #format mesusie results
      
      MESuSiE_cs<-rep(0,length(MESuSiE_res$pip))
      MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
      MESuSiE_res_dataframe<-data.frame(SNP = paintor_res$SNP,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )
      
      zfile_complete <- zfile_ori %>%
        left_join(MESuSiE_res_dataframe, by = "SNP")
      
      
      #format susie results
      
      susie_BB_cs<-rep(0,length(susie_BB$pip))
      susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
      susie_EU_cs<-rep(0,length(susie_EU$pip))
      susie_EU_cs[unlist(susie_EU$sets$cs)]<-1
      
      SuSiE_EU_res_dataframe<-data.frame(SNP = paintor_res$SNP,SuSiE_PIP_EU = susie_EU$pip,susie_EU_cs = susie_EU_cs )
      SuSiE_BB_res_dataframe<-data.frame(SNP = paintor_res$SNP,SuSiE_PIP_BB = susie_BB$pip,susie_BB_cs = susie_BB_cs )
      
      
      zfile_complete <- zfile_complete %>%
        left_join(SuSiE_EU_res_dataframe, by = "SNP")%>%
        left_join(SuSiE_BB_res_dataframe, by = "SNP")
      
      zfile_complete <- zfile_complete%>% replace(is.na(.), 0)
      zfile_complete <- zfile_complete%>%mutate(SuSiE_cs = ifelse(susie_EU_cs+susie_BB_cs==0,0,1),SuSiE_PIP_Either = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),SuSiE_PIP_Shared = pmin(SuSiE_PIP_EU,SuSiE_PIP_BB))
      
      #format susie weighted results
      susie_weighted_cs<-rep(0,length(susie_weighted_slalom$pip))
      susie_weighted_cs[unlist(susie_weighted_slalom$sets$cs)]<-1
      
      SuSiE_weighted_res_dataframe<-data.frame(SNP = metal_results_loci_sub_slalom$MarkerName,SuSiE_Weighted_PIP = susie_weighted_slalom$pip,SuSiE_weighted_cs = susie_weighted_cs )
      
      zfile_complete <- zfile_complete %>%
        left_join(SuSiE_weighted_res_dataframe, by = "SNP")
      
      #format susie merged results
      susie_merged_cs<-rep(0,length(susie_merged_slalom$pip))
      susie_merged_cs[unlist(susie_merged_slalom$sets$cs)]<-1
      
      SuSiE_merged_res_dataframe<-data.frame(SNP = metal_results_loci_sub_slalom$MarkerName,SuSiE_Merged_PIP = susie_merged_slalom$pip,SuSiE_merged_cs = susie_merged_cs )
      
      zfile_complete <- zfile_complete %>%
        left_join(SuSiE_merged_res_dataframe, by = "SNP")
      
      #format susiex results
      susiex_res_name  <- paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,
                                 "_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95_slalom.snp")
      susiex_cred_name <- paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,
                                 "_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95_slalom.cs")
      
      
      
      if (file.exists(susiex_res_name) && file.exists(susiex_cred_name)) {
        ## files exist: normal path
        susiex_res <- read.table(susiex_res_name, header = TRUE)
        susiex_cs  <- read.table(susiex_cred_name, header = TRUE)
        
        susiex_cs_num <- if (nrow(susiex_cs) > 0) max(susiex_cs$CS_ID) else 0L
        susiex_cred   <- which(paintor_res$POS %in% susiex_cs$BP)
        
        # df: data.frame with per-CS PIP columns for all SNPs
        susiex_pip_cols <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$",
                                names(susiex_res), value = TRUE)
        
        if (length(susiex_pip_cols) > 0) {
          # fast matrix calc with stability (row-wise)
          M <- as.matrix(susiex_res[, susiex_pip_cols, drop = FALSE])
          susiex_res$pip <- 1 - apply(1 - M, 1, prod)  # overall pip per SNP
        } else {
          # no PIP columns found → set all pips to 0
          susiex_res$pip <- rep(0, nrow(susiex_res))
        }
        
      } else {
        ## one or both files missing: define empty cred and pip = 0 for each SNP
        susiex_cs_num <- 0L
        susiex_cred   <- integer(0)
        
        # to obtain one PIP per SNP in paintor_res:
        susiex_res <- paintor_res    # or data.frame(pip = rep(0, nrow(paintor_res)))
        susiex_res$pip <- rep(0, nrow(paintor_res))
      }
      susiex_res_to_merge = susiex_res %>% select(SNP,pip) %>% mutate(SuSiEx_cs = 0) %>% dplyr::rename(SuSiEx_PIP_Either = "pip")
      
      susiex_res_to_merge$SuSiEx_cs[susiex_cred] <- 1
      
      zfile_complete <- zfile_complete %>%
        left_join(susiex_res_to_merge, by = "SNP")
      
      
      
      
      
      #format xmap results
      EU_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_slalom.LD1")))
      BB_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_slalom.LD2")))
      
      
      cs1 <- get_CS(xmap, Xcorr = EU_cov, coverage = 0.95, min_abs_corr = 0.5)
      cs2 <- get_CS(xmap, Xcorr = BB_cov, coverage = 0.95, min_abs_corr = 0.5)
      cs_xmap <- unique(unlist(cs1$cs[intersect(names(cs1$cs), names(cs2$cs))],use.names = FALSE))
      pip_xmap <- get_pip(xmap$gamma)
      
      
      xmap_res_dataframe<-data.frame(SNP = paintor_res$SNP,XMAP_PIP_Either = pip_xmap, XMAP_cs = 0)
      xmap_res_dataframe$XMAP_cs[cs_xmap] <- 1
      
      zfile_complete <- zfile_complete %>%
        left_join(xmap_res_dataframe, by = "SNP") %>% replace(is.na(.), 0)
      
      #format multisusie results
      multisusie_dir = paste0(result_dir,"multisusiex_result/")
      multisusie_snp_path <- paste0(
        multisusie_dir,
        "MultiSuSiE_CAUSAL_", causal_num,
        "_LOCI_",          LOCI_num,
        "_h2_",            h2,
        "_output_snp.txt"
      )
      
      multisusie_cs_path <- paste0(
        multisusie_dir,
        "MultiSuSiE_CAUSAL_", causal_num,
        "_LOCI_",          LOCI_num,
        "_h2_",            h2,
        "_output_cs.txt"
      )
      
      if (file.exists(multisusie_snp_path) && file.exists(multisusie_cs_path)) {
        ## Normal case: both files present
        multisusie_snp <- fread(multisusie_snp_path)
        multisusie_cs  <- fread(multisusie_cs_path)
        multisusie_cred <- which(paintor_res$POS %in% multisusie_cs$POS)
        pip_multisusie = multisusie_snp$pip
      } else {
        multisusie_cred <- integer(0)
        pip_multisusie = rep(0, nrow(paintor_res))
        
      }
      multisusie_snp_to_merge = multisusie_snp %>% select(RSID,pip) %>% dplyr::rename(MultiSuSiE_PIP_Either = "pip",SNP = RSID) %>% mutate(MultiSuSiE_cs = 0)
      multisusie_snp_to_merge$MultiSuSiE_cs[multisusie_cred] <- 1
      
      zfile_complete <- zfile_complete %>%
        left_join(multisusie_snp_to_merge, by = "SNP") %>% replace(is.na(.), 0)
      
      
      #format carmax results
      
      
      carmax_cs <- tryCatch({
        cs_list_1 <- carmax_results[[1]][["Credible set"]][[2]]
        cs_list_2 <- carmax_results[[2]][["Credible set"]][[2]]
        cs_1 <- if (is.list(cs_list_1) && length(cs_list_1) > 0) sort(unique(unlist(cs_list_1))) else integer(0)
        cs_2 <- if (is.list(cs_list_2) && length(cs_list_2) > 0) sort(unique(unlist(cs_list_2))) else integer(0)
        sort(intersect(cs_1, cs_2))
      }, error = function(e) integer(0))
      carmax_pip_eur <- carmax_results[[1]]$PIPs
      carmax_pip_afr <- carmax_results[[2]]$PIPs
      carmax_pip_either <- pmax(carmax_pip_eur, carmax_pip_afr, na.rm = TRUE)
      carmax_pip_shared <- pmin(carmax_pip_eur, carmax_pip_afr, na.rm = TRUE)
      
      carmax_res_dataframe<-data.frame(SNP = paintor_res$SNP,CARMAX_PIP_Either = carmax_pip_either, CARMAX_PIP_WB =carmax_pip_eur ,CARMAX_PIP_BB = carmax_pip_afr ,CARMAX_PIP_Shared =carmax_pip_shared, CARMAX_cs = 0)
      carmax_res_dataframe$CARMAX_cs[carmax_cs] <- 1
      zfile_complete <- zfile_complete %>%
        left_join(carmax_res_dataframe, by = "SNP") %>% replace(is.na(.), 0)
      
      
      #format PAINTOR results
      
      if(sum(paintor_res$Posterior_Prob)==0){
        paintor_cs<-rep(0,nrow(paintor_res))
      }else{
        paintor_cs_index<-find_index_until_threshold(paintor_res$Posterior_Prob/sum(paintor_res$Posterior_Prob),0.95)
        paintor_cs<-rep(0,nrow(paintor_res))
        paintor_cs[paintor_cs_index]<-1
      }
      
      paintor_res_dataframe<-data.frame(SNP = paintor_res$SNP,Paintor_PIP_Either = paintor_res$Posterior_Prob,Paintor_cs = paintor_cs)
      
      zfile_complete <- zfile_complete %>%
        left_join(paintor_res_dataframe, by = "SNP")%>% replace(is.na(.), 0)
      
      
      #Set Size/Power
      MESuSiE_SET = compute_set_susie(which(zfile_complete$MESuSiE_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))		
      SuSiE_SET = compute_set_susie(which(zfile_complete$SuSiE_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))	
      SuSiE_weighted_SET = compute_set_susie(which(zfile_complete$SuSiE_weighted_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))			
      SuSiE_merged_SET = compute_set_susie(which(zfile_complete$SuSiE_merged_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))			
      Paintor_SET = compute_set_pip(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob,"CAUSAL"=ifelse(paintor_res$Signal!=0,1,0)))	
      susiex_SET = compute_set_susie(which(zfile_complete$SuSiEx_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))	
      XMAP_SET = compute_set_susie(which(zfile_complete$XMAP_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))		
      multisusie_SET = compute_set_susie(which(zfile_complete$MultiSuSiE_cs!= 0),ifelse(zfile_complete$Signal!=0,1,0))	
      carmax_SET = compute_set_susie(which(zfile_complete$CARMAX_cs != 0),ifelse(zfile_complete$Signal!=0,1,0))	
      
      
      
      set_res_all<-rbind(set_res_all,rbind(c(MESuSiE_SET,"MESuSiE"),c( SuSiE_SET,"SuSiE"),c( SuSiE_weighted_SET,"SuSiE_meta_weighted"),c( SuSiE_merged_SET,"SuSiE_meta_merged"),c(Paintor_SET$Paintor_set,"Paintor"),c(susiex_SET,"SuSiEx"),c(XMAP_SET,"XMAP"),c(multisusie_SET,"MultiSuSiE"),c(carmax_SET,"CARMAX")))
      
      
      
      #ROC for shared/either/ancestry-specific
      MESuSiE_ROC = MESuSiE_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MESuSiE_PIP_Shared),ifelse(zfile_complete$Signal==3,1,0))
      SuSiE_ROC = SuSiE_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_PIP_Shared),ifelse(zfile_complete$Signal==3,1,0))
      SuSiE_weighted_ROC = SuSiE_weighted_ROC + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Weighted_PIP),ifelse(zfile_complete$Signal==3,1,0))
      SuSiE_merged_ROC = SuSiE_merged_ROC + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Merged_PIP),ifelse(zfile_complete$Signal==3,1,0))
      Paintor_ROC = Paintor_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$Paintor_PIP_Either),ifelse(zfile_complete$Signal==3,1,0))
      susiex_ROC = susiex_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiEx_PIP_Either),ifelse(zfile_complete$Signal==3,1,0))
      XMAP_ROC = XMAP_ROC + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$XMAP_PIP_Either),ifelse(zfile_complete$Signal==3,1,0))
      multisusie_ROC = multisusie_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MultiSuSiE_PIP_Either),ifelse(zfile_complete$Signal==3,1,0))
      carmax_ROC = carmax_ROC + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$CARMAX_PIP_Shared),ifelse(zfile_complete$Signal==3,1,0))
      
      
      
      MESuSiE_ROC_Either = MESuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MESuSiE_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      SuSiE_ROC_Either = SuSiE_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      SuSiE_weighted_ROC_Either = SuSiE_weighted_ROC_Either + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Weighted_PIP),ifelse(zfile_complete$Signal!=0,1,0))
      SuSiE_merged_ROC_Either = SuSiE_merged_ROC_Either + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Merged_PIP),ifelse(zfile_complete$Signal!=0,1,0))
      Paintor_ROC_Either = Paintor_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$Paintor_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      susiex_ROC_Either = susiex_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiEx_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      XMAP_ROC_Either = XMAP_ROC_Either + compute_ROC(data.frame("SNP" = zfile_complete$SNP ,"PIP" = zfile_complete$XMAP_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      multisusie_ROC_Either = multisusie_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MultiSuSiE_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      carmax_ROC_Either = carmax_ROC_Either + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$CARMAX_PIP_Either),ifelse(zfile_complete$Signal!=0,1,0))
      
      
      
      MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MESuSiE_PIP_WB),ifelse(zfile_complete$Signal==1,1,0))
      MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MESuSiE_PIP_BB),ifelse(zfile_complete$Signal==2,1,0))
      
      SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_PIP_EU),ifelse(zfile_complete$Signal==1,1,0))
      SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_PIP_BB),ifelse(zfile_complete$Signal==2,1,0))
      
      SuSiE_weighted_ROC_EU =SuSiE_weighted_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Weighted_PIP),ifelse(zfile_complete$Signal==1,1,0))
      SuSiE_weighted_ROC_BB =SuSiE_weighted_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Weighted_PIP),ifelse(zfile_complete$Signal==2,1,0))  
      
      SuSiE_merged_ROC_EU =SuSiE_merged_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Merged_PIP),ifelse(zfile_complete$Signal==1,1,0))
      SuSiE_merged_ROC_BB =SuSiE_merged_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiE_Merged_PIP),ifelse(zfile_complete$Signal==2,1,0))  
      
      Paintor_ROC_EU = Paintor_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$Paintor_PIP_Either),ifelse(zfile_complete$Signal==1,1,0))
      Paintor_ROC_BB = Paintor_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$Paintor_PIP_Either),ifelse(zfile_complete$Signal==2,1,0))
      
      susiex_ROC_EU = susiex_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiEx_PIP_Either),ifelse(zfile_complete$Signal==1,1,0))
      susiex_ROC_BB = susiex_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$SuSiEx_PIP_Either),ifelse(zfile_complete$Signal==2,1,0))
      
      XMAP_ROC_EU = XMAP_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$XMAP_PIP_Either),ifelse(zfile_complete$Signal==1,1,0))
      XMAP_ROC_BB = XMAP_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$XMAP_PIP_Either),ifelse(zfile_complete$Signal==2,1,0))
      
      multisusie_ROC_EU = multisusie_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MultiSuSiE_PIP_Either),ifelse(zfile_complete$Signal==1,1,0))
      multisusie_ROC_BB = multisusie_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$MultiSuSiE_PIP_Either),ifelse(zfile_complete$Signal==2,1,0))
      
      carmax_ROC_EU = carmax_ROC_EU + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$CARMAX_PIP_WB),ifelse(zfile_complete$Signal==1,1,0))
      carmax_ROC_BB = carmax_ROC_BB + compute_ROC(data.frame("SNP" =zfile_complete$SNP ,"PIP" = zfile_complete$CARMAX_PIP_BB),ifelse(zfile_complete$Signal==2,1,0))
      
      
      
      
      ##Sumstat
      
      
      
      locus_sumstat <- zfile_complete %>%
        dplyr::select(1:6,
                      # Paintor
                      Paintor_cs = Paintor_cs,
                      Paintor_PIP = Paintor_PIP_Either,
                      
                      # CARMAX
                      CARMAX_cs = CARMAX_cs,
                      CARMAX_Either = CARMAX_PIP_Either,
                      CARMAX_Shared = CARMAX_PIP_Shared,
                      CARMAX_WB = CARMAX_PIP_WB,
                      CARMAX_BB = CARMAX_PIP_BB,
                      
                      # SuSiEx
                      SuSiEx_cs = SuSiEx_cs,
                      SuSiEx_PIP = SuSiEx_PIP_Either,
                      
                      # XMAP
                      XMAP_cs = XMAP_cs,
                      XMAP_PIP = XMAP_PIP_Either,
                      
                      # MultiSuSiE
                      MultiSuSIE_cs = MultiSuSiE_cs,
                      MultiSuSiE_PIP = MultiSuSiE_PIP_Either,
                      
                      # MESuSiE
                      MESuSiE_Either = MESuSiE_PIP_Either,
                      MESuSiE_Shared = MESuSiE_PIP_Shared,
                      MESuSiE_WB = MESuSiE_PIP_WB,
                      MESuSiE_BB = MESuSiE_PIP_BB,
                      MESuSiE_cs = MESuSiE_cs,
                      
                      # SuSiE Standard
                      SuSiE_Either = SuSiE_PIP_Either,
                      SuSiE_Shared = SuSiE_PIP_Shared,
                      SuSiE_WB = SuSiE_PIP_EU,
                      SuSiE_BB = SuSiE_PIP_BB,
                      SuSiE_cs = SuSiE_cs,
                      
                      # SuSiE Weighted
                      SuSiE_weighted_PIP = SuSiE_Weighted_PIP,
                      SuSiE_weighted_cs = SuSiE_weighted_cs,
                      
                      # SuSiE Merged
                      SuSiE_merged_PIP = SuSiE_Merged_PIP,
                      SuSiE_merged_cs = SuSiE_merged_cs,
                      
                      # Sample Sizes (Mapping N_1/N_2 to WB/BB)
                      N_WB = N_1,
                      N_BB = N_2
        ) %>%
        # Add external simulation parameters
        mutate(
          h2 = h2,
          causal_num = causal_num,
          locus = LOCI_num
        )
      data_all<-rbind(data_all,locus_sumstat)
      cat(LOCI_num)
    }
    
    #######################################
    #
    #   credible set combined 
    #
    ####################################### 
    set_res_all_data = data.frame("Power" = as.numeric(set_res_all[,1]),"FDR" = as.numeric(set_res_all[,2]),"Size"=as.numeric(set_res_all[,4]),"Method"=set_res_all[,5])
    set_res_all_data$h2 = h2
    set_res_all_data$causal_num = causal_num
    set_res_summary<-set_res_all_data %>% group_by(Method)%>%summarise_at(vars(Power,Size),list(name = mean)) 
    
    set_res_summary$h2 = h2
    set_res_summary$causal_num = causal_num
    all_Set_data_list[[paste0("causal_",causal_num,"_h2_",h2)]]<-set_res_all_data		
    
    #######################################
    #
    #   PIP combined 
    #
    ####################################### 			
    pip_res_all$h2 = h2
    pip_res_all$causal_num = causal_num
    pip_res_all_list[[paste0("causal_",causal_num,"_h2_",h2)]]<-pip_res_all
    #######################################
    #
    #   ROC combined 
    #
    ####################################### 
    # Either signal part
    roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either, SuSiE_ROC_Either = SuSiE_ROC_Either,SuSiE_weighted_ROC_Either = SuSiE_weighted_ROC_Either,SuSiE_merged_ROC_Either = SuSiE_merged_ROC_Either,Paintor_ROC_Either = Paintor_ROC_Either,susiex_ROC_Either = susiex_ROC_Either,XMAP_ROC_Either = XMAP_ROC_Either, multisusie_ROC_Either = multisusie_ROC_Either, carmax_ROC_Either= carmax_ROC_Either)
    method_names <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", "Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")
    either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
    
    # Shared signal part
    roc_data_list <- list(MESuSiE_ROC = MESuSiE_ROC, SuSiE_ROC = SuSiE_ROC,  SuSiE_weighted_ROC = SuSiE_weighted_ROC,  SuSiE_merged_ROC = SuSiE_merged_ROC, Paintor_ROC = Paintor_ROC,susiex_ROC = susiex_ROC,XMAP_ROC = XMAP_ROC, multisusie_ROC = multisusie_ROC,carmax_ROC= carmax_ROC)
    method_names <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", "Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")
    shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
    
    # Ancestry specific part
    roc_data_list <- list(MESuSiE_ROC_EU = MESuSiE_ROC_EU,MESuSiE_ROC_BB = MESuSiE_ROC_BB, SuSiE_ROC_EU = SuSiE_ROC_EU, SuSiE_ROC_BB = SuSiE_ROC_BB,SuSiE_weighted_ROC_EU = SuSiE_weighted_ROC_EU, SuSiE_weighted_ROC_BB = SuSiE_weighted_ROC_BB,SuSiE_merged_ROC_EU = SuSiE_merged_ROC_EU, SuSiE_merged_ROC_BB = SuSiE_merged_ROC_BB,Paintor_ROC_EU = Paintor_ROC_EU,Paintor_ROC_BB = Paintor_ROC_BB, susiex_ROC_EU = susiex_ROC_EU, susiex_ROC_BB = susiex_ROC_BB, XMAP_ROC_EU = XMAP_ROC_EU, XMAP_ROC_BB = XMAP_ROC_BB,multisusie_ROC_EU = multisusie_ROC_EU,multisusie_ROC_BB = multisusie_ROC_BB, carmax_ROC_EU =carmax_ROC_EU, carmax_ROC_BB = carmax_ROC_BB)
    method_names <- c("MESuSiE_EU", "MESuSiE_BB","SuSiE_EU", "SuSiE_BB","SuSiE_weighted_EU", "SuSiE_weighted_BB","SuSiE_merged_EU", "SuSiE_merged_BB","Paintor_EU", "Paintor_BB","SuSiEx_EU","SuSiEx_BB","XMAP_EU","XMAP_BB","MultiSuSiE_EU","MultiSuSiE_BB","CARMAX_EU","CARMAX_BB")
    ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
    
  }
}
#######################################
#
#  Set Size and Power
#
####################################### 

all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
all_Set_data_dataframe<-all_Set_data_dataframe%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_meta_weighted","SuSiE_meta_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
#levels(all_Set_data_dataframe$Method)<-c("MeSuSiE","Paintor","SuSiE")
cols <- c("Method","h2","causal_num")
set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
colnames(set_power_summary)[4]<-"Power_name"

set_power_summary<-set_power_summary%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "SuSiE_meta_weighted","SuSiE_meta_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
#levels(set_power_summary$Method)<-c("MESuSiE","SuSiE","Paintor")


#######################################
#
#  PIP calibration
#
####################################### 

####PIP Either

PIP_calibration_either<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
PIP_calibration_either<-PIP_calibration_either%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_meta_weighted" = "SuSiE_weighted_PIP","SuSiE_meta_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor", "SuSiEx","XMAP","MultiSuSiE","CARMAX"))


####PIP Shared	
PIP_calibration_Shared<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared,SuSiE_weighted_PIP , SuSiE_merged_PIP ,Paintor_PIP, SuSiEx_PIP, XMAP_PIP, MultiSuSiE_PIP,CARMAX_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Shared"))
PIP_calibration_Shared<-PIP_calibration_Shared%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiE_meta_weighted" = "SuSiE_weighted_PIP","SuSiE_meta_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP","SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"))

###PIP ancestry
data_all<-data_all%>%mutate(SuSiE_weighted_WB = SuSiE_weighted_PIP,SuSiE_weighted_BB = SuSiE_weighted_PIP)
data_all<-data_all%>%mutate(SuSiE_merged_WB = SuSiE_merged_PIP,SuSiE_merged_BB = SuSiE_merged_PIP)
data_all<-data_all%>%mutate(Paintor_WB = Paintor_PIP,Paintor_BB = Paintor_PIP)
data_all<-data_all%>%mutate(SuSiEx_WB = SuSiEx_PIP,SuSiEx_BB = SuSiEx_PIP)
data_all<-data_all%>%mutate(XMAP_WB = XMAP_PIP,XMAP_BB = XMAP_PIP)
data_all<-data_all%>%mutate(MultiSuSiE_WB = MultiSuSiE_PIP,MultiSuSiE_BB = MultiSuSiE_PIP)


PIP_calibration_WB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_WB,SuSiE_WB,SuSiE_weighted_WB,SuSiE_merged_WB,Paintor_WB,SuSiEx_WB, XMAP_WB, MultiSuSiE_WB, CARMAX_WB),c(1),c("MESuSiE_WB","SuSiE_WB","SuSiE_weighted_WB","SuSiE_merged_WB","Paintor_WB","SuSiEx_WB","XMAP_WB", "MultiSuSiE_WB","CARMAX_WB"))
PIP_calibration_BB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_BB,SuSiE_BB,SuSiE_weighted_BB,SuSiE_merged_BB,Paintor_BB,SuSiEx_BB, XMAP_BB, MultiSuSiE_BB, CARMAX_BB),c(2),c("MESuSiE_BB","SuSiE_BB","SuSiE_weighted_BB","SuSiE_merged_BB","Paintor_BB","SuSiEx_BB","XMAP_BB", "MultiSuSiE_BB","CARMAX_BB"))

PIP_calibration_ancestry<-rbind(PIP_calibration_WB,PIP_calibration_BB)
PIP_calibration_ancestry<-PIP_calibration_ancestry%>%mutate(Method = fct_recode(Method, "MESuSiE_WB" = "MESuSiE_WB","MESuSiE_BB" = "MESuSiE_BB","SuSiE_WB" = "SuSiE_WB", "SuSiE_BB" = "SuSiE_BB","SuSiE_weighted_WB" = "SuSiE_weighted_WB","SuSiE_weighted_BB" = "SuSiE_weighted_BB","SuSiE_merged_WB" = "SuSiE_merged_WB","SuSiE_merged_BB" = "SuSiE_merged_BB","Paintor_WB" = "Paintor_WB","Paintor_BB" = "Paintor_BB","SuSiEx_WB" = "SuSiEx_WB","SuSiEx_BB" = "SuSiEx_BB","XMAP_WB" = "XMAP_WB","XMAP_BB" = "XMAP_BB", "MultiSuSiE_WB" = "MultiSuSiE_WB", "MultiSuSiE_BB" = "MultiSuSiE_BB","CARMAX_WB" = "CARMAX_WB","CARMAX_BB" = "CARMAX_BB"))%>%mutate(Method = fct_relevel(Method,"MESuSiE_WB","MESuSiE_BB","Paintor_WB", "Paintor_BB","SuSiEx_WB","SuSiEx_BB","XMAP_WB","XMAP_BB","MultiSuSiE_WB", "MultiSuSiE_BB","CARMAX_WB","CARMAX_BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","WB"),paste0("MESuSiE~","BB"),paste0("SuSiE~","WB"),paste0("SuSiE~","BB"),paste0("SuSiE_weighted~","WB"),paste0("SuSiE_weighted~","BB"),paste0("SuSiE_merged~","WB"),paste0("SuSiE_merged~","BB"),paste0("Paintor~","WB"),paste0("Paintor~","BB"),paste0("SuSiEx~","WB"),paste0("SuSiEx~","BB"),paste0("XMAP~","WB"),paste0("XMAP~","BB"),paste0("MultiSuSiE~","WB"),paste0("MultiSuSiE~","BB"),paste0("CARMAX~","WB"),paste0("CARMAX~","BB"))


#######################################
#
#   ROC data 
#
####################################### 
shared_all_ROC_data_dataframe<-ROC_list_to_data(shared_all_ROC_data_list)
either_all_ROC_data_dataframe<-ROC_list_to_data(either_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ROC_list_to_data(ancestry_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ancestry_all_ROC_data_dataframe%>%mutate(Method = fct_recode(Method,
                                                                                              "MESuSiE BB" = "MESuSiE_BB",
                                                                                              "MESuSiE WB" = "MESuSiE_EU",
                                                                                              "SuSiE BB" = "SuSiE_BB",
                                                                                              "SuSiE WB" = "SuSiE_EU",
                                                                                              "SuSiE_weighted BB" = "SuSiE_weighted_BB",
                                                                                              "SuSiE_weighted WB" = "SuSiE_weighted_EU",
                                                                                              "SuSiE_merged BB" = "SuSiE_merged_BB",
                                                                                              "SuSiE_merged WB" = "SuSiE_merged_EU",
                                                                                              "Paintor BB" = "Paintor_BB",
                                                                                              "Paintor WB" = "Paintor_EU",
                                                                                              "SuSiEx WB" = "SuSiEx_EU",
                                                                                              "SuSiEx BB" = "SuSiEx_BB",
                                                                                              "XMAP WB" = "XMAP_EU",
                                                                                              "XMAP BB" = "XMAP_BB",
                                                                                              "MultiSuSiE WB" = "MultiSuSiE_EU",
                                                                                              "MultiSuSiE BB" = "MultiSuSiE_BB",
                                                                                              "CARMAX WB" = "CARMAX_EU",
                                                                                              "CARMAX BB" = "CARMAX_BB",
))


########################################################################
#
#
#FDR Power based on threshold of 0.01,0.05,0.1,0.5
#
#
########################################################################

FDR_Power_shared<-FDR_Power(shared_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
FDR_Power_shared%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))


FDR_Power_either<-FDR_Power(either_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
FDR_Power_either%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

FDR_Power_ancestry<-FDR_Power(ancestry_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method,"MESuSiE BB","MESuSiE WB", "SuSiE BB" ,"SuSiE WB","SuSiE_weighted BB","SuSiE_weighted WB","SuSiE_merged WB","SuSiE_merged BB","Paintor BB" ,"Paintor WB","SuSiEx BB" ,"SuSiEx WB","XMAP WB","XMAP BB","MultiSuSiE WB","MultiSuSiE BB","CARMAX WB", "CARMAX BB"))
FDR_Power_ancestry%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))
#######################################
#
#   time data 
#
####################################### 
#time_all<-data.frame(time_all)
#colnames(time_all)<-c("MESuSiE","SuSiE","Paintor","XMAP","SuSiEx","CARMAX","MultiSuSiE","NUM_SNP")

#######################################
#
#   save the result
#
####################################### 

save(either_all_ROC_data_dataframe,shared_all_ROC_data_dataframe,ancestry_all_ROC_data_dataframe,FDR_Power_either,FDR_Power_shared,FDR_Power_ancestry,PIP_calibration_either,PIP_calibration_Shared,PIP_calibration_ancestry,all_Set_data_dataframe,set_power_summary,data_all,file = paste0(res_dir,res_out))

################################################################
#
#        Numbers to report
#
#
################################################################


##95% Credible set size/power/fdr
all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))
all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))


#Power with a FDR = 0.05 
FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))


#FDR and Power with a PIP threshold of 0.5	
shared_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
ancestry_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
