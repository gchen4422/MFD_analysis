
args<-as.numeric(commandArgs(TRUE))
causal_index<-args[1]
External_index<-args[2]
causal_index_name<-c("Both","One")[causal_index]
External_index_name<-c("External_","")[External_index]
simulation_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/",External_index_name,causal_index_name,"/")
plot_dir<-paste0(simulation_dir,"Figure/")
system(paste0("mkdir -p ",plot_dir))
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))
res_out<-paste0("Missing_non_causal_",External_index_name,causal_index_name,"_meta_updated.RData")
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(reshape)
library(data.table)
library(XMAP)
################################################
#
#   Used function 
#
#################################################

compute_ROC<-function(pip_data,true_causal){
  do.call(rbind,lapply(seq(1,0.01,-0.01),function(x){
    
    
    TP = sum(pip_data$SNP[pip_data$PIP>x]%in%pip_data$SNP[true_causal!=0])
    FP = sum(pip_data$PIP>x)-TP
    power = TP/sum(true_causal!=0)
    fdr = FP/(FP+TP)
    return(c(x,ifelse(is.na(power),0,power),ifelse(is.na(fdr),0,fdr)))
    
  }))
}

compute_ROC_uni_susie<-function(pip_data,true_causal){
  do.call(rbind,lapply(seq(1,0.01,-0.01),function(x){
    
    TP = sum(pip_data$SNP[pip_data$PIP_1>x&pip_data$PIP_2<x]%in%pip_data$SNP[true_causal!=0])
    FP = sum(pip_data$PIP_1>x&pip_data$PIP_2<x)-TP
    power = TP/sum(true_causal!=0)
    fdr = FP/(FP+TP)
    return(c(x,ifelse(is.na(power),0,power),ifelse(is.na(fdr),0,fdr)))
    
  }))
}


compute_set_pip<-function(pip_data,coverage = 0.95){
  
  pip_data$PIP = pip_data$PIP/sum(pip_data$PIP)
  snps_selected = pip_data[order(pip_data$PIP,decreasing = T),][which(cumsum(pip_data$PIP[order(pip_data$PIP,decreasing = T)]) < coverage),]
  set_length = nrow(snps_selected)
  total_signal = sum(pip_data$CAUSAL)
  identified_signal =sum(snps_selected$CAUSAL) 
  return(c(identified_signal/total_signal,identified_signal,set_length))
  
}
compute_set_susie<-function(snps_selected,true_causal,coverage = 0.95){
  
  
  set_length = length(snps_selected)
  total_signal = sum(true_causal)
  identified_signal =sum(true_causal[snps_selected]) 
  return(c(identified_signal/total_signal,identified_signal,set_length))
  
}


	create_data_frame <- function(roc_data_list, method_names, replicate_set, h2, causal_num) {
	  # Create data frames for each method, normalize by the number of replicates, and add method names
	  all_data_list <- lapply(1:length(roc_data_list), function(i) {
		data.frame(roc_data_list[[i]] / length(replicate_set), Method = method_names[i])
	  })

	  # Combine data frames into one and set column names
	  all_data_combine <- do.call(rbind, all_data_list)
	  colnames(all_data_combine) <- c("Cutoff", "Power", "FDR", "Method")

	  # Add h2 and causal_num columns
	  all_data_combine$h2 <- h2
	  all_data_combine$causal_num <- causal_num

	  return(all_data_combine)
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
data_all<-c()
for(causal_num in 1:3){
  for(h2 in 1:2){
    wrk_dir<-paste0(simulation_dir,"causal_num_",causal_num,"/")
	data_dir<-paste0(wrk_dir,"summary_data/")
	result_dir<-paste0(wrk_dir,"result/")
	
	setwd(result_dir)
	set1 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
	set2 = which(paste0("Missing_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".mcmc.paintor")%in%list.files())
	set3 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_additional.RData")%in%list.files())
	set4 = which(paste0("Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_updated.RData")%in%list.files())
	
	replicate_set = Reduce(intersect, list(set1,set2,set3,set4))
	
	empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
	MESuSiE_ROC <- SuSiE_ROC <- SuSiE_merged_ROC <-  Paintor_ROC <- susiex_ROC <- XMAP_ROC <- multisusie_ROC <- 
	  carmax_ROC <- MESuSiE_ROC_Either <- SuSiE_ROC_Either <- SuSiE_merged_ROC_Either<-Paintor_ROC_Either <- susiex_ROC_Either <- 
	  XMAP_ROC_Either <- multisusie_ROC_Either <- carmax_ROC_Either  <-
	  MESuSiE_ROC_EU <- MESuSiE_ROC_BB <- SuSiE_ROC_EU <- SuSiE_ROC_BB <- SuSiE_merged_ROC_EU <-SuSiE_merged_ROC_BB <-
	  Paintor_ROC_EU <- Paintor_ROC_BB <- susiex_ROC_EU <- susiex_ROC_BB <- XMAP_ROC_EU <- XMAP_ROC_BB <- 
	  multisusie_ROC_EU <- multisusie_ROC_BB <- carmax_ROC_EU <- carmax_ROC_BB <- 
	  empty_mat
	set_res_all<-c()
	pip_res_all<-c()  
	
for(LOCI_num in replicate_set){

	
	
	
	SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
	load(SuSiE_res_name)
	other_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
	load(other_res_name)
	SuSiE_res_name_metal<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_updated.RData")
	load(SuSiE_res_name_metal)
	XMAP_res_name<-paste0(result_dir,"XMAP_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
	load(XMAP_res_name)
	
	
	#format susiex results
	susiex_res_name<-paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.snp")
	susiex_cred_name<-paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.cs")
	susiex_res = read.table(susiex_res_name, header = T)
	susiex_cs = read.table(susiex_cred_name,header = T)
	#susiex_cs_num = max(susiex_cs$CS_ID)
	susiex_cred = which(zfile_MESuSiE_SuSiE_Paintor$POS %in% susiex_cs$BP)
	
	

	# df: data.frame with per-CS PIP columns for all SNPs
	susiex_pip_cols <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$", names(susiex_res), value = TRUE)
	
	# fast matrix calc with stability (row-wise)
	M <- as.matrix(susiex_res[, susiex_pip_cols, drop = FALSE])
	susiex_res$pip <- 1 - apply(1 - M, 1, prod)  # per row/SNP overall pip
	susiex_res_to_merge = susiex_res %>% dplyr::rename(RSID = "SNP") %>% select(RSID,pip) %>% mutate(SuSiEx_cs = 0) %>% dplyr::rename(SuSiEx_PIP_Either = "pip")
	
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
	  left_join(susiex_res_to_merge, by = "RSID")
	
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor%>% replace(is.na(.), 0)
	zfile_MESuSiE_SuSiE_Paintor$SuSiEx_cs[susiex_cred] <- 1
	
	
	EU_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD1")))
	BB_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD2")))
	
	
	if(causal_index==1){
	  missing_index = which(zfile_MESuSiE_SuSiE_Paintor$EU_missing!=0|zfile_MESuSiE_SuSiE_Paintor$BB_missing!=0)
	  zfile_MESuSiE_SuSiE_Paintor_subset<-zfile_MESuSiE_SuSiE_Paintor%>%filter(EU_missing==0,BB_missing==0)
	  zfile_MESuSiE_SuSiE_Paintor_subset_susie_merged<-zfile_MESuSiE_SuSiE_Paintor_subset
	  EU_cov_subset<-EU_cov[-missing_index,-missing_index]
	  BB_cov_subset<-BB_cov[-missing_index,-missing_index]
	}else if(causal_index==2){
	  missing_index = which(zfile_MESuSiE_SuSiE_Paintor$EU_missing!=0)
	  zfile_MESuSiE_SuSiE_Paintor_subset<-zfile_MESuSiE_SuSiE_Paintor%>%filter(EU_missing==0)
	  zfile_MESuSiE_SuSiE_Paintor_subset_susie_merged<-zfile_MESuSiE_SuSiE_Paintor
	  EU_cov_subset<-EU_cov[-missing_index,-missing_index]
	  BB_cov_subset<-BB_cov[-missing_index,-missing_index]
	}
	
	
	
	susie_merged_cs<-rep(0,length(susie_merged$pip))
	susie_merged_cs[unlist(susie_merged$sets$cs)]<-1
	SuSiE_merged_res_dataframe<-data.frame(RSID = zfile_MESuSiE_SuSiE_Paintor_subset_susie_merged$RSID,SuSiE_merged_PIP_Either = susie_merged$pip,SuSiE_merged_cs = susie_merged_cs )
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
	  left_join(SuSiE_merged_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
	
	
	
	
	
	
	cs1 <- get_CS(xmap, Xcorr = EU_cov_subset, coverage = 0.95, min_abs_corr = 0.5)
	cs2 <- get_CS(xmap, Xcorr = BB_cov_subset, coverage = 0.95, min_abs_corr = 0.5)
	cs_xmap <- unique(unlist(cs1$cs[intersect(names(cs1$cs), names(cs2$cs))],use.names = FALSE))
	pip_xmap <- get_pip(xmap$gamma)
	xmap_res_dataframe<-data.frame(RSID = zfile_MESuSiE_SuSiE_Paintor_subset$RSID,XMAP_PIP_Either = pip_xmap, XMAP_cs = 0)
	xmap_res_dataframe$XMAP_cs[cs_xmap] <- 1
	
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
	  left_join(xmap_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
	
	
	
	
	#format multisusie results
	multisusie_dir = paste0(result_dir,"multisusiex_result/")
	multisusie_snp = fread(paste0(multisusie_dir,"MultiSuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_snp.txt"))
	multisusie_path <- file.path(multisusie_dir, sprintf("MultiSuSiE_CAUSAL_%s_LOCI_%s_h2_%s_output_cs.txt", causal_num, LOCI_num, h2))
	
	multisusie_cs <- tryCatch({
	  if (file.exists(multisusie_path) && file.info(multisusie_path)$size > 0) {
	    dt <- data.table::fread(multisusie_path)
	    if (nrow(dt) > 0 && ncol(dt) > 0) dt else vector()
	  } else vector()
	}, error = function(e) vector())
	
	multisusie_cred <- if (identical(multisusie_cs, logical(0))) {
	  vector()
	} else {
	  which(zfile_MESuSiE_SuSiE_Paintor$POS %in% multisusie_cs$POS)
	}  
	
	multisusie_snp_to_merge = multisusie_snp %>% select(RSID,pip) %>% dplyr::rename(MultiSuSiE_PIP_Either = "pip") %>% mutate(MultiSuSiE_cs = 0)
	
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
	  left_join(multisusie_snp_to_merge, by = "RSID") %>% replace(is.na(.), 0)
	
	zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_cs[multisusie_cred] <- 1
	
	
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
	carmax_res_dataframe<-data.frame(RSID = zfile_MESuSiE_SuSiE_Paintor_subset$RSID,CARMAX_PIP_Either = carmax_pip_either, CARMAX_PIP_WB =carmax_pip_eur ,CARMAX_PIP_BB = carmax_pip_afr ,CARMAX_PIP_Shared =carmax_pip_shared, CARMAX_cs = 0)
	carmax_res_dataframe$CARMAX_cs[carmax_cs] <- 1
	zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
	  left_join(carmax_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
	
	
	
		
	MESuSiE_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(MESuSiE_Power = sum(MESuSiE_cs!=0&Signal!=0)/sum(Signal!=0),MESuSiE_Size = sum(MESuSiE_cs))	
	SuSiE_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(SuSiE_Power = sum(SuSiE_cs!=0&Signal!=0)/sum(Signal!=0),SuSiE_Size = sum(SuSiE_cs))	
	SuSiE_merged_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(SuSiE_merged_Power = sum(SuSiE_merged_cs!=0&Signal!=0)/sum(Signal!=0),SuSiE_Size = sum(susie_merged_cs))	
	Paintor_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(Paintor_Power = sum(Paintor_cs!=0&Signal!=0)/sum(Signal!=0),Paintor_Size = sum(Paintor_cs))
	susiex_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(SuSiEx_Power = sum(SuSiEx_cs!=0&Signal!=0)/sum(Signal!=0),SuSiE_Size = sum(SuSiEx_cs))	
	XMAP_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(XMAP_Power = sum(XMAP_cs!=0&Signal!=0)/sum(Signal!=0),XMAP_Size = sum(XMAP_cs))	
	multisusie_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(MultiSuSiE_Power = sum(MultiSuSiE_cs!=0&Signal!=0)/sum(Signal!=0),MultiSuSiE_Size = sum(MultiSuSiE_cs))	
	carmax_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(CARMAX_Power = sum(CARMAX_cs!=0&Signal!=0)/sum(Signal!=0),SuSiE_Size = sum(CARMAX_cs))	
	
	set_res_all <- rbind(
	  set_res_all,
	  rbind(
	    c(unlist(MESuSiE_SET),   "MESuSiE"),
	    c(unlist(SuSiE_SET),     "SuSiE"),
	    c(unlist(SuSiE_merged_SET),   "SuSiE_merged"),
	    c(unlist(Paintor_SET),   "Paintor"),
	    c(unlist(susiex_SET),    "SuSiEx"),
	    c(unlist(XMAP_SET),      "XMAP"),
	    c(unlist(multisusie_SET),"MultiSuSiE"),
	    c(unlist(carmax_SET),    "CARMAX")
	  )
	)
	
	
	pip_res_all <- zfile_MESuSiE_SuSiE_Paintor %>%
	  mutate(
	    # SuSiE: WB from EU PIP
	    SuSiE_PIP_WB = SuSiE_PIP_EU,
	    
	    # SuSiE merged: only Either
	    SuSiE_merged_PIP_WB      = SuSiE_merged_PIP_Either,
	    SuSiE_merged_PIP_BB      = SuSiE_merged_PIP_Either,
	    SuSiE_merged_PIP_Shared  = SuSiE_merged_PIP_Either,
	    
	    # Paintor: only has Either -> copy to WB/BB/Shared
	    Paintor_PIP_WB     = Paintor_PIP_Either,
	    Paintor_PIP_BB     = Paintor_PIP_Either,
	    Paintor_PIP_Shared = Paintor_PIP_Either,
	    
	    # SuSiEx: only Either
	    SuSiEx_PIP_WB      = SuSiEx_PIP_Either,
	    SuSiEx_PIP_BB      = SuSiEx_PIP_Either,
	    SuSiEx_PIP_Shared  = SuSiEx_PIP_Either,
	    
	    # XMAP: only Either
	    XMAP_PIP_WB        = XMAP_PIP_Either,
	    XMAP_PIP_BB        = XMAP_PIP_Either,
	    XMAP_PIP_Shared    = XMAP_PIP_Either,
	    
	    # MultiSuSiE: only Either
	    MultiSuSiE_PIP_WB      = MultiSuSiE_PIP_Either,
	    MultiSuSiE_PIP_BB      = MultiSuSiE_PIP_Either,
	    MultiSuSiE_PIP_Shared  = MultiSuSiE_PIP_Either
	  ) %>%
	  select(
	    Signal,
	    # Paintor
	    Paintor_PIP_Either, Paintor_PIP_WB, Paintor_PIP_BB, Paintor_PIP_Shared,
	    # MESuSiE
	    MESuSiE_PIP_Either, MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_PIP_Shared,
	    # SuSiE
	    SuSiE_PIP_Either,  SuSiE_PIP_WB, SuSiE_PIP_BB, SuSiE_PIP_Shared,
	    # SuSiE merged
	    SuSiE_merged_PIP_Either,  SuSiE_merged_PIP_WB, SuSiE_merged_PIP_BB, SuSiE_merged_PIP_Shared,
	    # SuSiEx
	    SuSiEx_PIP_Either, SuSiEx_PIP_WB, SuSiEx_PIP_BB, SuSiEx_PIP_Shared,
	    # XMAP
	    XMAP_PIP_Either,   XMAP_PIP_WB, XMAP_PIP_BB, XMAP_PIP_Shared,
	    # MultiSuSiE
	    MultiSuSiE_PIP_Either, MultiSuSiE_PIP_WB, MultiSuSiE_PIP_BB, MultiSuSiE_PIP_Shared,
	    # CARMAX (already has all four)
	    CARMAX_PIP_Either, CARMAX_PIP_WB, CARMAX_PIP_BB, CARMAX_PIP_Shared
	  )
	
	#########Shared part#############
	
	MESuSiE_ROC    = MESuSiE_ROC    + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_Shared),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	SuSiE_ROC      = SuSiE_ROC      + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Shared),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	SuSiE_merged_ROC      = SuSiE_merged_ROC      + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_merged_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	Paintor_ROC    = Paintor_ROC    + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$Paintor_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	susiex_ROC     = susiex_ROC     + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	XMAP_ROC       = XMAP_ROC       + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	multisusie_ROC = multisusie_ROC + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	carmax_ROC     = carmax_ROC     + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_Shared),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
	)
	
	
	
	
	#########Either part#############
	
	MESuSiE_ROC_Either   = MESuSiE_ROC_Either   + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	SuSiE_ROC_Either     = SuSiE_ROC_Either     + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	SuSiE_merged_ROC_Either     = SuSiE_merged_ROC_Either     + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_merged_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	Paintor_ROC_Either   = Paintor_ROC_Either   + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$Paintor_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	susiex_ROC_Either    = susiex_ROC_Either    + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	XMAP_ROC_Either      = XMAP_ROC_Either      + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	multisusie_ROC_Either = multisusie_ROC_Either + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	carmax_ROC_Either    = carmax_ROC_Either    + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
	)
	
	#########Ancestry part#############
	
	MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_WB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
	MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))

	SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC_uni_susie(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP_1" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU ,"PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
  SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC_uni_susie(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP_2" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU ,"PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))  
   	
  SuSiE_merged_ROC_EU =SuSiE_merged_ROC_EU + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_merged_PIP_Either),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
  SuSiE_merged_ROC_BB =SuSiE_merged_ROC_BB + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_merged_PIP_Either),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))
  
	Paintor_ROC_EU = Paintor_ROC_EU + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" =zfile_MESuSiE_SuSiE_Paintor$Paintor_PIP_Either),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
	Paintor_ROC_BB = Paintor_ROC_BB + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" =zfile_MESuSiE_SuSiE_Paintor$Paintor_PIP_Either),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))
	
	
	susiex_ROC_EU = susiex_ROC_EU + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("SuSiEx_PIP_WB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_WB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
	)
	
	susiex_ROC_BB = susiex_ROC_BB + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("SuSiEx_PIP_BB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_BB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
	)
	
	XMAP_ROC_EU = XMAP_ROC_EU + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("XMAP_PIP_WB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_WB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
	)
	
	XMAP_ROC_BB = XMAP_ROC_BB + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("XMAP_PIP_BB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_BB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
	)
	
	multisusie_ROC_EU = multisusie_ROC_EU + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("MultiSuSiE_PIP_WB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_WB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
	)
	
	multisusie_ROC_BB = multisusie_ROC_BB + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("MultiSuSiE_PIP_BB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_BB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$MultiSuSiE_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
	)
	
	carmax_ROC_EU = carmax_ROC_EU + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("CARMAX_PIP_WB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_WB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
	)
	
	carmax_ROC_BB = carmax_ROC_BB + compute_ROC(
	  data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
	             "PIP" = if ("CARMAX_PIP_BB" %in% colnames(zfile_MESuSiE_SuSiE_Paintor))
	               zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_BB
	             else
	               zfile_MESuSiE_SuSiE_Paintor$CARMAX_PIP_Either),
	  ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
	)
	
	
	zfile_MESuSiE_SuSiE_Paintor<-zfile_MESuSiE_SuSiE_Paintor%>%mutate(h2 = h2,causal_num  = causal_num,locus = LOCI_num)
	data_all<-rbind(data_all,zfile_MESuSiE_SuSiE_Paintor)
  cat(LOCI_num)
  
}   
  
 	#######################################
    #
    #   credible set combined 
    #
    ####################################### 
	set_res_all_data = data.frame("Power" = as.numeric(set_res_all[,1]),"Size"=as.numeric(set_res_all[,2]),"Method"=set_res_all[,3])
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
	roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either, SuSiE_ROC_Either = SuSiE_ROC_Either,SuSiE_merged_ROC_Either = SuSiE_merged_ROC_Either, Paintor_ROC_Either = Paintor_ROC_Either,susiex_ROC_Either = susiex_ROC_Either,XMAP_ROC_Either = XMAP_ROC_Either, multisusie_ROC_Either = multisusie_ROC_Either, carmax_ROC_Either= carmax_ROC_Either)
	method_names <- c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")
	either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
	
	# Shared signal part
	roc_data_list <- list(MESuSiE_ROC = MESuSiE_ROC, SuSiE_ROC = SuSiE_ROC,  SuSiE_merged_ROC = SuSiE_merged_ROC,Paintor_ROC = Paintor_ROC,susiex_ROC = susiex_ROC,XMAP_ROC = XMAP_ROC, multisusie_ROC = multisusie_ROC,carmax_ROC= carmax_ROC)
	method_names <- c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")
	shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
	
	# Ancestry specific part
	roc_data_list <- list(MESuSiE_ROC_EU = MESuSiE_ROC_EU,MESuSiE_ROC_BB = MESuSiE_ROC_BB, SuSiE_ROC_EU = SuSiE_ROC_EU, SuSiE_ROC_BB = SuSiE_ROC_BB,SuSiE_merged_ROC_EU = SuSiE_merged_ROC_EU, SuSiE_merged_ROC_BB = SuSiE_merged_ROC_BB,Paintor_ROC_EU = Paintor_ROC_EU,Paintor_ROC_BB = Paintor_ROC_BB, susiex_ROC_EU = susiex_ROC_EU, susiex_ROC_BB = susiex_ROC_BB, XMAP_ROC_EU = XMAP_ROC_EU, XMAP_ROC_BB = XMAP_ROC_BB,multisusie_ROC_EU = multisusie_ROC_EU,multisusie_ROC_BB = multisusie_ROC_BB, carmax_ROC_EU =carmax_ROC_EU, carmax_ROC_BB = carmax_ROC_BB)
	method_names <- c("MESuSiE_EU", "MESuSiE_BB","SuSiE_EU", "SuSiE_BB","SuSiE_merged_EU", "SuSiE_merged_BB","Paintor_EU", "Paintor_BB","SuSiEx_EU","SuSiEx_BB","XMAP_EU","XMAP_BB","MultiSuSiE_EU","MultiSuSiE_BB","CARMAX_EU","CARMAX_BB")
	ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
	
	}
	}
	
	#######################################
    #
    #  Set Size and Power
    #
    ####################################### 

	all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
  all_Set_data_dataframe<-all_Set_data_dataframe%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
  #levels(all_Set_data_dataframe$Method)<-c("MeSuSiE","Paintor","SuSiE")
	cols <- c("Method","h2","causal_num")
	set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
	colnames(set_power_summary)[4]<-"Power_name"

	set_power_summary<-set_power_summary%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
	#levels(set_power_summary$Method)<-c("MESuSiE","SuSiE","Paintor")

    #######################################
    #
    #  PIP calibration
    #
    ####################################### 
	
	####PIP Either
	

	PIP_calibration_either<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_Either,SuSiE_PIP_Either, SuSiE_merged_PIP_Either, Paintor_PIP_Either, SuSiEx_PIP_Either, XMAP_PIP_Either,MultiSuSiE_PIP_Either,CARMAX_PIP_Either),c(1,2,3),c("MESuSiE_PIP_Either","SuSiE_PIP_Either","SuSiE_merged_PIP_Either","Paintor_PIP_Either","SuSiEx_PIP_Either","XMAP_PIP_Either","MultiSuSiE_PIP_Either","CARMAX_PIP_Either"))
	PIP_calibration_either<-PIP_calibration_either%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_PIP_Either","SuSiE" = "SuSiE_PIP_Either","SuSiE_merged" = "SuSiE_merged_PIP_Either","Paintor" = "Paintor_PIP_Either", "SuSiEx" = "SuSiEx_PIP_Either","XMAP" = "XMAP_PIP_Either","MultiSuSiE" = "MultiSuSiE_PIP_Either","CARMAX" = "CARMAX_PIP_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_merged","Paintor", "SuSiEx","XMAP","MultiSuSiE","CARMAX"))
	
	
	####PIP Shared	

	PIP_calibration_Shared<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_Shared,SuSiE_PIP_Shared,SuSiE_merged_PIP_Either, Paintor_PIP_Either, SuSiEx_PIP_Either, XMAP_PIP_Either,MultiSuSiE_PIP_Either,CARMAX_PIP_Shared),c(3),c("MESuSiE_PIP_Shared","SuSiE_PIP_Shared","SuSiE_merged_PIP_Either","Paintor_PIP_Either","SuSiEx_PIP_Either","XMAP_PIP_Either","MultiSuSiE_PIP_Either","CARMAX_PIP_Shared"))
	PIP_calibration_Shared<-PIP_calibration_Shared%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_PIP_Shared","SuSiE" = "SuSiE_PIP_Shared","SuSiE_merged" = "SuSiE_merged_PIP_Either","Paintor" = "Paintor_PIP_Either", "SuSiEx" = "SuSiEx_PIP_Either","XMAP" = "XMAP_PIP_Either","MultiSuSiE" = "MultiSuSiE_PIP_Either","CARMAX" = "CARMAX_PIP_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_merged","Paintor", "SuSiEx","XMAP","MultiSuSiE","CARMAX"))

	
	###PIP ancestry
	data_all<-data_all%>%mutate(SuSiE_merged_WB = SuSiE_merged_PIP_Either,SuSiE_merged_BB = SuSiE_merged_PIP_Either)
	data_all<-data_all%>%mutate(Paintor_WB = Paintor_PIP_Either,Paintor_BB = Paintor_PIP_Either)
	data_all<-data_all%>%mutate(SuSiEx_WB = SuSiEx_PIP_Either,SuSiEx_BB = SuSiEx_PIP_Either)
	data_all<-data_all%>%mutate(XMAP_WB = XMAP_PIP_Either,XMAP_BB = XMAP_PIP_Either)
	data_all<-data_all%>%mutate(MultiSuSiE_WB = MultiSuSiE_PIP_Either,MultiSuSiE_BB = MultiSuSiE_PIP_Either)
	
	
	PIP_calibration_WB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_WB,SuSiE_PIP_EU,SuSiE_merged_WB,Paintor_WB,SuSiEx_WB, XMAP_WB, MultiSuSiE_WB, CARMAX_PIP_WB),c(1),c("MESuSiE_PIP_WB","SuSiE_PIP_EU","SuSiE_merged_WB","Paintor_WB","SuSiEx_WB","XMAP_WB", "MultiSuSiE_WB","CARMAX_PIP_WB"))
	PIP_calibration_BB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_BB,SuSiE_PIP_BB,SuSiE_merged_BB,Paintor_BB,SuSiEx_BB, XMAP_BB, MultiSuSiE_BB, CARMAX_PIP_BB),c(2),c("MESuSiE_PIP_BB","SuSiE_PIP_BB","SuSiE_merged_BB","Paintor_BB","SuSiEx_BB","XMAP_BB", "MultiSuSiE_BB","CARMAX_PIP_BB"))
	
	PIP_calibration_ancestry<-rbind(PIP_calibration_WB,PIP_calibration_BB)
	#PIP_calibration_ancestry<-PIP_calibration_ancestry%>%mutate(Method = fct_recode(Method, "MESuSiE_WB" = "MESuSiE_PIP_WB","MESuSiE_BB" = "MESuSiE_PIP_BB","Paintor_WB" = "Paintor_WB","Paintor_BB" = "Paintor_BB"))%>%mutate(Method = fct_relevel(Method,"MESuSiE_WB","MESuSiE_BB","Paintor_WB","Paintor_BB"))
	#levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","WB"),paste0("MESuSiE~","BB"),paste0("Paintor~","WB"),paste0("Paintor~","BB"))
    
	PIP_calibration_ancestry<-PIP_calibration_ancestry%>%mutate(Method = fct_recode(Method, "MESuSiE_WB" = "MESuSiE_PIP_WB","MESuSiE_BB" = "MESuSiE_PIP_BB","SuSiE_WB" = "SuSiE_PIP_EU", "SuSiE_BB" = "SuSiE_PIP_BB","SuSiE_merged_WB" = "SuSiE_merged_WB", "SuSiE_merged_BB" = "SuSiE_merged_BB","Paintor_WB" = "Paintor_WB","Paintor_BB" = "Paintor_BB","SuSiEx_WB" = "SuSiEx_WB","SuSiEx_BB" = "SuSiEx_BB","XMAP_WB" = "XMAP_WB","XMAP_BB" = "XMAP_BB", "MultiSuSiE_WB" = "MultiSuSiE_WB", "MultiSuSiE_BB" = "MultiSuSiE_BB","CARMAX_WB" = "CARMAX_PIP_WB","CARMAX_BB" = "CARMAX_PIP_BB"))%>%mutate(Method = fct_relevel(Method,"MESuSiE_WB","MESuSiE_BB","SuSiE_WB","SuSiE_BB","SuSiE_merged_WB","SuSiE_merged_BB","Paintor_WB", "Paintor_BB","SuSiEx_WB","SuSiEx_BB","XMAP_WB","XMAP_BB","MultiSuSiE_WB", "MultiSuSiE_BB","CARMAX_WB","CARMAX_BB"))
	levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","WB"),paste0("MESuSiE~","BB"),paste0("SuSiE~","WB"),paste0("SuSiE~","BB"),paste0("SuSiE_merged~","WB"),paste0("SuSiE_merged~","BB"),paste0("Paintor~","WB"),paste0("Paintor~","BB"),paste0("SuSiEx~","WB"),paste0("SuSiEx~","BB"),paste0("XMAP~","WB"),paste0("XMAP~","BB"),paste0("MultiSuSiE~","WB"),paste0("MultiSuSiE~","BB"),paste0("CARMAX~","WB"),paste0("CARMAX~","BB"))
	
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

	FDR_Power_shared<-FDR_Power(shared_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
	FDR_Power_shared%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))
	
	
	FDR_Power_either<-FDR_Power(either_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")))
	FDR_Power_either%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))
	
	FDR_Power_ancestry<-FDR_Power(ancestry_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method,"MESuSiE BB","MESuSiE WB", "SuSiE BB" ,"SuSiE WB","SuSiE_merged WB","SuSiE_merged BB","Paintor BB" ,"Paintor WB","SuSiEx BB" ,"SuSiEx WB","XMAP WB","XMAP BB","MultiSuSiE WB","MultiSuSiE BB","CARMAX WB", "CARMAX BB"))
	FDR_Power_ancestry%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

	################################################################
	#
	#        Numbers to report
	#
	#
	################################################################


	##95% Credible set size/power/fdr
	all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2))
	all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2))


	#Power with a FDR = 0.05 
				FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
				FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
				FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
					
				
	#FDR and Power with a PIP threshold of 0.5	
	shared_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
	ancestry_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
	
	save(either_all_ROC_data_dataframe,shared_all_ROC_data_dataframe,ancestry_all_ROC_data_dataframe,FDR_Power_either,FDR_Power_shared,FDR_Power_ancestry,PIP_calibration_either,PIP_calibration_Shared,PIP_calibration_ancestry,all_Set_data_dataframe,set_power_summary,data_all,file = paste0(res_dir,res_out))

	
