simulation_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/"
plot_dir<-paste0(simulation_dir,"Figure/")
system(paste0("mkdir -p ",plot_dir))
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/res_summary/")
system(paste0("mkdir -p ",res_dir))
res_out<-"shared_50_BB_50000_updated_xmap_meta.RData"
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
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/utility.R")
num_BB = 50000
num_EU = 300000
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
#susiex_runtime = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/susiex_result/susiex_runtime_log.txt", header = F)
for(causal_num in 1:3){
  for(h2 in 1:2){
    wrk_dir<-paste0(simulation_dir,"causal_num_",causal_num,"/")
    data_dir<-paste0(wrk_dir,"summary_data/")
    result_dir<-paste0(wrk_dir,"result/")
    
    setwd(result_dir)
    set1 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
    set2 = which(paste0("CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".mcmc.paintor")%in%list.files())
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
      paintor_res<-fread(paste0(result_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".mcmc.paintor"))%>%dplyr::rename(SNP = RSID)
      SuSiE_res_name_metal<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_updated.RData")
      load(SuSiE_res_name_metal)
      XMAP_res_name<-paste0(result_dir,"XMAP_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      load(XMAP_res_name)
      
      
      
      #format susiex results
      susiex_res_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.snp")
      susiex_cred_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.cs")
      susiex_res = read.table(susiex_res_name, header = T)
      susiex_cs = read.table(susiex_cred_name,header = T)
      #susiex_cs_num = max(susiex_cs$CS_ID)
      susiex_cred = which(paintor_res$POS %in% susiex_cs$BP)
      
      # df: data.frame with per-CS PIP columns for all SNPs
      susiex_pip_cols <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$", names(susiex_res), value = TRUE)
      
      # fast matrix calc with stability (row-wise)
      M <- as.matrix(susiex_res[, susiex_pip_cols, drop = FALSE])
      susiex_res$pip <- 1 - apply(1 - M, 1, prod)  # per row/SNP overall pip
      
      
      # if (susiex_cs_num == 0) {
      #   susiex_res$pip <- 1 / nrow(susiex_res)
      # } else {
      #   # Dynamically create a list of PIP columns for the given number of credible sets
      #   pip_columns <- paste0("PIP.CS", 1:susiex_cs_num,".")
      #   # Check if all the columns exist before proceeding
      #   missing_cols <- pip_columns[!pip_columns %in% colnames(susiex_res)]
      #   if (length(missing_cols) > 0) {
      #     stop(paste("The following columns are missing from the dataset:", paste(missing_cols, collapse = ", ")))
      #   }
      #   # Use pmax with do.call to calculate the element-wise maximum across the selected columns
      #   susiex_res$pip <- do.call(pmax, susiex_res[pip_columns])
      #   #!!!!!!!!!!!!!!
      #   #Need to think how to handle pips from susiex
      # }
      # 
      
      #format xmap results
      
      
      ld_dir = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/causal_num_",causal_num,"/")
      paintor_dir<-paste0(wrk_dir,"summary_data/")
      EU_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD1")))
      BB_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD2")))
      
      cs1 <- get_CS(xmap, Xcorr = EU_cov, coverage = 0.95, min_abs_corr = 0.5)
      cs2 <- get_CS(xmap, Xcorr = BB_cov, coverage = 0.95, min_abs_corr = 0.5)
      cs_xmap <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
      pip_xmap <- get_pip(xmap$gamma)
      
      
      #format multisusie results
      multisusie_dir = "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/multi_susie_result/"
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
        which(paintor_res$POS %in% multisusie_cs$POS)
      }  
      
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
      
      
      
      
      #Set Size/Power
      MESuSiE_SET = compute_set_susie(unique(unlist(MESuSiE_res$cs$cs)),ifelse(paintor_res$Signal!=0,1,0))		
      SuSiE_SET = compute_set_susie(unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs))),ifelse(paintor_res$Signal!=0,1,0))	
      SuSiE_weighted_SET = compute_set_susie(unique(unlist(susie_weighted$sets$cs)),ifelse(paintor_res$Signal!=0,1,0))			
      SuSiE_merged_SET = compute_set_susie(unique(unlist(susie_merged$sets$cs)),ifelse(paintor_res$Signal!=0,1,0))			
      Paintor_SET = compute_set_pip(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob,"CAUSAL"=ifelse(paintor_res$Signal!=0,1,0)))	
      susiex_SET = compute_set_susie(unique(susiex_cred),ifelse(paintor_res$Signal!=0,1,0))	
      XMAP_SET = compute_set_susie(unique(unlist(cs_xmap)),ifelse(paintor_res$Signal!=0,1,0))	
      multisusie_SET = compute_set_susie(unique(multisusie_cred),ifelse(paintor_res$Signal!=0,1,0))	
      carmax_SET = compute_set_susie(unique(carmax_cs),ifelse(paintor_res$Signal!=0,1,0))	
      
      
      
      set_res_all<-rbind(set_res_all,rbind(c(MESuSiE_SET,"MESuSiE"),c( SuSiE_SET,"SuSiE"),c( SuSiE_weighted_SET,"SuSiE_meta_weighted"),c( SuSiE_merged_SET,"SuSiE_meta_merged"),c(Paintor_SET$Paintor_set,"Paintor"),c(susiex_SET,"SuSiEx"),c(XMAP_SET,"XMAP"),c(multisusie_SET,"MultiSuSiE"),c(carmax_SET,"CARMAX")))
      
      
      #ROC for shared/either/ancestry-specific
      MESuSiE_ROC = MESuSiE_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = MESuSiE_res$pip_config[,3]),ifelse(paintor_res$Signal==3,1,0))
      SuSiE_ROC = SuSiE_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = pmin(susie_EU$pip,susie_BB$pip)),ifelse(paintor_res$Signal==3,1,0))
      SuSiE_weighted_ROC = SuSiE_weighted_ROC + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = susie_weighted$pip),ifelse(paintor_res$Signal==3,1,0))
      SuSiE_merged_ROC = SuSiE_merged_ROC + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = susie_merged$pip),ifelse(paintor_res$Signal==3,1,0))
      Paintor_ROC = Paintor_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob),ifelse(paintor_res$Signal==3,1,0))
      susiex_ROC = susiex_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susiex_res$pip),ifelse(paintor_res$Signal==3,1,0))
      XMAP_ROC = XMAP_ROC + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = pip_xmap),ifelse(paintor_res$Signal==3,1,0))
      multisusie_ROC = multisusie_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = multisusie_snp$pip),ifelse(paintor_res$Signal==3,1,0))
      carmax_ROC = carmax_ROC + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = carmax_pip_shared),ifelse(paintor_res$Signal==3,1,0))
      
      
      
      MESuSiE_ROC_Either = MESuSiE_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = MESuSiE_res$pip),ifelse(paintor_res$Signal!=0,1,0))
      SuSiE_ROC_Either = SuSiE_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = pmax(susie_EU$pip,susie_BB$pip)),ifelse(paintor_res$Signal!=0,1,0))
      SuSiE_weighted_ROC_Either = SuSiE_weighted_ROC_Either + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = susie_weighted$pip),ifelse(paintor_res$Signal!=0,1,0))
      SuSiE_merged_ROC_Either = SuSiE_merged_ROC_Either + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = susie_merged$pip),ifelse(paintor_res$Signal!=0,1,0))
      Paintor_ROC_Either = Paintor_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob),ifelse(paintor_res$Signal!=0,1,0))
      susiex_ROC_Either = susiex_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susiex_res$pip),ifelse(paintor_res$Signal!=0,1,0))
      XMAP_ROC_Either = XMAP_ROC_Either + compute_ROC(data.frame("SNP" = paintor_res$SNP ,"PIP" = pip_xmap),ifelse(paintor_res$Signal!=0,1,0))
      multisusie_ROC_Either = multisusie_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = multisusie_snp$pip),ifelse(paintor_res$Signal!=0,1,0))
      carmax_ROC_Either = carmax_ROC_Either + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = carmax_pip_either),ifelse(paintor_res$Signal!=0,1,0))
      
      
      
      
      MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = MESuSiE_res$pip_config[,1]),ifelse(paintor_res$Signal==1,1,0))
      MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = MESuSiE_res$pip_config[,2]),ifelse(paintor_res$Signal==2,1,0))
      
      SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC_uni_susie(data.frame("SNP" =paintor_res$SNP ,"PIP_1" =susie_EU$pip ,"PIP_2" = susie_BB$pip),ifelse(paintor_res$Signal==1,1,0))
      SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC_uni_susie(data.frame("SNP" =paintor_res$SNP ,"PIP_2" =susie_EU$pip ,"PIP_1" = susie_BB$pip),ifelse(paintor_res$Signal==2,1,0))  
      
      SuSiE_weighted_ROC_EU =SuSiE_weighted_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susie_weighted$pip),ifelse(paintor_res$Signal==1,1,0))
      SuSiE_weighted_ROC_BB =SuSiE_weighted_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susie_weighted$pip),ifelse(paintor_res$Signal==2,1,0))  
      
      SuSiE_merged_ROC_EU =SuSiE_merged_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susie_merged$pip),ifelse(paintor_res$Signal==1,1,0))
      SuSiE_merged_ROC_BB =SuSiE_merged_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susie_merged$pip),ifelse(paintor_res$Signal==2,1,0))  
      
      Paintor_ROC_EU = Paintor_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob),ifelse(paintor_res$Signal==1,1,0))
      Paintor_ROC_BB = Paintor_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = paintor_res$Posterior_Prob),ifelse(paintor_res$Signal==2,1,0))
      
      susiex_ROC_EU = susiex_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susiex_res$pip),ifelse(paintor_res$Signal==1,1,0))
      susiex_ROC_BB = susiex_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = susiex_res$pip),ifelse(paintor_res$Signal==2,1,0))
      
      XMAP_ROC_EU = XMAP_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = pip_xmap),ifelse(paintor_res$Signal==1,1,0))
      XMAP_ROC_BB = XMAP_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = pip_xmap),ifelse(paintor_res$Signal==2,1,0))
      
      multisusie_ROC_EU = multisusie_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = multisusie_snp$pip),ifelse(paintor_res$Signal==1,1,0))
      multisusie_ROC_BB = multisusie_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = multisusie_snp$pip),ifelse(paintor_res$Signal==2,1,0))
      
      carmax_ROC_EU = carmax_ROC_EU + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = carmax_pip_eur),ifelse(paintor_res$Signal==1,1,0))
      carmax_ROC_BB = carmax_ROC_BB + compute_ROC(data.frame("SNP" =paintor_res$SNP ,"PIP" = carmax_pip_afr),ifelse(paintor_res$Signal==2,1,0))
      
      
      
      ##Computation time
      time_out_name<-paste0(result_dir,"time_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".txt")
      time_multisusie = fread(paste0(multisusie_dir,"MultiSuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_runtime.txt")) %>% pull(V3)
      
      time_vec<-c(t(fread(time_out_name)),time_multisusie,nrow(paintor_res))
      time_all<-rbind(time_all,time_vec)
      
      ##Sumstat
      set_mat<-data.frame(matrix(0,ncol = 9,nrow = nrow(paintor_res)))
      colnames(set_mat)<-c("MESuSiE_cs","SuSiE_cs","SuSiE_weighted_cs","SuSiE_merged_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","MultiSuSiE_cs","CARMAX_cs")
      set_mat$MESuSiE_cs[unique(unlist(MESuSiE_res$cs$cs))]<-1
      set_mat$SuSiE_cs[unique(union(unlist(susie_EU$sets$cs),unlist(susie_BB$sets$cs)))]<-1
      set_mat$SuSiE_weighted_cs[unique(unlist(susie_weighted$sets$cs))]<-1
      set_mat$SuSiE_merged_cs[unique(unlist(susie_merged$sets$cs))]<-1
      set_mat$Paintor_cs[which(paintor_res$SNP%in%Paintor_SET$Paintor_cs)]<-1
      set_mat$SuSiEx_cs[susiex_cred]<-1
      set_mat$XMAP_cs[unique(unlist(cs_xmap))]<-1
      set_mat$MultiSuSiE_cs[multisusie_cred]<-1
      set_mat$CARMAX_cs[carmax_cs]<-1
      
      
      
      locus_sumstat<-paintor_res%>%mutate(Paintor_cs = set_mat$Paintor_cs,CARMAX_cs = set_mat$CARMAX_cs,CARMAX_Either = carmax_pip_either, CARMAX_Shared = carmax_pip_shared, CARMAX_WB = carmax_pip_eur, CARMAX_BB = carmax_pip_afr,
                                          SuSiEx_cs = set_mat$SuSiEx_cs,SuSiEx_PIP = susiex_res$pip,XMAP_cs = set_mat$XMAP_cs,XMAP_PIP =pip_xmap,MultiSuSIE_cs = set_mat$MultiSuSiE_cs,MultiSuSiE_PIP = multisusie_snp$pip,MESuSiE_Either = MESuSiE_res$pip , MESuSiE_Shared = MESuSiE_res$pip_config[,3],MESuSiE_WB = MESuSiE_res$pip_config[,1],MESuSiE_BB = MESuSiE_res$pip_config[,2],MESuSiE_cs = set_mat$MESuSiE_cs,
                                          SuSiE_Either = pmax(susie_EU$pip,susie_BB$pip),SuSiE_Shared = pmin(susie_EU$pip,susie_BB$pip),SuSiE_WB = susie_EU$pip,SuSiE_BB =susie_BB$pip, SuSiE_cs = set_mat$SuSiE_cs,SuSiE_weighted_PIP = susie_weighted$pip,SuSiE_weighted_cs = set_mat$SuSiE_weighted_cs,SuSiE_merged_PIP = susie_merged$pip,SuSiE_merged_cs = set_mat$SuSiE_merged_cs,
                                          N_WB = num_EU,N_BB = num_BB,h2 = h2,causal_num  = causal_num,locus = LOCI_num )%>%dplyr::rename(Paintor_PIP = Posterior_Prob)
      
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

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/res_summary/")
res_out<-"shared_50_BB_50000_updated_xmap_meta.RData"

load(paste0(res_dir,res_out))


format_supp_table <- function(data, group_col) {
  data %>%
    # 1. Calculate Mean Power
    group_by(across(all_of(group_col)), FDR) %>%
    summarise(mean_power = mean(Power, na.rm = TRUE), .groups = "drop") %>%
    
    # 2. Pivot to Wide Format
    pivot_wider(
      names_from = FDR, 
      values_from = mean_power,
      names_prefix = "Power_at_FDR_" # Adds a prefix so columns aren't just numbers
    ) %>%
    
    # 3. Round numeric columns to 3 decimal places for cleanliness
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

# --- Apply to the datasets ---

# Table 1: "Either" Configuration
supp_table_either <- format_supp_table(FDR_Power_either, "Method")

# Table 2: "Shared" Configuration
supp_table_shared <- format_supp_table(FDR_Power_shared, "Method")

# Table 3: Ancestry Specific (Grouping by Method used in the output)
# Ancestry output combines Method and Ancestry (e.g., "MESuSiE BB")
# If they are separate columns in the original data, use c("Method", "Ancestry")
supp_table_ancestry <- format_supp_table(FDR_Power_ancestry, "Method")

# Replace .RData with .xlsx for the excel output
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
system(paste0("mkdir -p ",res_dir_table))
res_out_xlsx <- str_replace(res_out, ".RData", ".xlsx")
full_path <- paste0(res_dir_table, res_out_xlsx)

sheet_name_raw <- str_remove(res_out_xlsx, ".xlsx")
sheet_name <- substr(sheet_name_raw, 1, 31) 

method_order_base <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                       "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")

# Generate Ancestry-Specific Levels (e.g., "MESuSiE WB", "MESuSiE BB")
method_order_ancestry <- unlist(lapply(method_order_base, function(m) paste0(m, c(" WB", " BB"))))
master_method_levels <- c(method_order_base, method_order_ancestry)
signal_type_order <- c("Either", "Shared", "Ancestry")


# 2. Process and Sort the Data
combined_data <- bind_rows(
  "Either"   = FDR_Power_either,
  "Shared"   = FDR_Power_shared,
  "Ancestry" = FDR_Power_ancestry,
  .id = "Signal_Type" 
) %>%
  select(Method, h2, causal_num, Signal_Type,FDR, Power) %>%
  
  # Apply Factors for ALL sorting columns
  mutate(
    Method = factor(Method, levels = master_method_levels),
    Signal_Type = factor(Signal_Type, levels = signal_type_order) # Enforces Either -> Shared -> Ancestry
  ) %>%
  
  # Order methods for reporting
  arrange(Signal_Type, Method)

# 3. View the result
print(head(combined_data))


wb <- createWorkbook()
addWorksheet(wb, sheet_name)

# Create a bold style for the Title
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# Write the Main Title at Row 1
writeData(wb, sheet_name, "Supplementary Table: FDR Controlled Power Analysis (Combined)", startRow = 1)
addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)

# Write the Combined Data Table starting at Row 2
writeData(wb, sheet_name, combined_data, startRow = 2)

# 5. Save
saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined file to:", full_path))
print(paste("Saved combined table to:", full_path))




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
