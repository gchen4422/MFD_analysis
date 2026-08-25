simulation_dir<-"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/"
plot_dir<-paste0(simulation_dir,"Figure/")
system(paste0("mkdir -p ",plot_dir))
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))
res_out<-"shared_50_baseline_updated_meta_correction_updated_Dec_17.RData"
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
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
num_BB = 300000
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
#susiex_runtime = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_result/susiex_runtime_log.txt", header = F)
for(causal_num in 1:3){
  for(h2 in 1:2){
    wrk_dir<-paste0(simulation_dir,"causal_num_",causal_num,"/")
    data_dir<-paste0(wrk_dir,"summary_data/")
    result_dir<-paste0(wrk_dir,"result/")
    
    setwd(result_dir)
    set1 = which(paste0("Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_updated.RData")%in%list.files())
    set2 = which(paste0("Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_correction_updated.RData")%in%list.files())
    replicate_set = Reduce(intersect, list(set1,set2))
    #replicate_set = rep(1:41)
    empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
    
    SuSiE_weighted_ROC        <- SuSiE_weighted_ROC_Either        <- SuSiE_weighted_ROC_EU        <- SuSiE_weighted_ROC_BB        <-
    SuSiE_merged_ROC          <- SuSiE_merged_ROC_Either          <- SuSiE_merged_ROC_EU          <- SuSiE_merged_ROC_BB          <-
    Het_weighted_ROC          <- Het_weighted_ROC_Either          <- Het_weighted_ROC_EU          <- Het_weighted_ROC_BB          <-
    Het_merged_ROC            <- Het_merged_ROC_Either            <- Het_merged_ROC_EU            <- Het_merged_ROC_BB            <-
    SuSiE_rss_weighted_ROC    <- SuSiE_rss_weighted_ROC_Either    <- SuSiE_rss_weighted_ROC_EU    <- SuSiE_rss_weighted_ROC_BB    <-
    SuSiE_rss_merged_ROC      <- SuSiE_rss_merged_ROC_Either      <- SuSiE_rss_merged_ROC_EU      <- SuSiE_rss_merged_ROC_BB      <-
    Slalom_weighted_ROC       <- Slalom_weighted_ROC_Either       <- Slalom_weighted_ROC_EU       <- Slalom_weighted_ROC_BB       <-
    Slalom_merged_ROC         <- Slalom_merged_ROC_Either         <- Slalom_merged_ROC_EU         <- Slalom_merged_ROC_BB         <- 
    empty_mat
    
    
    set_res_all<-c()
    pip_res_all<-c()  
    
    
    
    for(LOCI_num in replicate_set){
      
      paintor_res<-fread(paste0(result_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".mcmc.paintor"))%>%dplyr::rename(SNP = RSID)
      SuSiE_res_name_metal<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_updated.RData")
      load(SuSiE_res_name_metal)
      SuSiE_res_name_metal_correct<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_correction_updated.RData")
      load(SuSiE_res_name_metal_correct)
      
      
      zfile<-read.table(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2),header=T)
      metal_results = fread(paste0(result_dir,"Metal_",causal_num,"_h2_",h2,"_result1.tbl"))
      metal_results_loci <- metal_results %>%
        filter(MarkerName %in% zfile$RSID) %>%
        separate(MarkerName, into = c("CHR","POS"), sep = ":", remove = FALSE, convert = TRUE) %>% arrange(CHR, POS) %>% 
        mutate(Z = qnorm(log(`P-value`/2), lower.tail = FALSE, log.p = TRUE) * sign(Effect))
      
      metal_results_loci_sub = metal_results_loci %>% filter(MarkerName %in% zfile$RSID, HetPVal >0.05)
      
      
      #Format pip and cs for each approach
      
      susie_weighted_cs<-rep(0,length(susie_weighted$pip))
      susie_weighted_cs[unlist(susie_weighted$sets$cs)]<-1
      SuSiE_weighted_res_dataframe<-data.frame(RSID = zfile$RSID,SuSiE_weighted_PIP_Either = susie_weighted$pip,SuSiE_weighted_cs = susie_weighted_cs )
      zfile <- zfile %>%
        left_join(SuSiE_weighted_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      susie_merged_cs<-rep(0,length(susie_merged$pip))
      susie_merged_cs[unlist(susie_merged$sets$cs)]<-1
      SuSiE_merged_res_dataframe<-data.frame(RSID = zfile$RSID,SuSiE_merged_PIP_Either = susie_merged$pip,SuSiE_merged_cs = susie_merged_cs )
      zfile <- zfile %>%
        left_join(SuSiE_merged_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      
      het_weighted_cs<-rep(0,length(susie_weighted_het$pip))
      het_weighted_cs[unlist(susie_weighted_het$sets$cs)]<-1
      het_weighted_res_dataframe<-data.frame(RSID = metal_results_loci_sub$MarkerName,het_weighted_PIP_Either = susie_weighted_het$pip,het_weighted_cs = het_weighted_cs )
      zfile <- zfile %>%
        left_join(het_weighted_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      
      het_merged_cs<-rep(0,length(susie_merged_het$pip))
      het_merged_cs[unlist(susie_merged_het$sets$cs)]<-1
      het_merged_res_dataframe<-data.frame(RSID = metal_results_loci_sub$MarkerName,het_merged_PIP_Either = susie_merged_het$pip,het_merged_cs = het_merged_cs )
      zfile <- zfile %>%
        left_join(het_merged_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      SuSiE_rss_weighted_cs<-rep(0,length(susie_weighted_susie_rss$pip))
      SuSiE_rss_weighted_cs[unlist(susie_weighted_susie_rss$sets$cs)]<-1
      SuSiE_rss_weighted_res_dataframe<-data.frame(RSID = zfile_sub_susie_rss_weighted$RSID,SuSiE_rss_weighted_PIP_Either = susie_weighted_susie_rss$pip,SuSiE_rss_weighted_cs = SuSiE_rss_weighted_cs )
      zfile <- zfile %>%
        left_join(SuSiE_rss_weighted_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)

      
      SuSiE_rss_merged_cs<-rep(0,length(susie_merged_susie_rss$pip))
      SuSiE_rss_merged_cs[unlist(susie_merged_susie_rss$sets$cs)]<-1
      SuSiE_rss_merged_res_dataframe<-data.frame(RSID = zfile_sub_susie_rss_merged$RSID,SuSiE_rss_merged_PIP_Either = susie_merged_susie_rss$pip,SuSiE_rss_merged_cs = SuSiE_rss_merged_cs )
      zfile <- zfile %>%
        left_join(SuSiE_rss_merged_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      
      Slalom_weighted_cs<-rep(0,length(susie_weighted_slalom$pip))
      Slalom_weighted_cs[unlist(susie_weighted_slalom$sets$cs)]<-1
      Slalom_weighted_res_dataframe<-data.frame(RSID = metal_results_loci_sub_slalom$MarkerName,Slalom_weighted_PIP_Either = susie_weighted_slalom$pip,Slalom_weighted_cs = Slalom_weighted_cs )
      zfile <- zfile %>%
        left_join(Slalom_weighted_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      
      Slalom_merged_cs<-rep(0,length(susie_merged_slalom$pip))
      Slalom_merged_cs[unlist(susie_merged_slalom$sets$cs)]<-1
      Slalom_merged_res_dataframe<-data.frame(RSID = metal_results_loci_sub_slalom$MarkerName,Slalom_merged_PIP_Either = susie_merged_slalom$pip,Slalom_merged_cs = Slalom_merged_cs )
      zfile <- zfile %>%
        left_join(Slalom_merged_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)
      
      
      #Set Size/Power
      
      
      SuSiE_weighted_SET      <- zfile %>% summarise(SuSiE_weighted_Power    = sum(SuSiE_weighted_cs != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(SuSiE_weighted_cs))
      SuSiE_merged_SET        <- zfile %>% summarise(SuSiE_merged_Power      = sum(SuSiE_merged_cs   != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(SuSiE_merged_cs))
      Het_weighted_SET        <- zfile %>% summarise(Het_weighted_Power      = sum(het_weighted_cs   != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(het_weighted_cs))
      Het_merged_SET          <- zfile %>% summarise(Het_merged_Power        = sum(het_merged_cs     != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(het_merged_cs))
      SuSiE_rss_weighted_SET  <- zfile %>% summarise(SuSiE_rss_weighted_Power= sum(SuSiE_rss_weighted_cs != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(SuSiE_rss_weighted_cs))
      SuSiE_rss_merged_SET    <- zfile %>% summarise(SuSiE_rss_merged_Power  = sum(SuSiE_rss_merged_cs   != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(SuSiE_rss_merged_cs))
      Slalom_weighted_SET     <- zfile %>% summarise(Slalom_weighted_Power   = sum(Slalom_weighted_cs    != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(Slalom_weighted_cs))
      Slalom_merged_SET       <- zfile %>% summarise(Slalom_merged_Power     = sum(Slalom_merged_cs      != 0 & Signal != 0) / sum(Signal != 0), SuSiE_Size = sum(Slalom_merged_cs))
      

      set_res_all <- rbind(
        set_res_all,
        rbind(
          c(SuSiE_weighted_SET,      "SuSiE_weighted"),
          c(SuSiE_merged_SET,        "SuSiE_merged"),
          c(Het_weighted_SET,        "Het_weighted"),
          c(Het_merged_SET,          "Het_merged"),
          c(SuSiE_rss_weighted_SET,  "SuSiE_rss_weighted"),
          c(SuSiE_rss_merged_SET,    "SuSiE_rss_merged"),
          c(Slalom_weighted_SET,     "Slalom_weighted"),
          c(Slalom_merged_SET,       "Slalom_merged")
        )
      )

            
      
      #ROC for shared/either/ancestry-specific
      
      SuSiE_weighted_ROC      = SuSiE_weighted_ROC      + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_weighted_PIP_Either),      ifelse(zfile$Signal == 3, 1, 0))
      SuSiE_merged_ROC        = SuSiE_merged_ROC        + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_merged_PIP_Either),        ifelse(zfile$Signal == 3, 1, 0))
      Het_weighted_ROC        = Het_weighted_ROC        + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$het_weighted_PIP_Either),        ifelse(zfile$Signal == 3, 1, 0))
      Het_merged_ROC          = Het_merged_ROC          + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$het_merged_PIP_Either),          ifelse(zfile$Signal == 3, 1, 0))
      SuSiE_rss_weighted_ROC  = SuSiE_rss_weighted_ROC  + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_rss_weighted_PIP_Either),  ifelse(zfile$Signal == 3, 1, 0))
      SuSiE_rss_merged_ROC    = SuSiE_rss_merged_ROC    + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_rss_merged_PIP_Either),    ifelse(zfile$Signal == 3, 1, 0))
      Slalom_weighted_ROC     = Slalom_weighted_ROC     + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$Slalom_weighted_PIP_Either),     ifelse(zfile$Signal == 3, 1, 0))
      Slalom_merged_ROC       = Slalom_merged_ROC       + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$Slalom_merged_PIP_Either),       ifelse(zfile$Signal == 3, 1, 0))
      
      
      
      SuSiE_weighted_ROC_Either      = SuSiE_weighted_ROC_Either      + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_weighted_PIP_Either),      ifelse(zfile$Signal != 0, 1, 0))
      SuSiE_merged_ROC_Either        = SuSiE_merged_ROC_Either        + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_merged_PIP_Either),        ifelse(zfile$Signal != 0, 1, 0))
      Het_weighted_ROC_Either        = Het_weighted_ROC_Either        + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$het_weighted_PIP_Either),        ifelse(zfile$Signal != 0, 1, 0))
      Het_merged_ROC_Either          = Het_merged_ROC_Either          + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$het_merged_PIP_Either),          ifelse(zfile$Signal != 0, 1, 0))
      SuSiE_rss_weighted_ROC_Either  = SuSiE_rss_weighted_ROC_Either  + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_rss_weighted_PIP_Either),  ifelse(zfile$Signal != 0, 1, 0))
      SuSiE_rss_merged_ROC_Either    = SuSiE_rss_merged_ROC_Either    + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$SuSiE_rss_merged_PIP_Either),    ifelse(zfile$Signal != 0, 1, 0))
      Slalom_weighted_ROC_Either     = Slalom_weighted_ROC_Either     + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$Slalom_weighted_PIP_Either),     ifelse(zfile$Signal != 0, 1, 0))
      Slalom_merged_ROC_Either       = Slalom_merged_ROC_Either       + compute_ROC(data.frame("SNP" = zfile$RSID, "PIP" = zfile$Slalom_merged_PIP_Either),       ifelse(zfile$Signal != 0, 1, 0))
      
      
      SuSiE_weighted_ROC_EU     = SuSiE_weighted_ROC_EU     + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_weighted_PIP_Either),    ifelse(zfile$Signal == 1, 1, 0))
      SuSiE_weighted_ROC_BB     = SuSiE_weighted_ROC_BB     + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_weighted_PIP_Either),    ifelse(zfile$Signal == 2, 1, 0))
      
      SuSiE_merged_ROC_EU       = SuSiE_merged_ROC_EU       + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_merged_PIP_Either),      ifelse(zfile$Signal == 1, 1, 0))
      SuSiE_merged_ROC_BB       = SuSiE_merged_ROC_BB       + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_merged_PIP_Either),      ifelse(zfile$Signal == 2, 1, 0))
      
      Het_weighted_ROC_EU       = Het_weighted_ROC_EU       + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$het_weighted_PIP_Either),      ifelse(zfile$Signal == 1, 1, 0))
      Het_weighted_ROC_BB       = Het_weighted_ROC_BB       + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$het_weighted_PIP_Either),      ifelse(zfile$Signal == 2, 1, 0))
      
      Het_merged_ROC_EU         = Het_merged_ROC_EU         + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$het_merged_PIP_Either),        ifelse(zfile$Signal == 1, 1, 0))
      Het_merged_ROC_BB         = Het_merged_ROC_BB         + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$het_merged_PIP_Either),        ifelse(zfile$Signal == 2, 1, 0))
      
      SuSiE_rss_weighted_ROC_EU = SuSiE_rss_weighted_ROC_EU + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_rss_weighted_PIP_Either), ifelse(zfile$Signal == 1, 1, 0))
      SuSiE_rss_weighted_ROC_BB = SuSiE_rss_weighted_ROC_BB + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_rss_weighted_PIP_Either), ifelse(zfile$Signal == 2, 1, 0))
      
      SuSiE_rss_merged_ROC_EU   = SuSiE_rss_merged_ROC_EU   + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_rss_merged_PIP_Either),   ifelse(zfile$Signal == 1, 1, 0))
      SuSiE_rss_merged_ROC_BB   = SuSiE_rss_merged_ROC_BB   + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$SuSiE_rss_merged_PIP_Either),   ifelse(zfile$Signal == 2, 1, 0))
      
      Slalom_weighted_ROC_EU    = Slalom_weighted_ROC_EU    + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$Slalom_weighted_PIP_Either),   ifelse(zfile$Signal == 1, 1, 0))
      Slalom_weighted_ROC_BB    = Slalom_weighted_ROC_BB    + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$Slalom_weighted_PIP_Either),   ifelse(zfile$Signal == 2, 1, 0))
      
      Slalom_merged_ROC_EU      = Slalom_merged_ROC_EU      + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$Slalom_merged_PIP_Either),     ifelse(zfile$Signal == 1, 1, 0))
      Slalom_merged_ROC_BB      = Slalom_merged_ROC_BB      + compute_ROC(data.frame(SNP = zfile$RSID, PIP = zfile$Slalom_merged_PIP_Either),     ifelse(zfile$Signal == 2, 1, 0))
      
    
      ##Sumstat
      set_mat <- zfile %>%
        dplyr::transmute(
          SuSiE_weighted_cs       = as.integer(SuSiE_weighted_cs),
          SuSiE_merged_cs         = as.integer(SuSiE_merged_cs),
          Het_weighted_cs         = as.integer(het_weighted_cs),
          Het_merged_cs           = as.integer(het_merged_cs),
          SuSiE_rss_weighted_cs   = as.integer(SuSiE_rss_weighted_cs),
          SuSiE_rss_merged_cs     = as.integer(SuSiE_rss_merged_cs),
          Slalom_weighted_cs      = as.integer(Slalom_weighted_cs),
          Slalom_merged_cs        = as.integer(Slalom_merged_cs)
        )
      
      
      
      locus_sumstat<-zfile%>%mutate(h2 = h2,causal_num  = causal_num,locus = LOCI_num ) %>% dplyr::rename(N_WB = N_1,N_BB = N_2)
      
      data_all<-rbind(data_all,locus_sumstat)
      cat(LOCI_num)
    }
    
    #######################################
    #
    #   credible set combined 
    #
    ####################################### 
    clean_mat <- cbind(
      Power  = as.numeric(set_res_all[, 1]),
      Size   = as.numeric(set_res_all[, 2]),
      Method = c("SuSiE_weighted","SuSiE_merged",
                 "Het_weighted","Het_merged",
                 "SuSiE_rss_weighted","SuSiE_rss_merged",
                 "Slalom_weighted","Slalom_merged")
    )
    
    set_res_all_data <- as.data.frame(clean_mat, stringsAsFactors = FALSE)
    set_res_all_data$Power <- as.numeric(set_res_all_data$Power)
    set_res_all_data$Size  <- as.numeric(set_res_all_data$Size)
    
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
    
    ## ---------------- Either signal part ----------------
    roc_data_list <- list(
      SuSiE_weighted_ROC_Either      = SuSiE_weighted_ROC_Either,
      SuSiE_merged_ROC_Either        = SuSiE_merged_ROC_Either,
      Het_weighted_ROC_Either        = Het_weighted_ROC_Either,
      Het_merged_ROC_Either          = Het_merged_ROC_Either,
      SuSiE_rss_weighted_ROC_Either  = SuSiE_rss_weighted_ROC_Either,
      SuSiE_rss_merged_ROC_Either    = SuSiE_rss_merged_ROC_Either,
      Slalom_weighted_ROC_Either     = Slalom_weighted_ROC_Either,
      Slalom_merged_ROC_Either       = Slalom_merged_ROC_Either
    )
    
    method_names <- c(
      "SuSiE_weighted",
      "SuSiE_merged",
      "Het_weighted",
      "Het_merged",
      "SuSiE_rss_weighted",
      "SuSiE_rss_merged",
      "Slalom_weighted",
      "Slalom_merged"
    )
    
    either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <-
      create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
    
    
    ## ---------------- Shared signal part ----------------
    roc_data_list <- list(
      SuSiE_weighted_ROC      = SuSiE_weighted_ROC,
      SuSiE_merged_ROC        = SuSiE_merged_ROC,
      Het_weighted_ROC        = Het_weighted_ROC,
      Het_merged_ROC          = Het_merged_ROC,
      SuSiE_rss_weighted_ROC  = SuSiE_rss_weighted_ROC,
      SuSiE_rss_merged_ROC    = SuSiE_rss_merged_ROC,
      Slalom_weighted_ROC     = Slalom_weighted_ROC,
      Slalom_merged_ROC       = Slalom_merged_ROC
    )
    
    method_names <- c(
      "SuSiE_weighted",
      "SuSiE_merged",
      "Het_weighted",
      "Het_merged",
      "SuSiE_rss_weighted",
      "SuSiE_rss_merged",
      "Slalom_weighted",
      "Slalom_merged"
    )
    
    shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <-
      create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
    
    
    ## ---------------- Ancestry-specific part ----------------
    roc_data_list <- list(
      SuSiE_weighted_ROC_EU      = SuSiE_weighted_ROC_EU,
      SuSiE_weighted_ROC_BB      = SuSiE_weighted_ROC_BB,
      SuSiE_merged_ROC_EU        = SuSiE_merged_ROC_EU,
      SuSiE_merged_ROC_BB        = SuSiE_merged_ROC_BB,
      Het_weighted_ROC_EU        = Het_weighted_ROC_EU,
      Het_weighted_ROC_BB        = Het_weighted_ROC_BB,
      Het_merged_ROC_EU          = Het_merged_ROC_EU,
      Het_merged_ROC_BB          = Het_merged_ROC_BB,
      SuSiE_rss_weighted_ROC_EU  = SuSiE_rss_weighted_ROC_EU,
      SuSiE_rss_weighted_ROC_BB  = SuSiE_rss_weighted_ROC_BB,
      SuSiE_rss_merged_ROC_EU    = SuSiE_rss_merged_ROC_EU,
      SuSiE_rss_merged_ROC_BB    = SuSiE_rss_merged_ROC_BB,
      Slalom_weighted_ROC_EU     = Slalom_weighted_ROC_EU,
      Slalom_weighted_ROC_BB     = Slalom_weighted_ROC_BB,
      Slalom_merged_ROC_EU       = Slalom_merged_ROC_EU,
      Slalom_merged_ROC_BB       = Slalom_merged_ROC_BB
    )
    
    method_names <- c(
      "SuSiE_weighted_EU",      "SuSiE_weighted_BB",
      "SuSiE_merged_EU",        "SuSiE_merged_BB",
      "Het_weighted_EU",        "Het_weighted_BB",
      "Het_merged_EU",          "Het_merged_BB",
      "SuSiE_rss_weighted_EU",  "SuSiE_rss_weighted_BB",
      "SuSiE_rss_merged_EU",    "SuSiE_rss_merged_BB",
      "Slalom_weighted_EU",     "Slalom_weighted_BB",
      "Slalom_merged_EU",       "Slalom_merged_BB"
    )
    
    ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <-
      create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)
  }
}
#######################################
#
#  Set Size and Power
#
####################################### 

all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
all_Set_data_dataframe <- all_Set_data_dataframe %>% mutate(Method = fct_relevel(Method, "SuSiE_weighted","SuSiE_merged","Het_weighted","Het_merged","SuSiE_rss_weighted","SuSiE_rss_merged","Slalom_weighted","Slalom_merged"))
#levels(all_Set_data_dataframe$Method)<-c("MeSuSiE","Paintor","SuSiE")
cols <- c("Method","h2","causal_num")
set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
colnames(set_power_summary)[4]<-"Power_name"

set_power_summary <- set_power_summary %>% mutate(Method = fct_relevel(Method, "SuSiE_weighted","SuSiE_merged","Het_weighted","Het_merged","SuSiE_rss_weighted","SuSiE_rss_merged","Slalom_weighted","Slalom_merged"))
#levels(set_power_summary$Method)<-c("MESuSiE","SuSiE","Paintor")




#######################################
#
#   ROC data 
#
####################################### 
shared_all_ROC_data_dataframe<-ROC_list_to_data(shared_all_ROC_data_list)
either_all_ROC_data_dataframe<-ROC_list_to_data(either_all_ROC_data_list)
ancestry_all_ROC_data_dataframe<-ROC_list_to_data(ancestry_all_ROC_data_list)
ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%
  mutate(
    Method = fct_recode(
      Method,
      "SuSiE_weighted WB"       = "SuSiE_weighted_EU",
      "SuSiE_weighted BB"       = "SuSiE_weighted_BB",
      "SuSiE_merged WB"         = "SuSiE_merged_EU",
      "SuSiE_merged BB"         = "SuSiE_merged_BB",
      "Het_weighted WB"         = "Het_weighted_EU",
      "Het_weighted BB"         = "Het_weighted_BB",
      "Het_merged WB"           = "Het_merged_EU",
      "Het_merged BB"           = "Het_merged_BB",
      "SuSiE_rss_weighted WB"   = "SuSiE_rss_weighted_EU",
      "SuSiE_rss_weighted BB"   = "SuSiE_rss_weighted_BB",
      "SuSiE_rss_merged WB"     = "SuSiE_rss_merged_EU",
      "SuSiE_rss_merged BB"     = "SuSiE_rss_merged_BB",
      "Slalom_weighted WB"      = "Slalom_weighted_EU",
      "Slalom_weighted BB"      = "Slalom_weighted_BB",
      "Slalom_merged WB"        = "Slalom_merged_EU",
      "Slalom_merged BB"        = "Slalom_merged_BB"
    )
  )


########################################################################
#
#
#FDR Power based on threshold of 0.01,0.05,0.1,0.5
#
#
########################################################################


FDR_Power_shared   <- FDR_Power(shared_all_ROC_data_dataframe)   %>% mutate(Method = fct_relevel(Method, "SuSiE_weighted","SuSiE_merged","Het_weighted","Het_merged","SuSiE_rss_weighted","SuSiE_rss_merged","Slalom_weighted","Slalom_merged"))
FDR_Power_shared   %>% group_by(Method, FDR) %>% summarise(mean_power = mean(Power))

FDR_Power_either   <- FDR_Power(either_all_ROC_data_dataframe)   %>% mutate(Method = fct_relevel(Method, "SuSiE_weighted","SuSiE_merged","Het_weighted","Het_merged","SuSiE_rss_weighted","SuSiE_rss_merged","Slalom_weighted","Slalom_merged"))
FDR_Power_either   %>% group_by(Method, FDR) %>% summarise(mean_power = mean(Power))

FDR_Power_ancestry <- FDR_Power(ancestry_all_ROC_data_dataframe) %>% mutate(Method = fct_relevel(Method, "SuSiE_weighted BB","SuSiE_weighted WB","SuSiE_merged BB","SuSiE_merged WB","Het_weighted BB","Het_weighted WB","Het_merged BB","Het_merged WB","SuSiE_rss_weighted BB","SuSiE_rss_weighted WB","SuSiE_rss_merged BB","SuSiE_rss_merged WB","Slalom_weighted BB","Slalom_weighted WB","Slalom_merged BB","Slalom_merged WB"))
FDR_Power_ancestry %>% group_by(Method, FDR) %>% summarise(mean_power = mean(Power))




#######################################
#
#   save the result
#
####################################### 

save(either_all_ROC_data_dataframe,shared_all_ROC_data_dataframe,ancestry_all_ROC_data_dataframe,FDR_Power_either,FDR_Power_shared,FDR_Power_ancestry,all_Set_data_dataframe,set_power_summary,data_all,file = paste0(res_dir,res_out))

################################################################
#
#        Numbers to report
#
#
################################################################


##95% Credible set size/power/fdr
test =all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2))
all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2))


#Power with a FDR = 0.05 
FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method)%>%summarise(round(mean(Power),2))


#FDR and Power with a PIP threshold of 0.5	
shared_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
ancestry_all_ROC_data_dataframe%>%filter(Cutoff==0.5)%>%group_by(Method)%>%summarise(round(mean(Power),2),round(mean(FDR),2))
