
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
res_out<-paste0("Missing_non_causal_",External_index_name,causal_index_name,"_MFD_endpoint_sensitivity.RData")
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
    set_susiex = which(paste0("susiex_result/SuSiEx_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_output_cs95.snp")%in%list.files(recursive = TRUE))
    set_xmap = which(paste0("XMAP_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_additional.RData")%in%list.files())

    replicate_set = Reduce(intersect, list(set1, set_susiex, set_xmap))
    mf_decision = fread(paste0(result_dir, "Finemapping_Method_Decisions_h2_",h2, ".csv"))


    empty_mat<-matrix(0,ncol =3,nrow = length(seq(1,0.01,-0.01)))
    MESuSiE_ROC <- SuSiE_ROC <- MFD_MESuSiE_ROC <- MFD_SuSiEx_ROC <- MFD_XMAP_ROC <-
      MESuSiE_ROC_Either <- SuSiE_ROC_Either <- MFD_MESuSiE_ROC_Either <- MFD_SuSiEx_ROC_Either <- MFD_XMAP_ROC_Either <-
      MESuSiE_ROC_EU <- MESuSiE_ROC_BB <- SuSiE_ROC_EU <- SuSiE_ROC_BB <-
      MFD_MESuSiE_ROC_EU <- MFD_MESuSiE_ROC_BB <-
      MFD_SuSiEx_ROC_EU <- MFD_SuSiEx_ROC_BB <-
      MFD_XMAP_ROC_EU <- MFD_XMAP_ROC_BB <- empty_mat


    set_res_all<-c()
    pip_res_all<-c()

    for(LOCI_num in replicate_set){




      SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      load(SuSiE_res_name)


      ########## Load SuSiEx results ##########
      susiex_res_name<-paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.snp")
      susiex_cred_name<-paste0(result_dir,"susiex_result/","SuSiEx_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_output_cs95.cs")
      susiex_res = read.table(susiex_res_name, header = T)
      susiex_cs = read.table(susiex_cred_name, header = T)
      susiex_cred = which(zfile_MESuSiE_SuSiE_Paintor$POS %in% susiex_cs$BP)

      susiex_pip_cols <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$", names(susiex_res), value = TRUE)
      M <- as.matrix(susiex_res[, susiex_pip_cols, drop = FALSE])
      susiex_res$pip <- 1 - apply(1 - M, 1, prod)
      susiex_res_to_merge = susiex_res %>% dplyr::rename(RSID = "SNP") %>% select(RSID, pip) %>% mutate(SuSiEx_cs = 0) %>% dplyr::rename(SuSiEx_PIP_Either = "pip")

      zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
        left_join(susiex_res_to_merge, by = "RSID")
      zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>% replace(is.na(.), 0)
      zfile_MESuSiE_SuSiE_Paintor$SuSiEx_cs[susiex_cred] <- 1


      ########## Load XMAP results ##########
      XMAP_res_name<-paste0(result_dir,"XMAP_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
      load(XMAP_res_name)

      EU_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD1")))
      BB_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD2")))

      if(causal_index==1){
        missing_index = which(zfile_MESuSiE_SuSiE_Paintor$EU_missing!=0|zfile_MESuSiE_SuSiE_Paintor$BB_missing!=0)
        zfile_MESuSiE_SuSiE_Paintor_subset<-zfile_MESuSiE_SuSiE_Paintor%>%filter(EU_missing==0,BB_missing==0)
        EU_cov_subset<-EU_cov[-missing_index,-missing_index]
        BB_cov_subset<-BB_cov[-missing_index,-missing_index]
      }else if(causal_index==2){
        missing_index = which(zfile_MESuSiE_SuSiE_Paintor$EU_missing!=0)
        zfile_MESuSiE_SuSiE_Paintor_subset<-zfile_MESuSiE_SuSiE_Paintor%>%filter(EU_missing==0)
        EU_cov_subset<-EU_cov[-missing_index,-missing_index]
        BB_cov_subset<-BB_cov[-missing_index,-missing_index]
      }

      cs1 <- get_CS(xmap, Xcorr = EU_cov_subset, coverage = 0.95, min_abs_corr = 0.5)
      cs2 <- get_CS(xmap, Xcorr = BB_cov_subset, coverage = 0.95, min_abs_corr = 0.5)
      cs_xmap <- unique(unlist(cs1$cs[intersect(names(cs1$cs), names(cs2$cs))], use.names = FALSE))
      pip_xmap <- get_pip(xmap$gamma)
      xmap_res_dataframe<-data.frame(RSID = zfile_MESuSiE_SuSiE_Paintor_subset$RSID, XMAP_PIP_Either = pip_xmap, XMAP_cs = 0)
      xmap_res_dataframe$XMAP_cs[cs_xmap] <- 1

      zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor %>%
        left_join(xmap_res_dataframe, by = "RSID") %>% replace(is.na(.), 0)


      ########## Credible Sets ##########

      MESuSiE_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(MESuSiE_Power = sum(MESuSiE_cs!=0&Signal!=0)/sum(Signal!=0),MESuSiE_Size = sum(MESuSiE_cs))
      SuSiE_SET = zfile_MESuSiE_SuSiE_Paintor%>%summarise(SuSiE_Power = sum(SuSiE_cs!=0&Signal!=0)/sum(Signal!=0),SuSiE_Size = sum(SuSiE_cs))

      # MFD_MESuSiE: original MFD logic (MESuSiE vs SuSiE)
      if (mf_decision$Method[LOCI_num] == "MESuSiE") {
        MFD_MESuSiE_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_MESuSiE_Power = sum(MESuSiE_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_MESuSiE_Size  = sum(MESuSiE_cs != 0)
          )
      } else {
        MFD_MESuSiE_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_MESuSiE_Power = sum(SuSiE_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_MESuSiE_Size  = sum(SuSiE_cs != 0)
          )
      }

      # MFD_SuSiEx: SuSiEx as joint endpoint (SuSiEx vs SuSiE)
      if (mf_decision$Method[LOCI_num] == "MESuSiE") {
        MFD_SuSiEx_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_SuSiEx_Power = sum(SuSiEx_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_SuSiEx_Size  = sum(SuSiEx_cs != 0)
          )
      } else {
        MFD_SuSiEx_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_SuSiEx_Power = sum(SuSiE_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_SuSiEx_Size  = sum(SuSiE_cs != 0)
          )
      }

      # MFD_XMAP: XMAP as joint endpoint (XMAP vs SuSiE)
      if (mf_decision$Method[LOCI_num] == "MESuSiE") {
        MFD_XMAP_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_XMAP_Power = sum(XMAP_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_XMAP_Size  = sum(XMAP_cs != 0)
          )
      } else {
        MFD_XMAP_SET <- zfile_MESuSiE_SuSiE_Paintor %>%
          summarise(
            MFD_XMAP_Power = sum(SuSiE_cs != 0 & Signal != 0) / sum(Signal != 0),
            MFD_XMAP_Size  = sum(SuSiE_cs != 0)
          )
      }

      set_res_all <- rbind(
        set_res_all,
        rbind(
          c(unlist(MESuSiE_SET),       "MESuSiE"),
          c(unlist(SuSiE_SET),         "SuSiE"),
          c(unlist(MFD_MESuSiE_SET),   "MFD_MESuSiE"),
          c(unlist(MFD_SuSiEx_SET),    "MFD_SuSiEx"),
          c(unlist(MFD_XMAP_SET),      "MFD_XMAP")
        )
      )


      ########## PIP columns ##########
      is_mesusie <- mf_decision$Method[LOCI_num] == "MESuSiE"

      zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>%
        mutate(
          SuSiE_PIP_WB = SuSiE_PIP_EU,

          # MFD_MESuSiE PIP columns
          MFD_MESuSiE_PIP_Either = if (is_mesusie) MESuSiE_PIP_Either else SuSiE_PIP_Either,
          MFD_MESuSiE_PIP_WB     = if (is_mesusie) MESuSiE_PIP_WB     else SuSiE_PIP_WB,
          MFD_MESuSiE_PIP_BB     = if (is_mesusie) MESuSiE_PIP_BB     else SuSiE_PIP_BB,
          MFD_MESuSiE_PIP_Shared = if (is_mesusie) MESuSiE_PIP_Shared else SuSiE_PIP_Shared,

          # MFD_SuSiEx PIP columns (SuSiEx only has Either, use it for all)
          MFD_SuSiEx_PIP_Either = if (is_mesusie) SuSiEx_PIP_Either else SuSiE_PIP_Either,
          MFD_SuSiEx_PIP_WB     = if (is_mesusie) SuSiEx_PIP_Either else SuSiE_PIP_WB,
          MFD_SuSiEx_PIP_BB     = if (is_mesusie) SuSiEx_PIP_Either else SuSiE_PIP_BB,
          MFD_SuSiEx_PIP_Shared = if (is_mesusie) SuSiEx_PIP_Either else SuSiE_PIP_Shared,

          # MFD_XMAP PIP columns (XMAP only has Either, use it for all)
          MFD_XMAP_PIP_Either = if (is_mesusie) XMAP_PIP_Either else SuSiE_PIP_Either,
          MFD_XMAP_PIP_WB     = if (is_mesusie) XMAP_PIP_Either else SuSiE_PIP_WB,
          MFD_XMAP_PIP_BB     = if (is_mesusie) XMAP_PIP_Either else SuSiE_PIP_BB,
          MFD_XMAP_PIP_Shared = if (is_mesusie) XMAP_PIP_Either else SuSiE_PIP_Shared
        )



      pip_res_all <- zfile_MESuSiE_SuSiE_Paintor %>%
        select(
          Signal,
          # MESuSiE
          MESuSiE_PIP_Either, MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_PIP_Shared,
          # SuSiE
          SuSiE_PIP_Either, SuSiE_PIP_WB, SuSiE_PIP_BB, SuSiE_PIP_Shared,
          # MFD_MESuSiE
          MFD_MESuSiE_PIP_Either, MFD_MESuSiE_PIP_WB, MFD_MESuSiE_PIP_BB, MFD_MESuSiE_PIP_Shared,
          # MFD_SuSiEx
          MFD_SuSiEx_PIP_Either, MFD_SuSiEx_PIP_WB, MFD_SuSiEx_PIP_BB, MFD_SuSiEx_PIP_Shared,
          # MFD_XMAP
          MFD_XMAP_PIP_Either, MFD_XMAP_PIP_WB, MFD_XMAP_PIP_BB, MFD_XMAP_PIP_Shared
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

      if (is_mesusie) {
        MFD_MESuSiE_ROC = MFD_MESuSiE_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_Shared),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
        MFD_SuSiEx_ROC = MFD_SuSiEx_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
        MFD_XMAP_ROC = MFD_XMAP_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
      } else {
        MFD_MESuSiE_ROC = MFD_MESuSiE_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Shared),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
        MFD_SuSiEx_ROC = MFD_SuSiEx_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Shared),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
        MFD_XMAP_ROC = MFD_XMAP_ROC + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Shared),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 3, 1, 0)
        )
      }




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

      if (is_mesusie) {
        MFD_MESuSiE_ROC_Either = MFD_MESuSiE_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
        MFD_SuSiEx_ROC_Either = MFD_SuSiEx_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
        MFD_XMAP_ROC_Either = MFD_XMAP_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
      } else {
        MFD_MESuSiE_ROC_Either = MFD_MESuSiE_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
        MFD_SuSiEx_ROC_Either = MFD_SuSiEx_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
        MFD_XMAP_ROC_Either = MFD_XMAP_ROC_Either + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal != 0, 1, 0)
        )
      }


      #########Ancestry part#############

      MESuSiE_ROC_EU = MESuSiE_ROC_EU + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_WB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
      MESuSiE_ROC_BB = MESuSiE_ROC_BB + compute_ROC(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))

      SuSiE_ROC_EU =SuSiE_ROC_EU + compute_ROC_uni_susie(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP_1" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU ,"PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==1,1,0))
      SuSiE_ROC_BB =SuSiE_ROC_BB + compute_ROC_uni_susie(data.frame("SNP" =zfile_MESuSiE_SuSiE_Paintor$RSID ,"PIP_2" =zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU ,"PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal==2,1,0))

      if (is_mesusie) {

        # --- MFD_MESuSiE: EU/BB Decision: MESuSiE ---
        MFD_MESuSiE_ROC_EU = MFD_MESuSiE_ROC_EU + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_WB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_MESuSiE_ROC_BB = MFD_MESuSiE_ROC_BB + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$MESuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

        # --- MFD_SuSiEx: EU/BB Decision: SuSiEx (only has Either) ---
        MFD_SuSiEx_ROC_EU = MFD_SuSiEx_ROC_EU + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_SuSiEx_ROC_BB = MFD_SuSiEx_ROC_BB + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$SuSiEx_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

        # --- MFD_XMAP: EU/BB Decision: XMAP (only has Either) ---
        MFD_XMAP_ROC_EU = MFD_XMAP_ROC_EU + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_XMAP_ROC_BB = MFD_XMAP_ROC_BB + compute_ROC(
          data.frame("SNP" = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP" = zfile_MESuSiE_SuSiE_Paintor$XMAP_PIP_Either),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

      } else {

        # --- MFD_MESuSiE: EU/BB Decision: SuSiE Post-Hoc ---
        MFD_MESuSiE_ROC_EU = MFD_MESuSiE_ROC_EU + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_MESuSiE_ROC_BB = MFD_MESuSiE_ROC_BB + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

        # --- MFD_SuSiEx: EU/BB Decision: SuSiE Post-Hoc ---
        MFD_SuSiEx_ROC_EU = MFD_SuSiEx_ROC_EU + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_SuSiEx_ROC_BB = MFD_SuSiEx_ROC_BB + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

        # --- MFD_XMAP: EU/BB Decision: SuSiE Post-Hoc ---
        MFD_XMAP_ROC_EU = MFD_XMAP_ROC_EU + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 1, 1, 0)
        )
        MFD_XMAP_ROC_BB = MFD_XMAP_ROC_BB + compute_ROC_uni_susie(
          data.frame("SNP"   = zfile_MESuSiE_SuSiE_Paintor$RSID,
                     "PIP_2" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_EU,
                     "PIP_1" = zfile_MESuSiE_SuSiE_Paintor$SuSiE_PIP_BB),
          ifelse(zfile_MESuSiE_SuSiE_Paintor$Signal == 2, 1, 0)
        )

      }


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
    roc_data_list <- list(MESuSiE_ROC_Either = MESuSiE_ROC_Either, SuSiE_ROC_Either = SuSiE_ROC_Either, MFD_MESuSiE_ROC_Either = MFD_MESuSiE_ROC_Either, MFD_SuSiEx_ROC_Either = MFD_SuSiEx_ROC_Either, MFD_XMAP_ROC_Either = MFD_XMAP_ROC_Either)
    method_names <- c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")
    either_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)

    # Shared signal part
    roc_data_list <- list(MESuSiE_ROC = MESuSiE_ROC, SuSiE_ROC = SuSiE_ROC, MFD_MESuSiE_ROC = MFD_MESuSiE_ROC, MFD_SuSiEx_ROC = MFD_SuSiEx_ROC, MFD_XMAP_ROC = MFD_XMAP_ROC)
    method_names <- c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")
    shared_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)

    # Ancestry specific part
    roc_data_list <- list(
      MESuSiE_ROC_EU = MESuSiE_ROC_EU, MESuSiE_ROC_BB = MESuSiE_ROC_BB,
      SuSiE_ROC_EU = SuSiE_ROC_EU, SuSiE_ROC_BB = SuSiE_ROC_BB,
      MFD_MESuSiE_ROC_EU = MFD_MESuSiE_ROC_EU, MFD_MESuSiE_ROC_BB = MFD_MESuSiE_ROC_BB,
      MFD_SuSiEx_ROC_EU = MFD_SuSiEx_ROC_EU, MFD_SuSiEx_ROC_BB = MFD_SuSiEx_ROC_BB,
      MFD_XMAP_ROC_EU = MFD_XMAP_ROC_EU, MFD_XMAP_ROC_BB = MFD_XMAP_ROC_BB
    )
    method_names <- c("MESuSiE_EU", "MESuSiE_BB", "SuSiE_EU", "SuSiE_BB",
                      "MFD_MESuSiE_EU", "MFD_MESuSiE_BB",
                      "MFD_SuSiEx_EU", "MFD_SuSiEx_BB",
                      "MFD_XMAP_EU", "MFD_XMAP_BB")
    ancestry_all_ROC_data_list[[paste0("causal_", causal_num, "_h2_", h2)]] <- create_data_frame(roc_data_list, method_names, replicate_set, h2, causal_num)

  }
}


#######################################
#
#  Set Size and Power
#
#######################################

all_Set_data_dataframe<-ROC_list_to_data(all_Set_data_list)
all_Set_data_dataframe<-all_Set_data_dataframe%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")))
cols <- c("Method","h2","causal_num")
set_power_summary<-data.frame(all_Set_data_dataframe %>% group_by(across(all_of(cols))) %>% summarize_at(vars(Power),list(name = mean)))
colnames(set_power_summary)[4]<-"Power_name"

set_power_summary<-set_power_summary%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")))

#######################################
#
#  PIP calibration
#
#######################################

####PIP Either


PIP_calibration_either<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_Either,SuSiE_PIP_Either, MFD_MESuSiE_PIP_Either, MFD_SuSiEx_PIP_Either, MFD_XMAP_PIP_Either),c(1,2,3),c("MESuSiE_PIP_Either","SuSiE_PIP_Either","MFD_MESuSiE_PIP_Either","MFD_SuSiEx_PIP_Either","MFD_XMAP_PIP_Either"))
PIP_calibration_either<-PIP_calibration_either%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_PIP_Either","SuSiE" = "SuSiE_PIP_Either","MFD_MESuSiE" = "MFD_MESuSiE_PIP_Either","MFD_SuSiEx" = "MFD_SuSiEx_PIP_Either","MFD_XMAP" = "MFD_XMAP_PIP_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP"))


####PIP Shared

PIP_calibration_Shared<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_Shared,SuSiE_PIP_Shared,MFD_MESuSiE_PIP_Shared,MFD_SuSiEx_PIP_Shared,MFD_XMAP_PIP_Shared),c(3),c("MESuSiE_PIP_Shared","SuSiE_PIP_Shared","MFD_MESuSiE_PIP_Shared","MFD_SuSiEx_PIP_Shared","MFD_XMAP_PIP_Shared"))
PIP_calibration_Shared<-PIP_calibration_Shared%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_PIP_Shared","SuSiE" = "SuSiE_PIP_Shared","MFD_MESuSiE" = "MFD_MESuSiE_PIP_Shared","MFD_SuSiEx" = "MFD_SuSiEx_PIP_Shared","MFD_XMAP" = "MFD_XMAP_PIP_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP"))


###PIP ancestry

PIP_calibration_WB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_WB,SuSiE_PIP_EU,MFD_MESuSiE_PIP_WB,MFD_SuSiEx_PIP_WB,MFD_XMAP_PIP_WB),c(1),c("MESuSiE_PIP_WB","SuSiE_PIP_EU","MFD_MESuSiE_PIP_WB","MFD_SuSiEx_PIP_WB","MFD_XMAP_PIP_WB"))
PIP_calibration_BB<-create_obs_frq(data_all%>%select(Signal,h2,causal_num,MESuSiE_PIP_BB,SuSiE_PIP_BB,MFD_MESuSiE_PIP_BB,MFD_SuSiEx_PIP_BB,MFD_XMAP_PIP_BB),c(2),c("MESuSiE_PIP_BB","SuSiE_PIP_BB","MFD_MESuSiE_PIP_BB","MFD_SuSiEx_PIP_BB","MFD_XMAP_PIP_BB"))

PIP_calibration_ancestry<-rbind(PIP_calibration_WB,PIP_calibration_BB)

PIP_calibration_ancestry<-PIP_calibration_ancestry%>%mutate(Method = fct_recode(Method, "MESuSiE_WB" = "MESuSiE_PIP_WB","MESuSiE_BB" = "MESuSiE_PIP_BB","SuSiE_WB" = "SuSiE_PIP_EU", "SuSiE_BB" = "SuSiE_PIP_BB","MFD_MESuSiE_WB" = "MFD_MESuSiE_PIP_WB","MFD_MESuSiE_BB" = "MFD_MESuSiE_PIP_BB","MFD_SuSiEx_WB" = "MFD_SuSiEx_PIP_WB","MFD_SuSiEx_BB" = "MFD_SuSiEx_PIP_BB","MFD_XMAP_WB" = "MFD_XMAP_PIP_WB","MFD_XMAP_BB" = "MFD_XMAP_PIP_BB"))%>%mutate(Method = fct_relevel(Method,"MESuSiE_WB","MESuSiE_BB","SuSiE_WB","SuSiE_BB","MFD_MESuSiE_WB","MFD_MESuSiE_BB","MFD_SuSiEx_WB","MFD_SuSiEx_BB","MFD_XMAP_WB","MFD_XMAP_BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","WB"),paste0("MESuSiE~","BB"),paste0("SuSiE~","WB"),paste0("SuSiE~","BB"),paste0("MFD_MESuSiE~","WB"),paste0("MFD_MESuSiE~","BB"),paste0("MFD_SuSiEx~","WB"),paste0("MFD_SuSiEx~","BB"),paste0("MFD_XMAP~","WB"),paste0("MFD_XMAP~","BB"))

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
                                                                                              "MFD_MESuSiE BB" = "MFD_MESuSiE_BB",
                                                                                              "MFD_MESuSiE WB" = "MFD_MESuSiE_EU",
                                                                                              "MFD_SuSiEx BB" = "MFD_SuSiEx_BB",
                                                                                              "MFD_SuSiEx WB" = "MFD_SuSiEx_EU",
                                                                                              "MFD_XMAP BB" = "MFD_XMAP_BB",
                                                                                              "MFD_XMAP WB" = "MFD_XMAP_EU"
))

########################################################################
#
#
#FDR Power based on threshold of 0.01,0.05,0.1,0.5
#
#
########################################################################

FDR_Power_shared<-FDR_Power(shared_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")))
FDR_Power_shared%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))


FDR_Power_either<-FDR_Power(either_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method, c("MESuSiE", "SuSiE", "MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")))
FDR_Power_either%>%group_by(Method,FDR)%>%summarise(mean_power= mean(Power))

FDR_Power_ancestry<-FDR_Power(ancestry_all_ROC_data_dataframe)%>%mutate(Method = fct_relevel(Method,"MESuSiE BB","MESuSiE WB", "SuSiE BB" ,"SuSiE WB","MFD_MESuSiE WB","MFD_MESuSiE BB","MFD_SuSiEx WB","MFD_SuSiEx BB","MFD_XMAP WB","MFD_XMAP BB"))
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
