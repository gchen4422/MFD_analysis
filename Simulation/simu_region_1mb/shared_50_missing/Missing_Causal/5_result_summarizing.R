args<-as.numeric(commandArgs(TRUE))
causal_index = args[1]
External_index = args[2]
causal_index_name = c("Both","One")[causal_index]
External_index_name = c("","External_")[External_index]
simulation_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/")
plot_dir<-paste0(simulation_dir,"Figure/")
system(paste0("mkdir -p ",plot_dir))
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))

res_out_name<-paste0("Missing_causal_",External_index_name,causal_index_name,"_updated")
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
library(tidyverse)
library(tidyr)
library(reshape)
library(data.table)
library(XMAP)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
################################################################
#
#        Result Summarizing
#
#
################################################################
all_ld_true_signal<-c()
data_all<-c()
for(causal_num in 1){
  for(h2 in 1:2){
    wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",causal_num,"/")	
	  data_dir<-paste0(wrk_dir,"summary_data/")
	  data_dir_ld = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/",External_index_name,causal_index_name,"/causal_num_",causal_num,"/summary_data/")
    result_dir<-paste0(wrk_dir,"result/")    
    setwd(result_dir)
    
    
    set1 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
    set2 = which(paste0("Missing_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".mcmc.paintor")%in%list.files())
    set3 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_additional.RData")%in%list.files())
    set4 = which(paste0("Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_updated.RData")%in%list.files())
    
    
    replicate_set = Reduce(intersect, list(set1,set2,set3,set4))
    
    ld_true_signal<-c()
    data_all_h2<-c()
    
    for(LOCI_num in replicate_set){
       
		SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
		load(SuSiE_res_name)
		other_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
		load(other_res_name)
		xmap_res_name<-paste0(result_dir,"XMAP_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
		load(xmap_res_name)
		SuSiE_res_name_metal<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_updated.RData")
		load(SuSiE_res_name_metal)
		
		
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
		
		
		EU_cov<-as.matrix(fread(paste0(data_dir_ld,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD1")))
		BB_cov<-as.matrix(fread(paste0(data_dir_ld,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD2")))
		
		
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
		
		
		

		zfile_MESuSiE_SuSiE_Paintor$Region = LOCI_num
		data_all_h2<-rbind(data_all_h2,zfile_MESuSiE_SuSiE_Paintor)

		which(zfile_MESuSiE_SuSiE_Paintor$Signal!=0)
		true_signal_pos<-which(zfile_MESuSiE_SuSiE_Paintor$Signal!=0)

		SuSiE_signal_pos_either<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(pmax(SuSiE_PIP_BB,SuSiE_PIP_EU)))%>%pull(max_index)
		SuSiE_signal_pos_share<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(pmin(SuSiE_PIP_BB,SuSiE_PIP_EU)))%>%pull(max_index)

		MESuSiE_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(MESuSiE_PIP_Shared))%>%pull(max_index)
		SuSiE_merged_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(SuSiE_merged_PIP_Either))%>%pull(max_index)
		Paintor_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(Paintor_PIP_Either))%>%pull(max_index)
		SuSiEx_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(SuSiEx_PIP_Either))%>%pull(max_index)
		XMAP_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(XMAP_PIP_Either))%>%pull(max_index)
		MultiSuSiE_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(MultiSuSiE_PIP_Either))%>%pull(max_index)
		##For CARMAX think again which pip should be used
		CARMAX_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(CARMAX_PIP_Either))%>%pull(max_index) 

		#EU_cov<-as.matrix(fread(paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/LD_data/Region_",gene_num,".LD1")))	
		#BB_cov<-as.matrix(fread(paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/LD_data/Region_",gene_num,".LD2")))
		max_cov<-pmax(abs(EU_cov),abs(BB_cov))
		ld_vec<-c(max_cov[true_signal_pos,c(MESuSiE_signal_pos,SuSiE_signal_pos_either,SuSiE_signal_pos_share,SuSiE_merged_signal_pos,Paintor_signal_pos,SuSiEx_signal_pos,XMAP_signal_pos,MultiSuSiE_signal_pos,CARMAX_signal_pos)])
		ld_true_signal<-rbind(ld_true_signal,ld_vec)
		cat(LOCI_num)
		
    }   
   
    ld_true_signal<-data.frame(ld_true_signal)
    colnames(ld_true_signal)<-c("MESuSiE","SuSiE_Either","SuSiE_Share","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")
    ld_true_signal$h2 = h2    
    all_ld_true_signal<-rbind(all_ld_true_signal, ld_true_signal)
  
    data_all_h2<-data.frame(data_all_h2)
    data_all_h2$h2 = h2    
    data_all<-rbind(data_all,data_all_h2)
  }
}

data_all_ori = data_all
####################################################################################
#
#
#		 Check the PIP threshold of 0.5
#
#
#
####################################################################################

###Check Number of Signals and Potentially Compute the proprotion of shared and ancestry-specific ones
data_all%>%summarise(sum(MESuSiE_PIP_Shared>0.5),sum(MESuSiE_PIP_WB>0.5),sum(MESuSiE_PIP_BB>0.5))/200
data_all%>%summarise(sum(SuSiE_PIP_Shared>0.5),sum(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5),sum(SuSiE_PIP_BB>0.5&SuSiE_PIP_EU<0.5))/200
data_all%>%summarise(sum(Paintor_PIP_Either>0.5))/200
####################################################################################
#
#
#			Distance
#
#
#
####################################################################################
###Signal Distance and LD (MESuSiE, SuSiE Either, Paintor)
Signal_POS<-data_all%>%filter(Signal!=0)%>%mutate(signal_pos = POS)%>%select(h2,Region,signal_pos)

MESuSiE_distance<-data_all%>%filter(MESuSiE_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "MESuSiE")
SuSiE_distance<-data_all%>%filter(SuSiE_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "SuSiE")
SuSiE_merged_distance<-data_all%>%filter(SuSiE_merged_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "SuSiE_merged")
Paintor_distance<-data_all%>%filter(Paintor_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "Paintor")
SuSiEx_distance<-data_all%>%filter(SuSiEx_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "SuSiEx")
XMAP_distance<-data_all%>%filter(XMAP_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "XMAP")
MultiSuSiE_distance<-data_all%>%filter(MultiSuSiE_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "MultiSuSiE")
CARMAX_distance<-data_all%>%filter(CARMAX_PIP_Either>0.5)%>% left_join(Signal_POS, by = c("h2","Region"))%>%mutate(distance_to_signal =abs(signal_pos-POS))%>%group_by(h2)%>%summarise( Distance = mean(distance_to_signal))%>%mutate(Method = "CARMAX")





Signal_distance_data<-rbind(MESuSiE_distance,SuSiE_distance,SuSiE_merged_distance,Paintor_distance,SuSiEx_distance,XMAP_distance,MultiSuSiE_distance,CARMAX_distance)

Signal_distance_data$h2<-factor(Signal_distance_data$h2)
levels(Signal_distance_data$h2)=c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
Signal_distance_data$Method<-factor(Signal_distance_data$Method)

Signal_distance_data <- Signal_distance_data %>%
  mutate(Method = fct_relevel(Method,
                             "MESuSiE" ,
                             "SuSiE",
                             "SuSiE_merged",
                             "Paintor",
                             "SuSiEx",
                             "XMAP",
                             "MultiSuSiE",
                             "CARMAX"))

Signal_distance_data%>%group_by(Method)%>%summarise(mean(Distance)/1024)	
						 
p_signal_distance = ggplot(data=Signal_distance_data,aes(x = Method, y = Distance,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e"))+ylab("Distance to true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
p_signal_distance = p_signal_distance + theme(axis.text.x = element_blank(),
                                  axis.text.y = element_text( size = 18),  
                                  axis.title.x = element_text( size = 18,face="bold"),
                                  axis.title.y = element_text( size = 18,face="bold"),
                                  strip.text.x = element_text(size = 18),
                                  strip.text.y= element_text(size = 18),
                                  legend.text=element_text(size=18),
                                  legend.title=element_text(size=22,face="bold"),
                                  plot.title = element_text(size=16,hjust = 0.5)) 
####################################################################################
#
#
#			LD
#
#
#
####################################################################################
all_ld_true_signal$Region = rep(seq_len(100),2)
all_ld_true_signal<-data.frame(all_ld_true_signal)

MESuSiE_LD<-data_all%>%filter(MESuSiE_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(MESuSiE), Method = "MESuSiE" )%>%select(h2,Method,Cor )
SuSiE_LD<-data_all%>%filter(SuSiE_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(SuSiE_Either), Method = "SuSiE" )%>%select(h2,Method,Cor )
SuSiE_merged_LD<-data_all%>%filter(SuSiE_merged_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(SuSiE_merged), Method = "SuSiE_merged" )%>%select(h2,Method,Cor )
Paintor_LD<-data_all%>%filter(Paintor_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(Paintor), Method = "Paintor" )%>%select(h2,Method,Cor )
SuSiEx_LD<-data_all%>%filter(SuSiEx_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(SuSiEx), Method = "SuSiEx" )%>%select(h2,Method,Cor )
XMAP_LD<-data_all%>%filter(XMAP_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(XMAP), Method = "XMAP" )%>%select(h2,Method,Cor )
MultiSuSiE_LD<-data_all%>%filter(MultiSuSiE_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(MultiSuSiE), Method = "MultiSuSiE" )%>%select(h2,Method,Cor )
CARMAX_LD<-data_all%>%filter(CARMAX_PIP_Either>0.5)%>% left_join(all_ld_true_signal, by = c("h2","Region"))%>%mutate(Cor = abs(CARMAX), Method = "CARMAX" )%>%select(h2,Method,Cor )




Signal_ld_data<-rbind(MESuSiE_LD,SuSiE_LD,SuSiE_merged_LD,Paintor_LD,SuSiEx_LD,XMAP_LD,MultiSuSiE_LD,CARMAX_LD)


Signal_ld_data$h2<-factor(Signal_ld_data$h2)
levels(Signal_ld_data$h2)=c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
Signal_ld_data$Method<-factor(Signal_ld_data$Method)
					  
Signal_ld_data <- Signal_ld_data %>%
   mutate(Method = fct_relevel(Method,
                               "MESuSiE" ,
                               "SuSiE",
                               "SuSiE_merged",
                               "Paintor",
                               "SuSiEx",
                               "XMAP",
                               "MultiSuSiE",
                               "CARMAX"))



p_size_signal_box = ggplot(data=Signal_ld_data,aes(x = Method, y =abs(Cor),fill=Method))+geom_boxplot( aes(fill=Method))+coord_cartesian(ylim = c(0.5, 1))+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e"))+ylab("Correlation with true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
p_size_signal_box = p_size_signal_box + theme(axis.text.x = element_blank(),
                                axis.text.y = element_text( size = 18),  
                                axis.title.x = element_text( size = 18,face="bold"),
                                axis.title.y = element_text( size = 18,face="bold"),
                                strip.text.x = element_text(size = 18),
                                strip.text.y= element_text(size = 18),
                                legend.text=element_text(size=18),
                                legend.title=element_text(size=22,face="bold"),
                                plot.title = element_text(size=16,hjust = 0.5)) 

p_out<-p_size_signal_box/p_signal_distance+ plot_layout(guides = "collect") & theme(legend.position = 'bottom')	

ggsave(paste0(plot_dir,"cor_distance_signal_PIP_",External_index_name,causal_index_name,"_meta.pdf"),  p_out ,height =10,width =18,dpi=300)

save(data_all,all_ld_true_signal,Signal_ld_data,Signal_distance_data,file = paste0(res_dir,res_out_name,"_meta.RData"))


data_all%>%summarise(sum(MESuSiE_PIP_Shared>0.5),sum(MESuSiE_PIP_WB>0.5),sum(MESuSiE_PIP_BB>0.5))/200
data_all%>%summarise(sum(SuSiE_PIP_Shared>0.5),sum(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5),sum(SuSiE_PIP_BB>0.5&SuSiE_PIP_EU<0.5))/200
data_all%>%summarise(sum(Paintor_PIP_Either>0.5))/200
Signal_distance_data%>%group_by(Method)%>%summarise(mean(Distance)/1024)
Signal_ld_data%>%group_by(Method)%>%summarise(mean(abs(Cor)))
