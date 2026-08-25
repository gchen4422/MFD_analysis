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

res_out_name<-paste0("Missing_causal_",External_index_name,causal_index_name)
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
    #set2 = which(paste0("Missing_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".mcmc.paintor")%in%list.files())
    #set3 = which(paste0("MESuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,"_additional.RData")%in%list.files())
    #set4 = which(paste0("Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",1:100,"_h2_",h2,".RData")%in%list.files())
    
    
    #replicate_set = Reduce(intersect, list(set1,set2,set3,set4))
    replicate_set = Reduce(intersect, list(set1))
    mf_decision = fread(paste0(result_dir, "Finemapping_Method_Decisions_h2_",h2, ".csv"))
    
    ld_true_signal<-c()
    data_all_h2<-c()
    
    for(LOCI_num in replicate_set){
      
      SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      load(SuSiE_res_name)
      #other_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,"_additional.RData")
      #load(other_res_name)
      #SuSiE_res_name_metal<-paste0(result_dir,"Metal_SuSiE_CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".RData")
      #load(SuSiE_res_name_metal)
      
      
      is_mesusie <- mf_decision$Method[LOCI_num] == "MESuSiE"
      
      zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>%
        mutate(
          # 1. Define SuSiE_PIP_WB first so it is available for use below
          SuSiE_PIP_WB = SuSiE_PIP_EU,
          
          # 2. Conditionally assign the MFD columns line-by-line
          # Use standard 'if' because 'is_mesusie' is a single value, not a vector.
          MFD_PIP_Either = if (is_mesusie) MESuSiE_PIP_Either else SuSiE_PIP_Either,
          MFD_PIP_WB     = if (is_mesusie) MESuSiE_PIP_WB     else SuSiE_PIP_WB,
          MFD_PIP_BB     = if (is_mesusie) MESuSiE_PIP_BB     else SuSiE_PIP_BB,
          MFD_PIP_Shared = if (is_mesusie) MESuSiE_PIP_Shared else SuSiE_PIP_Shared
          
        ) 
      
      EU_cov<-as.matrix(fread(paste0(data_dir_ld,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD1")))
      BB_cov<-as.matrix(fread(paste0(data_dir_ld,"CAUSAL_",causal_num,"_LOCI_",LOCI_num,"_h2_",h2,".LD2")))
      
      
      zfile_MESuSiE_SuSiE_Paintor$Region = LOCI_num
      data_all_h2<-rbind(data_all_h2,zfile_MESuSiE_SuSiE_Paintor)
      
      which(zfile_MESuSiE_SuSiE_Paintor$Signal!=0)
      true_signal_pos<-which(zfile_MESuSiE_SuSiE_Paintor$Signal!=0)
      
      SuSiE_signal_pos_either<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(pmax(SuSiE_PIP_BB,SuSiE_PIP_EU)))%>%pull(max_index)
      #SuSiE_signal_pos_share<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(pmin(SuSiE_PIP_BB,SuSiE_PIP_EU)))%>%pull(max_index)
      
      MESuSiE_signal_pos<-zfile_MESuSiE_SuSiE_Paintor%>%summarise(max_index = which.max(MESuSiE_PIP_Shared))%>%pull(max_index)

      
      MFD_signal_pos = ifelse(is_mesusie,MESuSiE_signal_pos,SuSiE_signal_pos_either)
      
      #EU_cov<-as.matrix(fread(paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/LD_data/Region_",gene_num,".LD1")))	
      #BB_cov<-as.matrix(fread(paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/LD_data/Region_",gene_num,".LD2")))
      max_cov<-pmax(abs(EU_cov),abs(BB_cov))
      ld_vec<-c(max_cov[true_signal_pos,c(MESuSiE_signal_pos,SuSiE_signal_pos_either,MFD_signal_pos)])
      ld_true_signal<-rbind(ld_true_signal,ld_vec)
      cat(LOCI_num)
      
    }   
    
    ld_true_signal<-data.frame(ld_true_signal)
    colnames(ld_true_signal)<-c("MESuSiE","SuSiE_Either","MFD")
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
data_all%>%summarise(sum(MFD_PIP_Shared >0.5),sum(MFD_PIP_WB>0.5&MFD_PIP_BB<0.5),sum(MFD_PIP_BB>0.5&MFD_PIP_WB<0.5))/200

#data_all%>%summarise(sum(Paintor_PIP_Either>0.5))/200
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
MFD_distance <- data_all %>% 
  filter(MFD_PIP_Either > 0.5) %>% 
  left_join(Signal_POS, by = c("h2", "Region")) %>% 
  mutate(distance_to_signal = abs(signal_pos - POS)) %>% 
  group_by(h2) %>% 
  summarise(Distance = mean(distance_to_signal)) %>% 
  mutate(Method = "MFD")




Signal_distance_data<-rbind(MESuSiE_distance,SuSiE_distance,MFD_distance)

Signal_distance_data$h2<-factor(Signal_distance_data$h2)
levels(Signal_distance_data$h2)=c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
Signal_distance_data$Method<-factor(Signal_distance_data$Method)

Signal_distance_data <- Signal_distance_data %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE" ,
                              "SuSiE",
                              "MFD"))

Signal_distance_data%>%group_by(Method)%>%summarise(mean(Distance)/1024)	

p_signal_distance = ggplot(data=Signal_distance_data,aes(x = Method, y = Distance,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e","MFD" = "#377eb8"))+ylab("Distance to true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
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
MFD_LD <- data_all %>% 
  filter(MFD_PIP_Either > 0.5) %>% 
  left_join(all_ld_true_signal, by = c("h2", "Region")) %>% 
  mutate(Cor = abs(MFD), Method = "MFD") %>% # Use the 'MFD' column from all_ld_true_signal
  select(h2, Method, Cor)


Signal_ld_data<-rbind(MESuSiE_LD,SuSiE_LD,MFD_LD)


Signal_ld_data$h2<-factor(Signal_ld_data$h2)
levels(Signal_ld_data$h2)=c(paste0("~h^2 == 10^-4"),paste0("~h^2 == 2%*%10^-4"))
Signal_ld_data$Method<-factor(Signal_ld_data$Method)

Signal_ld_data <- Signal_ld_data %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE" ,
                              "SuSiE",
                              "MFD"))



p_size_signal_box = ggplot(data=Signal_ld_data,aes(x = Method, y =abs(Cor),fill=Method))+geom_boxplot( aes(fill=Method))+coord_cartesian(ylim = c(0.5, 1))+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e","MFD" = "#377eb8"))+ylab("Correlation with true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
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

save(data_all,all_ld_true_signal,Signal_ld_data,Signal_distance_data,file = paste0(res_dir,res_out_name,"_MFD.RData"))


data_all%>%summarise(sum(MESuSiE_PIP_Shared>0.5),sum(MESuSiE_PIP_WB>0.5),sum(MESuSiE_PIP_BB>0.5))/200
data_all%>%summarise(sum(SuSiE_PIP_Shared>0.5),sum(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5),sum(SuSiE_PIP_BB>0.5&SuSiE_PIP_EU<0.5))/200
data_all%>%summarise(sum(MFD_PIP_Shared >0.5),sum(MFD_PIP_WB>0.5&MFD_PIP_BB<0.5),sum(MFD_PIP_BB>0.5&MFD_PIP_WB<0.5))/200
Signal_distance_data%>%group_by(Method)%>%summarise(mean(Distance)/1024)
Signal_ld_data%>%group_by(Method)%>%summarise(mean(abs(Cor)))
