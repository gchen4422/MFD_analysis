#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
res_out_name<-c("shared_50_baseline_updated_meta_correction_updated_Dec_17.RData")
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
library(ggpubr)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/metal/utility.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))
#for(res_out in res_out_name){
  
  res_out =res_out_name
  
  plot_dir<-paste0(res_dir,"Figure/",unlist(strsplit(res_out,".RData")),"/")
  system(paste0("mkdir -p ",plot_dir))
  ########################################################################
  #
  #           50% causal variants are shared
  #
  ########################################################################
  ###################
  #
  #Set Size&Power
  #
  ###################
  load(paste0(res_dir,res_out))   
  upper_limit<-round(all_Set_data_dataframe%>%filter(Method == "SuSiE_weighted",causal_num=="Num~Causal  == 5 ",h2=="~h^2 == 10^-4")%>%summarise(upper = quantile(Size,0.75))%>%pull(upper))+50
  p_size_box<-Set_Size_fun_order(all_Set_data_dataframe%>%mutate(Size = log2(Size+1)),upper_limit = log2(upper_limit))
  p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")+
    # Remove X-axis labels/ticks and suppress legend completely
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) 

  
  p_power_bar<-Set_Power_fun(set_power_summary)+
    # Remove X-axis labels/ticks and suppress legend completely
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) 

  
  #p_power_bar<-p_power_bar+ylim(0,min(max(set_power_summary$Power_name)+0.05,1))
  size_power<-p_size_box/p_power_bar+plot_annotation(tag_levels = 'a')& 
    theme(plot.tag = element_text(size = 7,face="bold"))
  
  ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=180, width=210, units = "mm",dpi=600)
  
  
  ##########################################################
  #
  # Either ancestry 
  # Set Size&Power | ROC | FDR Power 
  #
  ##########################################################
  
  ###################
  #
  #ROC
  #
  ###################
  either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor"))
  p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
  ###################
  #
  #FDR&Power
  #
  ###################
  power_upper_limit<-FDR_Power_either%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
  p_FDR_Power_Either<-FDR_Power_shared_fun(FDR_Power_either%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
  #ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  

  
  ##########################################################
  #
  # Shared Signal 
  #  ROC | FDR Power | PIP calibration | PIP scatter plot 
  #
  ##########################################################
  ###################
  #
  #ROC
  #
  ###################
  shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"))
  p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
  ###################
  #
  #FDR&Power
  #
  ###################
  
  power_upper_limit<-FDR_Power_shared%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
  
  #  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared)+ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1))
  #  ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  #  ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
  #  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=180, width=180, units = "mm",dpi=500)
  

  
  
  # p_calibration_shared<-PIP_calibration_shared_fun(PIP_calibration_Shared)
  # ggsave(paste0(plot_dir,"PIP_Calibration_Shared.pdf"),p_calibration_shared,height = 121, width = 121,units = "mm",dpi=500)
  
#}
