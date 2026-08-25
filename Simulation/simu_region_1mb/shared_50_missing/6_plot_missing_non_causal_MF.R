#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
#res_out_name<-c("shared_50_baseline_updated_meta.RData")
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
library(stringr)
library(ggpubr)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility_missing.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))





Compare_Plot_Missing <- function(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) {
  
  
  load(paste0(res_dir, res_out_1))  
  all_Set_data_dataframe_1<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_1)
  set_power_summary_1<-load_data_with_pattern(set_power_summary, pattern_name_1)
  either_all_ROC_data_dataframe_1 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_1)
  shared_all_ROC_data_dataframe_1 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_1)
  ancestry_all_ROC_data_dataframe_1 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_1)
  
  FDR_Power_either_1 <- load_data_with_pattern(FDR_Power_either, pattern_name_1)
  FDR_Power_shared_1 <- load_data_with_pattern(FDR_Power_shared, pattern_name_1)
  FDR_Power_ancestry_1 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_1)
  
  # Load the second dataset
  load(paste0(res_dir, res_out_2)) 
  
  all_Set_data_dataframe_2<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_2)
  set_power_summary_2<-load_data_with_pattern(set_power_summary, pattern_name_2)
  either_all_ROC_data_dataframe_2 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_2)
  shared_all_ROC_data_dataframe_2 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_2)
  ancestry_all_ROC_data_dataframe_2 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_2)
  
  FDR_Power_either_2 <- load_data_with_pattern(FDR_Power_either, pattern_name_2)
  FDR_Power_shared_2 <- load_data_with_pattern(FDR_Power_shared, pattern_name_2)
  FDR_Power_ancestry_2 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_2)
  
  # Load the third dataset
  load(paste0(res_dir, "shared_50_baseline_updated.RData"))
  set_power_summary <- set_power_summary %>%
    filter(Method %in% c("MESuSiE","SuSiE"))
  set_power_summary_updated <- set_power_summary %>%
    bind_rows(
      set_power_summary %>%
        filter(Method == "MESuSiE") %>%   # rows 1–6 in the example
        mutate(Method = "MFD")            # change method to MFD
    )
  set_power_summary_3<-load_data_with_pattern(set_power_summary_updated, pattern_name_3)
  all_Set_data_dataframe <- all_Set_data_dataframe %>%
    filter(Method %in% c("MESuSiE","SuSiE"))
  all_Set_data_dataframe_updated <- all_Set_data_dataframe %>%
    bind_rows(
      all_Set_data_dataframe %>%
        filter(Method == "MESuSiE") %>%   # rows 1–6 in the example
        mutate(Method = "MFD")            # change method to MFD
    )
  all_Set_data_dataframe_3<-load_data_with_pattern(all_Set_data_dataframe_updated, pattern_name_3)%>%dplyr::select(-FDR)
  either_all_ROC_data_dataframe_updated <- either_all_ROC_data_dataframe %>%
    filter(Method %in% c("MESuSiE", "SuSiE")) %>%
        bind_rows(
      filter(., Method == "MESuSiE") %>% mutate(Method = "MFD")
    )
  either_all_ROC_data_dataframe_3 <- load_data_with_pattern(either_all_ROC_data_dataframe_updated, pattern_name_3)
  shared_all_ROC_data_dataframe_updated <- shared_all_ROC_data_dataframe %>%
    filter(Method %in% c("MESuSiE", "SuSiE")) %>%
    bind_rows(
      filter(., Method == "MESuSiE") %>% mutate(Method = "MFD")
    )
  shared_all_ROC_data_dataframe_3 <- load_data_with_pattern(shared_all_ROC_data_dataframe_updated, pattern_name_3)
  ancestry_all_ROC_data_dataframe_updated <- ancestry_all_ROC_data_dataframe %>%
    filter(Method %in% c("MESuSiE WB", "MESuSiE BB","SuSiE WB","SuSiE BB")) %>%
    bind_rows(
      filter(., Method %in% c("MESuSiE WB", "MESuSiE BB")) %>%
        mutate(Method = gsub("MESuSiE", "MFD", Method))
    )
  ancestry_all_ROC_data_dataframe_3 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe_updated, pattern_name_3)
  
  FDR_Power_either_updated <- FDR_Power_either %>%
    filter(Method %in% c("MESuSiE", "SuSiE")) %>%
    bind_rows(
      filter(., Method == "MESuSiE") %>% mutate(Method = "MFD")
    )
  FDR_Power_either_3 <- load_data_with_pattern(FDR_Power_either_updated, pattern_name_3)
  FDR_Power_shared_updated <- FDR_Power_shared %>%
    filter(Method %in% c("MESuSiE", "SuSiE")) %>%
    bind_rows(
      filter(., Method == "MESuSiE") %>% mutate(Method = "MFD")
    )
  FDR_Power_shared_3 <- load_data_with_pattern(FDR_Power_shared_updated, pattern_name_3)
  FDR_Power_ancestry_updated <- FDR_Power_ancestry %>%
    filter(Method %in% c("MESuSiE WB", "MESuSiE BB", "SuSiE WB", "SuSiE BB")) %>%
    bind_rows(
      filter(., Method %in% c("MESuSiE WB", "MESuSiE BB")) %>%
        mutate(Method = gsub("MESuSiE", "MFD", Method))
    )
  FDR_Power_ancestry_3 <- load_data_with_pattern(FDR_Power_ancestry_updated, pattern_name_3)
  
  # Combine Data Together
  Set_Data <- rbind(all_Set_data_dataframe_1,all_Set_data_dataframe_2,all_Set_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  Set_Power_Data<-rbind(set_power_summary_1,set_power_summary_2,set_power_summary_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_either <- rbind(either_all_ROC_data_dataframe_1, either_all_ROC_data_dataframe_2, either_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_shared <- rbind(shared_all_ROC_data_dataframe_1, shared_all_ROC_data_dataframe_2, shared_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_ancestry <- rbind(ancestry_all_ROC_data_dataframe_1, ancestry_all_ROC_data_dataframe_2, ancestry_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  
  FDR_Power_either <- rbind(FDR_Power_either_1, FDR_Power_either_2, FDR_Power_either_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  FDR_Power_shared <- rbind(FDR_Power_shared_1, FDR_Power_shared_2, FDR_Power_shared_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  FDR_Power_ancestry <- rbind(FDR_Power_ancestry_1, FDR_Power_ancestry_2, FDR_Power_ancestry_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  
  # Define a directory to save the plots
  plot_dir <- paste0(res_dir, "Figure/",plot_dir_name,"/")
  system(paste0("mkdir -p ",plot_dir))
  # Define patterns for plotting
  #pattern_input_Power <- c("stripe","weave", "none")
  pattern_input_Power <- c("stripe","circle", "none")
  
  names(pattern_input_Power) <- c(pattern_name_1, pattern_name_2, pattern_name_3)
  
  pattern_input_line <- c("dashed","dotted", "solid")
  names(pattern_input_line) <- c(pattern_name_1, pattern_name_2, pattern_name_3)
  
  #####################################################################################
  Set_Data$Pattern<-factor(Set_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
  p_size_box<-Set_Size_Compare_fun(Set_Data%>%mutate(Size = log2(Size+1)),upper_limit = log2(max(Set_Data$Size))+1,pattern_input_Power)
  p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")
  
  Set_Power_Data$Pattern<-factor(Set_Power_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
  p_power_bar<-Set_Power_Compare_fun(Set_Power_Data,pattern_input_Power)
  p_size_power <- (p_size_box / p_power_bar) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_size_power<-p_size_power+ plot_layout(guides = "collect",heights=c(2,1))&theme(legend.position = 'bottom')+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  ggsave(paste0(plot_dir,"Set_Size_Power_Compare.pdf"),p_size_power,height=180, width=210, units = "mm",dpi=600)
  
  # Function to generate FDR power comparison plot (defined elsewhere)
  p_FDR_Power_either <- FDR_Power_Compare_fun_supp(FDR_Power_either, pattern_input_Power)
  
  # Function to generate ROC ancestry comparison plot (defined elsewhere)
  p_ROC_either <- ROC_ancestry_Compare_fun(ROC_either%>%mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")), pattern_input_line)
  
  # Create and save the combined plot for either PIP
  p_ROC_FDR_Power_either <- (p_ROC_either + p_FDR_Power_either) + 
    plot_annotation(tag_levels = 'a') & 
    theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_either <- p_ROC_FDR_Power_either + plot_layout(heights = c(2, 1))
  
  ggsave(paste0(plot_dir, "PIP_ROC_FDR_Power_Either_Compare.pdf"), p_ROC_FDR_Power_either, height = 180, width = 210, units = "mm", dpi = 600)
  
  ###Plot for Shared PIP
  
  p_FDR_Power_shared<-FDR_Power_Compare_fun_supp(FDR_Power_shared,pattern_input_Power)
  p_ROC_shared<-ROC_ancestry_Compare_fun(ROC_shared%>%mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")),pattern_input_line )
  
  p_ROC_FDR_Power_shared<-(p_ROC_shared + p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_shared<-p_ROC_FDR_Power_shared+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared_Compare.pdf"),p_ROC_FDR_Power_shared,height=180, width=210, units = "mm",dpi=600)
  
  ###Plot for WB specific
  p_FDR_Power_WB<-FDR_Power_Compare_fun_supp(FDR_Power_ancestry[grep("WB",FDR_Power_ancestry$Method),] %>%
                                          mutate(Method = case_when(
                                            str_detect(Method, "^MESuSiE(\\b|_)") ~ "MESuSiE",
                                            str_detect(Method, "^MultiSuSiE(\\b|_)") ~ "MultiSuSiE",
                                            str_detect(Method, "^SuSiEx(\\b|_)") ~ "SuSiEx",
                                            str_detect(Method, "^SuSiE_merged(\\b|_)")  ~ "SuSiE_merged",  # Match merged method before base method
                                            str_detect(Method, "^SuSiE(\\b|_)") ~ "SuSiE",
                                            str_detect(Method, "^Paintor(\\b|_)") ~ "Paintor",
                                            str_detect(Method, "^XMAP(\\b|_)") ~ "XMAP",
                                            str_detect(Method, "^CARMAX(\\b|_)") ~ "CARMAX",
                                            str_detect(Method, "^MFD(\\b|_)") ~ "MFD",
                                            TRUE ~ Method
                                          ))%>%  ungroup() %>%
                                          mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")),pattern_input_Power)
  
  p_ROC_WB<-ROC_ancestry_Compare_fun(ROC_ancestry[grep("WB",ROC_ancestry$Method),]%>%
                                       mutate(Method = case_when(
                                         str_detect(Method, "^MESuSiE(\\b|_)") ~ "MESuSiE",
                                         str_detect(Method, "^MultiSuSiE(\\b|_)") ~ "MultiSuSiE",
                                         str_detect(Method, "^SuSiEx(\\b|_)") ~ "SuSiEx",
                                         str_detect(Method, "^SuSiE_merged(\\b|_)")  ~ "SuSiE_merged",  # Match merged method before base method
                                         str_detect(Method, "^SuSiE(\\b|_)") ~ "SuSiE",
                                         str_detect(Method, "^Paintor(\\b|_)") ~ "Paintor",
                                         str_detect(Method, "^XMAP(\\b|_)") ~ "XMAP",
                                         str_detect(Method, "^CARMAX(\\b|_)") ~ "CARMAX",
                                         str_detect(Method, "^MFD(\\b|_)") ~ "MFD",
                                         TRUE ~ Method
                                       ))%>% ungroup() %>%
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")),pattern_input_line )
  
  p_ROC_FDR_Power_WB<-(p_ROC_WB + p_FDR_Power_WB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_WB<-p_ROC_FDR_Power_WB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_WB_Compare.pdf"),p_ROC_FDR_Power_WB,height=180, width=210, units = "mm",dpi=600)
  
  
  ###Plot for BB specific
  p_FDR_Power_BB<-FDR_Power_Compare_fun_supp(FDR_Power_ancestry[grep("BB",FDR_Power_ancestry$Method),] %>%
                                          mutate(Method = case_when(
                                            str_detect(Method, "^MESuSiE(\\b|_)") ~ "MESuSiE",
                                            str_detect(Method, "^MultiSuSiE(\\b|_)") ~ "MultiSuSiE",
                                            str_detect(Method, "^SuSiEx(\\b|_)") ~ "SuSiEx",
                                            str_detect(Method, "^SuSiE_merged(\\b|_)")  ~ "SuSiE_merged",  # Match merged method before base method
                                            str_detect(Method, "^SuSiE(\\b|_)") ~ "SuSiE",
                                            str_detect(Method, "^Paintor(\\b|_)") ~ "Paintor",
                                            str_detect(Method, "^XMAP(\\b|_)") ~ "XMAP",
                                            str_detect(Method, "^CARMAX(\\b|_)") ~ "CARMAX",
                                            str_detect(Method, "^MFD(\\b|_)") ~ "MFD",
                                            TRUE ~ Method
                                          ))%>% ungroup() %>%
                                          mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")),pattern_input_Power)
  
  p_ROC_BB<-ROC_ancestry_Compare_fun(ROC_ancestry[grep("BB",ROC_ancestry$Method),]%>%
                                       mutate(Method = case_when(
                                         str_detect(Method, "^MESuSiE(\\b|_)") ~ "MESuSiE",
                                         str_detect(Method, "^MultiSuSiE(\\b|_)") ~ "MultiSuSiE",
                                         str_detect(Method, "^SuSiEx(\\b|_)") ~ "SuSiEx",
                                         str_detect(Method, "^SuSiE_merged(\\b|_)")  ~ "SuSiE_merged",  # Match merged method before base method
                                         str_detect(Method, "^SuSiE(\\b|_)") ~ "SuSiE",
                                         str_detect(Method, "^Paintor(\\b|_)") ~ "Paintor",
                                         str_detect(Method, "^XMAP(\\b|_)") ~ "XMAP",
                                         str_detect(Method, "^CARMAX(\\b|_)") ~ "CARMAX",
                                         str_detect(Method, "^MFD(\\b|_)") ~ "MFD",
                                         TRUE ~ Method
                                       ))%>% ungroup() %>%
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MESuSiE","SuSiE","MFD")),pattern_input_line )
  
  p_ROC_FDR_Power_BB<-(p_ROC_BB + p_FDR_Power_BB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_BB<-p_ROC_FDR_Power_BB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_BB_Compare.pdf"),p_ROC_FDR_Power_BB,height=180, width=210, units = "mm",dpi=600)
  
}



res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_MFD"
res_out_1<-"Missing_non_causal_One_MFD.RData"
res_out_2<-"Missing_non_causal_Both_MFD.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 


res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_Ref_MFD"
res_out_1<-"Missing_non_causal_External_One_MFD.RData"
res_out_2<-"Missing_non_causal_External_Both_MFD.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 








##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
for (External_index in 1:2) {
  External_index_name <- c("", "_External")[External_index]
  
  plot_dir <- paste0(
    res_dir,
    "Figure/Missing_Causal", External_index_name, "_MFD/"
  )
  system(paste0("mkdir -p ", plot_dir))
  
  for (causal_index in 2) {
    causal_index_name <- c("Both", "One")[causal_index]
    
    rdata_file <- paste0(
      res_dir,
      "Missing_causal_",
      ifelse(External_index_name == "_External", "External_", ""),
      causal_index_name,
      "_MFD.RData"
    )
    
    # load safely into fresh env each iteration
    env_tmp <- new.env(parent = emptyenv())
    load(rdata_file, envir = env_tmp)
    
    Signal_ld_data <- env_tmp$Signal_ld_data
    Signal_distance_data <- env_tmp$Signal_distance_data
    
    p_LD_signal_box <- Signal_LD_Fun(Signal_ld_data)
    p_distance_signal_box <- Signal_Distance_Fun(Signal_distance_data)
    
    p_out <- p_LD_signal_box / p_distance_signal_box +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    ggsave(
      filename = paste0(plot_dir, "Missing_Causal_", causal_index_name, ".pdf"),
      plot = p_out,
      dpi = 600, units = "mm", height = 180, width = 210
    )
  }
}
