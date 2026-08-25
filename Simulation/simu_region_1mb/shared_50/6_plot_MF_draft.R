#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")

########################################################################
#
#           BASELINE simulation Figure
#
########################################################################

res_out_name<-c("shared_50_baseline_updated_meta_updated_xmap.RData")
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
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))

custom_theme <- function() {
  theme_bw() +  # Start with the standard black & white theme
    theme(
      # 1. Remove all grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Configure text sizes and styles
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = element_text(size = 9, face = "bold"),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 9, hjust = 0.5)
    )
}


#for(res_out in res_out_name){
  
  res_out = res_out_name
  plot_dir<-paste0(res_dir,"Figure/Draft_main_figure/")
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
  method_levels_base <- c("MESuSiE", "SuSiE", "SuSiE_weighted","SuSiE_merged", 
                         "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
  
  load(paste0(res_dir,res_out))   
  all_Set_data_dataframe <- all_Set_data_dataframe %>%
    mutate(Method = recode(Method,
                           "SuSiE_meta_weighted" = "SuSiE-weighted",
                           "SuSiE_meta_merged"   = "SuSiE-merged"))
  upper_limit<-round(all_Set_data_dataframe%>%filter(Method == "Paintor",causal_num=="Num~Causal  == 5 ",h2=="~h^2 == 10^-4")%>%summarise(upper = quantile(Size,0.75))%>%pull(upper))+50
  p_size_box<-Set_Size_fun_order(all_Set_data_dataframe%>%mutate(Size = log2(Size+1)),upper_limit = log2(upper_limit))
  p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")
  
  set_power_summary <- set_power_summary %>%
    mutate(Method = recode(Method,
                           "SuSiE_meta_weighted" = "SuSiE_weighted",
                           "SuSiE_meta_merged"   = "SuSiE_merged")) %>% mutate(Method  = factor(Method, levels = method_levels_base))


  p_power_bar<-Set_Power_fun(set_power_summary)

  all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))
  
  ####################
  #
  #PIP calibration
  #
  ####################
  # ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  #  ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
  
  PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
  PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,method_levels_base))
  
  p_calibration_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either_byh2)
  

  
  
  ##############################################################################################################################################################
  
  
  
  
  ##########################################################
  #
  #
  #             Summarized result
  #
  #
  ##########################################################
  

  # 1. PREPARE DATA
  # Aggregate 'Either' data
  either_mean <- either_all_ROC_data_dataframe %>%
    group_by(Method, Cutoff) %>%
    summarise(Power = mean(Power, na.rm = TRUE), 
              FDR = mean(FDR, na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(Signal = "Either") # Add label
  
  # Aggregate 'Shared' data
  shared_mean <- shared_all_ROC_data_dataframe %>%
    group_by(Method, Cutoff) %>%
    summarise(Power = mean(Power, na.rm = TRUE), 
              FDR = mean(FDR, na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(Signal = "Shared") # Add label
  
  # Combine into one dataframe
  combined_roc_data <- bind_rows(either_mean, shared_mean)
  
  # Set factor levels to control facet order (Either left, Shared right)
  combined_roc_data$Signal <- factor(combined_roc_data$Signal, levels = c("Either", "Shared"))
  
  
  # 2. PLOT
  p_roc_combined <- ggplot(combined_roc_data, aes(x = Power, y = 1 - FDR)) +
    
    # Main Lines
    geom_line(aes(color = Method), size = 0.5) +
    
    # PIP > 0.5 Points
    geom_point(data = combined_roc_data %>% filter(Cutoff == 0.5),
               aes(shape = 'PIP > 0.5', color = Method),
               size = 2.5, stroke = 0.8, fill = "white") +
    
    # Faceting: Creates the two panels side-by-side
    facet_wrap(~ Signal, scales = "fixed") +
    
    # Colors
    scale_color_manual(values = c(
      "MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
      "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
      "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMAX"="#f2b56e"
    )) +
    
    # Shape for PIP point
    scale_shape_manual(name = " ", values = c('PIP > 0.5' = 21)) +
    
    # Scales & Labels
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    xlab("Recall") +
    ylab("Precision") +
    
    # Theme
    custom_theme() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = NA), # Optional: distinctive strip header
      strip.text = element_text(size = 8, face = "bold")
    ) +
    
    # Legend Guides
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    )
  
  # 3. RENDER
  print(p_roc_combined)
  
  
  
  FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
  split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
  FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
    Method = sapply(split_list, `[`, 1),
    Ancestry = sapply(split_list, `[`, 2)
  )%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))
  
  
  FDR_Power_Combined<-rbind( FDR_Power_either%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                               filter(FDR==0.05)%>%mutate(Cat = "Either"),
                             FDR_Power_shared%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%
                               filter(FDR==0.05)%>%mutate(Cat = "Shared"),
                             FDR_Power_ancestry%>%group_by(FDR,Method)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.05)%>%mutate(Cat = "Ancestry"))
  
  
  
  FDR_Power_Combined$Method <- factor(
    FDR_Power_Combined$Method,
    levels = method_levels_base
  )
  FDR_Power_Combined$Cat<-factor(	 FDR_Power_Combined$Cat,levels = c("Either","Shared","Ancestry"))
  #FDR_Power_Combined$Method<-factor(FDR_Power_Combined$Method,levels = c("SuSiE","Paintor","Paintor_fun","MESuSiE","MESuSiE_fun"))
  p_FDR_Power <- ggplot(FDR_Power_Combined, aes(Cat, Power, fill = Method)) +
    # 1. Set explicit dodge width for bars (standard is 0.9)
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    
    # 2. Match that width exactly in geom_text
    geom_text(
      aes(x = Cat, group = Method, y = Power, label = Power),
      position = position_dodge(width = 0.9), # MUST match the bar width
      vjust = -0.5,
      size = 5 * 5 / 10
    ) +
    
    scale_fill_manual(values = c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5",
                                 "SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95",
                                 "Paintor"="#fc8d62","MultiSuSiE"="#e78ac3",
                                 "SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e")) +
    ylab("Power") + 
    xlab("FDR = 0.05") + 
    ylim(0, 0.3) +
    theme_bw() + 
    custom_theme() + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2) +
    theme(legend.position = 'bottom')

  # Define the method order
  method_order <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
  
  # 2. Force p_roc_combined to use this specific order
  # Re-add scale_color_manual with the 'breaks' argument.
  p_roc_combined <- p_roc_combined + 
    scale_color_manual(
      values = c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                 "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                 "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMAX"="#f2b56e"),
      breaks = method_order # Preserve legend order
    ) +
    theme(legend.position = "bottom")
  
  
  # Create placeholder plots for the target layout
  p_a <- p_size_box +  custom_theme()
  p_b <- p_power_bar +  custom_theme()
  p_c <- p_roc_combined+ custom_theme()
  p_d <- p_calibration_byh2 + custom_theme()
  p_e <- p_FDR_Power + custom_theme()
  
  p_a <- p_a + 
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()  # Removes the little tick marks
    ) + 
    guides(fill = "none", color = "none")
  
  p_b <- p_b + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()  # Removes the little tick marks
    ) + 
    guides(fill = "none", color = "none")
  
  p_d <- p_d + guides(fill = "none", color = "none")
  p_e <- p_e + guides(fill = "none", color = "none")
  

  p_roc_combined <- p_roc_combined + theme(legend.position = "bottom")
  
  top <- p_a + p_b + p_roc_combined + plot_layout(widths = c(1, 1, 1.3))
  
  # Bottom Ratio: 2 vs 1.2   ->  Matches Top!
  bottom <- p_d + p_e + plot_layout(widths = c(2, 1.2))
  
  # 4. Combine
  final <- top / bottom + 
    plot_annotation(
      tag_levels = 'a',
      theme = theme(plot.tag = element_text(face = "bold", size = 12)) # Scaled tags up slightly to match
    ) +
    plot_layout(guides = "collect") & 
    # Apply theme settings to ALL subplots
    theme(
      legend.position = "bottom"
    ) &
    # Force single-row legend
    guides(
      color = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(nrow = 1, order = 2)
    )
  
  print(final)
  
  #print(final)
  
  ggsave(paste0(plot_dir,"Baseline_simulation_summary_updated_test.pdf"),final,height=300, width=600, units = "mm", dpi = 600)
  
#}      



########################################################################
#
#           Missing SNP simulation Figure
#
########################################################################






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

dyn.load("/apps/spack/negishi/apps/curl/7.85.0-gcc-12.2.0-4phaxqw/lib/libcurl.so.4",
         local = FALSE)

dyn.load("/apps/spack/negishi/apps/r-sf/1.0-9-gcc-12.2.0-sftozdw/rlib/R/library/sf/libs/sf.so")




#####################################################################################################################
#
#
#                     Missing Non-causal Part
#
#
#####################################################################################################################

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Draft_main_figure"
res_out_1<-"Missing_non_causal_One_meta_updated.RData"
res_out_2<-"Missing_non_causal_Both_meta_updated.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 



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
  load(paste0(res_dir, "shared_50_baseline_updated_meta_updated_xmap.RData"))
  set_power_summary <- set_power_summary %>%
    filter(Method != "SuSiE_meta_weighted") %>%
    mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))
  set_power_summary_3<-load_data_with_pattern(set_power_summary, pattern_name_3)
  all_Set_data_dataframe <- all_Set_data_dataframe %>%
    filter(Method != "SuSiE_meta_weighted") %>%
    mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))
  all_Set_data_dataframe_3<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_3)%>%dplyr::select(-FDR)
  either_all_ROC_data_dataframe <- either_all_ROC_data_dataframe %>%
    filter(Method != "SuSiE_weighted")
  either_all_ROC_data_dataframe_3 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_3)
  shared_all_ROC_data_dataframe <- shared_all_ROC_data_dataframe %>%
    filter(Method != "SuSiE_weighted")
  shared_all_ROC_data_dataframe_3 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_3)
  ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%
    filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))
  ancestry_all_ROC_data_dataframe_3 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_3)
  
  FDR_Power_either <- FDR_Power_either %>%
    dplyr::filter(Method != "SuSiE_weighted")
  FDR_Power_either_3 <- load_data_with_pattern(FDR_Power_either, pattern_name_3)
  FDR_Power_shared <- FDR_Power_shared %>%
    dplyr::filter(Method != "SuSiE_weighted")
  FDR_Power_shared_3 <- load_data_with_pattern(FDR_Power_shared, pattern_name_3)
  FDR_Power_ancestry <- FDR_Power_ancestry %>%
    filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))
  FDR_Power_ancestry_3 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_3)
  
  # Combine Data Together
  Set_Data <- rbind(all_Set_data_dataframe_1,all_Set_data_dataframe_2,all_Set_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  Set_Power_Data<-rbind(set_power_summary_1,set_power_summary_2,set_power_summary_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_either <- rbind(either_all_ROC_data_dataframe_1, either_all_ROC_data_dataframe_2, either_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_shared <- rbind(shared_all_ROC_data_dataframe_1, shared_all_ROC_data_dataframe_2, shared_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  ROC_ancestry <- rbind(ancestry_all_ROC_data_dataframe_1, ancestry_all_ROC_data_dataframe_2, ancestry_all_ROC_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  
  FDR_Power_either <- rbind(FDR_Power_either_1, FDR_Power_either_2, FDR_Power_either_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  FDR_Power_shared <- rbind(FDR_Power_shared_1, FDR_Power_shared_2, FDR_Power_shared_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  FDR_Power_ancestry <- rbind(FDR_Power_ancestry_1, FDR_Power_ancestry_2, FDR_Power_ancestry_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  
  
  Set_Data%>%group_by(Method,Pattern)%>%summarise(round(median(Size)),round(mean(Power),2))
  
  
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
  #p_size_power <- (p_size_box / p_power_bar) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  #p_size_power<-p_size_power+ plot_layout(guides = "collect",heights=c(2,1))&theme(legend.position = 'bottom')
  
  
  # Function to generate FDR power comparison plot (defined elsewhere)
  FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
  split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
  FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
    Method = sapply(split_list, `[`, 1),
    Ancestry = sapply(split_list, `[`, 2)
  )%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))
  
  FDR_Power_Combined<-rbind( FDR_Power_either%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%
                               filter(FDR==0.05)%>%mutate(Cat = "Either"),
                             FDR_Power_shared%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%
                               filter(FDR==0.05)%>%mutate(Cat = "Shared"),
                             FDR_Power_ancestry%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.05)%>%mutate(Cat = "Ancestry"))
  
  
  
  FDR_Power_Combined$Method <- factor(
    FDR_Power_Combined$Method,
    levels = c(
      "MESuSiE", "SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"
    )
  )
  
  FDR_Power_Combined$Cat<-factor(	 FDR_Power_Combined$Cat,levels = c("Either","Shared","Ancestry"))

  p_FDR_Power<-FDR_Power_Compare_fun(FDR_Power_Combined,pattern_input_Power)
  






res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Draft_main_figure"
res_out_1<-"Missing_non_causal_External_One_meta_updated.RData"
res_out_2<-"Missing_non_causal_External_Both_meta_updated.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 








##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
#for (External_index in 1:2) {

External_index = 1
External_index_name <- c("", "_External")[External_index]

# Create a list to store the column plots for each iteration
plot_list <- list()

for (causal_index in 1:2) {
  causal_index_name <- c("Both", "One")[causal_index]
  
  rdata_file <- paste0(
    res_dir,
    "Missing_causal_",
    ifelse(External_index_name == "_External", "External_", ""), 
    causal_index_name,
    "_updated_meta.RData"
  )
  
  # Load safely into fresh env each iteration
  env_tmp <- new.env(parent = emptyenv())
  load(rdata_file, envir = env_tmp)
  
  Signal_ld_data <- env_tmp$Signal_ld_data
  Signal_distance_data <- env_tmp$Signal_distance_data
  
  # Generate the individual subplots
  p_LD_signal_box <- Signal_LD_Fun(Signal_ld_data)+ 
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()  # Removes the little tick marks
    ) 
  p_distance_signal_box <- Signal_Distance_Fun(Signal_distance_data)+ 
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()  # Removes the little tick marks
    ) 
  
  # Add a title to the top plot to distinguish the scenario
  p_LD_signal_box <- p_LD_signal_box + ggtitle(paste("Scenario: Missing Causal SNP in", causal_index_name, "Ancestry group"))
  
  # Combine vertically (LD on top of Distance) for this specific causal index
  # Do not add layout/guides here yet; add them at the final stage
  plot_list[[causal_index]] <- p_LD_signal_box / p_distance_signal_box
}

# Combine the two columns side-by-side
# index 1 is "Both", index 2 is "One"
final_figure <- (plot_list[[1]] | plot_list[[2]]) +
  plot_layout(guides = "collect") &    # Collects legends into one
  theme(legend.position = "bottom")

# View the final combined object
final_figure

####################################################################################
#
#
#			Combine all the missing SNP figures
#
#
#
####################################################################################



p_a = p_size_box+custom_theme()
p_b = p_power_bar+custom_theme()
p_c = p_FDR_Power
p_d = final_figure+custom_theme()

blank_x <- theme(
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

drop_legends <- theme(legend.position = "none")

# Global theme: DO NOT touch axis.text / axis.text.x here
base_theme <- theme(
  axis.text.y  = element_text(size = 9),
  axis.title.y = element_text(size = 9, face = "bold"),
  strip.text.x = element_text(size = 9),
  strip.text.y = element_text(size = 9)
)

p_a <- p_a + blank_x + drop_legends
p_c <- p_c + drop_legends
p_d <- p_d + drop_legends

# p_b is the ONLY legend source; force single-row legends HERE
p_b <- p_b + blank_x +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    fill  = guide_legend(nrow = 1, byrow = TRUE, order = 2),
    shape = guide_legend(nrow = 1, byrow = TRUE, order = 3)
    # Disable the ggpattern legend if present:
    # pattern = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

row_top    <- p_a + p_b
row_bottom <- p_c + p_d + plot_layout(widths = c(1, 1))

final <- row_top / row_bottom +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    )
  ) &
  base_theme

ggsave(
  paste0(plot_dir, "Missing_simulation_summary_updated.pdf"),
  final,
  height = 300, width = 500, units = "mm", dpi = 600
)




########################################################################
#
#           LD mismatch simulation Figure
#
########################################################################

## In this script /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/10_combined_plot_draft_MF.R

########################################################################
#
#           MFD simulation Figure
#
########################################################################


library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# --- 1. Define Nodes with Detailed Parameters ---
# Adjusted Y-coordinates for better vertical spacing
nodes <- tribble(
  ~id, ~type,      ~x,   ~y,    ~fill_color, ~label_text,
  1,   "box",      0,    5.5,   "#F2F9FF",   "Start: Locus Analysis\n(Input: GWAS SumStats)",
  2,   "diamond",  0,    4.2,   "#FFF5E6",   "Are AS-Vs significant?\n(P < 5e-8 in specific ancestry)",
  3,   "diamond",  0,    2.4,   "#FFF5E6",   "In High LD with Shared SNPs?\n(r2 > 0.6 with lead Shared Signal)",
  4,   "diamond",  0,    0.6,   "#FFF5E6",   "Shared SNPs significant in ALL?\n(P < 5e-8 in both/all pops)",
  5,   "result",  -3.2,  2.4,   "#FFE6E6",   "Use SuSiE post-hoc\n(Distinct Signals or\nAncestry-Specific Drivers)",
  6,   "result",   3.2,  3.0,   "#E6FCE6",   "Use MESuSiE\n(Shared Causal Signal\nor No valid AS-V)"
)

# --- 2. Helper Functions for Shapes ---

# Generate Diamond Coordinates (Widened to w=2.6 to fit text)
make_diamond <- function(cx, cy, w=2.6, h=1.0) {
  tibble(x = c(cx, cx + w/2, cx, cx - w/2),
         y = c(cy + h/2, cy, cy - h/2, cy))
}

# Generate Box Coordinates
make_box <- function(cx, cy, w=2.4, h=1.0) {
  tibble(xmin = cx - w/2, xmax = cx + w/2,
         ymin = cy - h/2, ymax = cy + h/2)
}

# --- 3. Generate Shape Data Frames ---

# Generate Diamonds
diamond_nodes <- nodes %>% filter(type == "diamond")
diamond_data <- bind_rows(lapply(1:nrow(diamond_nodes), function(i) {
  row <- diamond_nodes[i, ]
  coords <- make_diamond(row$x, row$y, w=2.6, h=1.0)
  coords$id <- row$id
  coords$fill <- row$fill_color
  return(coords)
}))

# Generate Boxes (Start and Results)
box_nodes <- nodes %>% filter(type != "diamond")
box_data <- bind_rows(lapply(1:nrow(box_nodes), function(i) {
  row <- box_nodes[i, ]
  # Make result boxes slightly taller/wider
  w_box <- ifelse(row$type == "result", 2.8, 2.4)
  h_box <- ifelse(row$type == "result", 1.2, 0.8)
  
  coords <- make_box(row$x, row$y, w=w_box, h=h_box)
  coords$id <- row$id
  coords$fill <- row$fill_color
  return(coords)
}))

# --- 4. Define Arrows with Manhattan (Right-Angle) Paths ---
# H = half height of diamond (0.5), W = half width of diamond (1.3)
arrows <- tribble(
  ~x,   ~y,   ~xend, ~yend, ~label, ~lbl_x, ~lbl_y,
  
  # 1. Start -> Q1
  0,    5.1,  0,     4.7,   NA,     NA,     NA,
  
  # 2. Q1 (No) -> MESuSiE (Right Result)
  1.3,  4.2,  3.2,   4.2,   "No",   1.6,    4.35, # Horizontal out
  3.2,  4.2,  3.2,   3.6,   NA,     NA,     NA,   # Vertical down to box top
  
  # 3. Q1 (Yes) -> Q2 (Down)
  0,    3.7,  0,     2.9,   "Yes",  0.2,    3.3,
  
  # 4. Q2 (No) -> SuSiE (Left Result)
  -1.3, 2.4,  -1.8,  2.4,   "No",   -1.5,   2.55,
  
  # 5. Q2 (Yes) -> Q3 (Down)
  0,    1.9,  0,     1.1,   "Yes",  0.2,    1.5,
  
  # 6. Q3 (No) -> SuSiE (Left connection)
  -1.3, 0.6,  -3.2,  0.6,   "No",   -1.6,   0.75, # Horizontal out left
  -3.2, 0.6,  -3.2,  1.8,   NA,     NA,     NA,   # Vertical up to box bottom
  
  # 7. Q3 (Yes) -> MESuSiE (Right connection)
  1.3,  0.6,  3.2,   0.6,   "Yes",  1.6,    0.75, # Horizontal out right
  3.2,  0.6,  3.2,   2.4,   NA,     NA,     NA    # Vertical up to box bottom
)

# --- 5. Build the Plot ---
p_flowchart <- ggplot() +
  # 1. Draw Edges (Arrows)
  geom_segment(data = arrows, aes(x=x, y=y, xend=xend, yend=yend),
               arrow = arrow(length = unit(0.25, "cm"), type="closed"), 
               color="gray30", size=0.8, lineend = "round", linejoin = "mitre") +
  
  # 2. Draw Arrow Labels (Yes/No)
  geom_text(data = arrows, aes(x=lbl_x, y=lbl_y, label=label), 
            size=3.5, fontface="italic", color="gray20") +
  
  # 3. Draw Shapes (Diamonds and Boxes)
  geom_polygon(data = diamond_data, aes(x=x, y=y, group=id, fill=fill), 
               color="black", size=0.4) +
  geom_rect(data = box_data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
            color="black", size=0.4) +
  
  # 4. Draw Node Labels (with string wrapping)
  geom_text(data = nodes, aes(x=x, y=y, label=label_text), 
            size=3, fontface="bold", lineheight=1.1, color="black") +
  
  # 5. Scales and Theme
  scale_fill_identity() +
  coord_fixed(ratio = 1, xlim = c(-4.2, 4.2), ylim = c(0, 6)) +
  theme_void() +
  theme(
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill="white", color=NA)
  )

# Print
print(p_flowchart)





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



#####################################################################################################################
#
#
#                     Missing Non-causal Part
#
#
#####################################################################################################################

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Draft_main_figure"
res_out_1<-"Missing_non_causal_One_MFD.RData"
res_out_2<-"Missing_non_causal_Both_MFD.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 



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


#####################################################################################
Set_Data$Pattern<-factor(Set_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
p_size_box<-Set_Size_Compare_fun(Set_Data%>%mutate(Size = log2(Size+1)),upper_limit = log2(max(Set_Data$Size))+1,pattern_input_Power)
p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")

Set_Power_Data$Pattern<-factor(Set_Power_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
p_power_bar<-Set_Power_Compare_fun(Set_Power_Data,pattern_input_Power)
#p_size_power <- (p_size_box / p_power_bar) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#p_size_power<-p_size_power+ plot_layout(guides = "collect",heights=c(2,1))&theme(legend.position = 'bottom')


##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
#for (External_index in 1:2) {

External_index = 1
External_index_name <- c("", "_External")[External_index]

# Create a list to store the column plots for each iteration
plot_list <- list()

causal_index = 2

causal_index_name <- c("Both", "One")[causal_index]
  
rdata_file <- paste0(
    res_dir,
    "Missing_causal_",
    ifelse(External_index_name == "_External", "External_", ""), 
    causal_index_name,
    "_MFD.RData"
)
  
  # Load safely into fresh env each iteration
env_tmp <- new.env(parent = emptyenv())
load(rdata_file, envir = env_tmp)

Signal_ld_data <- env_tmp$Signal_ld_data
Signal_distance_data <- env_tmp$Signal_distance_data

  # Generate the individual subplots
p_LD_signal_box <- Signal_LD_Fun(Signal_ld_data)+custom_theme()
  
p_distance_signal_box <- Signal_Distance_Fun(Signal_distance_data)+custom_theme()
  
  
  # Add a title to the top plot to distinguish the scenario
#p_LD_signal_box <- p_LD_signal_box + ggtitle(paste("Scenario: Missing Causal SNP in", causal_index_name, "Ancestry group"))
  
p_out <- p_LD_signal_box + p_distance_signal_box +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

####################################################################################
#
#
#			Combine all the missing SNP figures
#
#
####################################################################################



p_a = p_flowchart
p_b = p_size_box
p_c = p_power_bar
p_d = p_out+custom_theme()



# --- 1. Clean up individual plots ---
# Remove legends to prevent duplicates
p_a <- p_a + theme(legend.position = "none") 
p_b <- p_b + theme(legend.position = "none") + custom_theme() 
p_d <- p_d + theme(legend.position = "none") 

# Configure master legend on p_c
p_c <- p_c +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    fill  = guide_legend(nrow = 1, byrow = TRUE, order = 2),
    shape = guide_legend(nrow = 1, byrow = TRUE, order = 3)
  ) +
  theme(legend.position = "bottom") + custom_theme()

# --- 2. Define Layout Topology ---

# Right Column: Stack B, C, D vertically
col_right <- (p_b / p_c / p_d) + 
  plot_layout(heights = c(0.8, 0.8, 0.5)) 

# --- 3. Combine and Finalize ---

# Left Side: p_a
# Right Side: col_right
# Use 'widths' to control the ratio. c(1.8, 1) gives p_a nearly 2x the width of the right column.
final <- (p_a | col_right) +
  plot_layout(
    guides = "collect",
    widths = c(1.6, 1)  # Increased from 1.6 to 1.8 to allow p_a to grow larger
  ) +
  plot_annotation(
    tag_levels = 'a',        
    tag_prefix = '(',        
    tag_suffix = ')',       
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold")
    )
  ) & 
  theme(
    # Legend settings
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(t = 0, b = 0)
  )

# View
final

# When saving, ensure the width is wide enough to support the 1.8:1 ratio
ggsave(
  paste0(plot_dir, "MFD_summary_final_updated.pdf"),
  final,
  height = 250, width = 500, units = "mm", dpi = 600
)









########################################################################
#
#           Missing SNP simulation Figure
#
########################################################################






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

#dyn.load("/apps/spack/negishi/apps/curl/7.85.0-gcc-12.2.0-4phaxqw/lib/libcurl.so.4",
#         local = FALSE)

#dyn.load("/apps/spack/negishi/apps/r-sf/1.0-9-gcc-12.2.0-sftozdw/rlib/R/library/sf/libs/sf.so")




#####################################################################################################################
#
#
#                     Missing Non-causal Part
#
#
#####################################################################################################################

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Draft_main_figure"
res_out_1<-"Missing_non_causal_One_MFD.RData"
res_out_2<-"Missing_non_causal_Both_MFD.RData"
res_out_1_external = "Missing_non_causal_One_meta_updated.RData"
res_out_2_external = "Missing_non_causal_Both_meta_updated.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 



pattern_name_1 = "Missing in One Ancestry"

# ---------------------------------------------------------
# 1. LOAD MFD DATA (Extract only MFD rows)
# ---------------------------------------------------------
load(paste0(res_dir, res_out_1))  

# Filter and store MFD data temporarily
all_Set_data_mfd      <- load_data_with_pattern(all_Set_data_dataframe, pattern_name_1) %>% filter(Method == "MFD")
set_power_summary_mfd <- load_data_with_pattern(set_power_summary, pattern_name_1)      %>% filter(Method == "MFD")

# ROC Data
either_all_ROC_mfd   <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_1)   %>% filter(Method == "MFD")
shared_all_ROC_mfd   <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_1)   %>% filter(Method == "MFD")
ancestry_all_ROC_mfd <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_1) %>% filter(Method %in% c("MFD BB", "MFD WB"))

# FDR Power Data
FDR_Power_either_mfd   <- load_data_with_pattern(FDR_Power_either, pattern_name_1)   %>% filter(Method == "MFD")
FDR_Power_shared_mfd   <- load_data_with_pattern(FDR_Power_shared, pattern_name_1)   %>% filter(Method == "MFD")
FDR_Power_ancestry_mfd <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_1) %>% filter(Method %in% c("MFD BB", "MFD WB"))

# ---------------------------------------------------------
# 2. LOAD MAIN DATA & COMBINE
# ---------------------------------------------------------
load(paste0(res_dir, res_out_1_external))

# Combine Main Data (loaded above) with extracted MFD Data
all_Set_data_dataframe_1 <- bind_rows(
  load_data_with_pattern(all_Set_data_dataframe, pattern_name_1),
  all_Set_data_mfd
)

set_power_summary_1 <- bind_rows(
  load_data_with_pattern(set_power_summary, pattern_name_1),
  set_power_summary_mfd
)

# Combine ROC
either_all_ROC_data_dataframe_1 <- bind_rows(
  load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_1),
  either_all_ROC_mfd
)

shared_all_ROC_data_dataframe_1 <- bind_rows(
  load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_1),
  shared_all_ROC_mfd
)

ancestry_all_ROC_data_dataframe_1 <- bind_rows(
  load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_1),
  ancestry_all_ROC_mfd
)

# Combine FDR Power
FDR_Power_either_1 <- bind_rows(
  load_data_with_pattern(FDR_Power_either, pattern_name_1),
  FDR_Power_either_mfd
)

FDR_Power_shared_1 <- bind_rows(
  load_data_with_pattern(FDR_Power_shared, pattern_name_1),
  FDR_Power_shared_mfd
)

FDR_Power_ancestry_1 <- bind_rows(
  load_data_with_pattern(FDR_Power_ancestry, pattern_name_1),
  FDR_Power_ancestry_mfd
)



pattern_name_2 = "Missing in Both Ancestry"

# ---------------------------------------------------------
# 1. LOAD MFD DATA (Extract only MFD rows)
# ---------------------------------------------------------
# Load the MFD specific results
load(paste0(res_dir, res_out_2))

# Filter and temporarily store MFD data
all_Set_data_mfd_2      <- load_data_with_pattern(all_Set_data_dataframe, pattern_name_2) %>% filter(Method == "MFD")
set_power_summary_mfd_2 <- load_data_with_pattern(set_power_summary, pattern_name_2)      %>% filter(Method == "MFD")

# ROC Data
either_all_ROC_mfd_2   <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_2)   %>% filter(Method == "MFD")
shared_all_ROC_mfd_2   <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_2)   %>% filter(Method == "MFD")
ancestry_all_ROC_mfd_2 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_2) %>% filter(Method %in% c("MFD BB", "MFD WB"))

# FDR Power Data
FDR_Power_either_mfd_2   <- load_data_with_pattern(FDR_Power_either, pattern_name_2)   %>% filter(Method == "MFD")
FDR_Power_shared_mfd_2   <- load_data_with_pattern(FDR_Power_shared, pattern_name_2)   %>% filter(Method == "MFD")
FDR_Power_ancestry_mfd_2 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_2) %>% filter(Method %in% c("MFD BB", "MFD WB"))

# ---------------------------------------------------------
# 2. LOAD MAIN DATA & COMBINE
# ---------------------------------------------------------
# Load the standard "meta_updated" results
load(paste0(res_dir, res_out_2_external)) 

# Combine Main Data with extracted MFD Data
all_Set_data_dataframe_2 <- bind_rows(
  load_data_with_pattern(all_Set_data_dataframe, pattern_name_2),
  all_Set_data_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

set_power_summary_2 <- bind_rows(
  load_data_with_pattern(set_power_summary, pattern_name_2),
  set_power_summary_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

# Combine ROC
either_all_ROC_data_dataframe_2 <- bind_rows(
  load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_2),
  either_all_ROC_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

shared_all_ROC_data_dataframe_2 <- bind_rows(
  load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_2),
  shared_all_ROC_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

ancestry_all_ROC_data_dataframe_2 <- bind_rows(
  load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_2),
  ancestry_all_ROC_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

# Combine FDR Power
FDR_Power_either_2 <- bind_rows(
  load_data_with_pattern(FDR_Power_either, pattern_name_2),
  FDR_Power_either_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

FDR_Power_shared_2 <- bind_rows(
  load_data_with_pattern(FDR_Power_shared, pattern_name_2),
  FDR_Power_shared_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)

FDR_Power_ancestry_2 <- bind_rows(
  load_data_with_pattern(FDR_Power_ancestry, pattern_name_2),
  FDR_Power_ancestry_mfd_2
) %>% arrange(h2, causal_num, Pattern, Method)






# Load the third dataset

pattern_name_3 = "None Missing"



# Load the dataset
load(paste0(res_dir, "shared_50_baseline_updated_meta_updated_xmap.RData"))

# ---------------------------------------------------------
# 1. Set Size and Power Summary
# ---------------------------------------------------------

# Clean original
set_power_summary <- set_power_summary %>%
  filter(Method != "SuSiE_meta_weighted") %>%
  mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))

# Create MFD (Clone MESuSiE)
set_power_mfd <- set_power_summary %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

# Combine and Apply Pattern
set_power_summary_3 <- bind_rows(set_power_summary, set_power_mfd) %>%
  load_data_with_pattern(pattern_name_3) %>%
  arrange(h2, causal_num, Method)

# ---------------------------------------------------------
# 2. All Set Data DataFrame
# ---------------------------------------------------------

# Clean original
all_Set_data_dataframe <- all_Set_data_dataframe %>%
  filter(Method != "SuSiE_meta_weighted") %>%
  mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))

# Create MFD (Clone MESuSiE)
all_Set_data_mfd <- all_Set_data_dataframe %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

# Combine, Apply Pattern, and Remove FDR column
all_Set_data_dataframe_3 <- bind_rows(all_Set_data_dataframe, all_Set_data_mfd) %>%
  load_data_with_pattern(pattern_name_3) %>%
  dplyr::select(-FDR) %>%
  arrange(h2, causal_num, Method)

# ---------------------------------------------------------
# 3. ROC Data (Either & Shared)
# ---------------------------------------------------------

# --- Either ---
either_all_ROC_data_dataframe <- either_all_ROC_data_dataframe %>%
  filter(Method != "SuSiE_weighted")

either_mfd <- either_all_ROC_data_dataframe %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

either_all_ROC_data_dataframe_3 <- bind_rows(either_all_ROC_data_dataframe, either_mfd) %>%
  load_data_with_pattern(pattern_name_3)

# --- Shared ---
shared_all_ROC_data_dataframe <- shared_all_ROC_data_dataframe %>%
  filter(Method != "SuSiE_weighted")

shared_mfd <- shared_all_ROC_data_dataframe %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

shared_all_ROC_data_dataframe_3 <- bind_rows(shared_all_ROC_data_dataframe, shared_mfd) %>%
  load_data_with_pattern(pattern_name_3)

# ---------------------------------------------------------
# 4. ROC Data (Ancestry)
# ---------------------------------------------------------

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%
  filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))

# Create MFD (Clone MESuSiE WB/BB -> MFD WB/BB)
ancestry_mfd <- ancestry_all_ROC_data_dataframe %>%
  filter(Method %in% c("MESuSiE WB", "MESuSiE BB")) %>%
  mutate(Method = case_when(
    Method == "MESuSiE WB" ~ "MFD WB",
    Method == "MESuSiE BB" ~ "MFD BB",
    TRUE ~ Method
  ))

ancestry_all_ROC_data_dataframe_3 <- bind_rows(ancestry_all_ROC_data_dataframe, ancestry_mfd) %>%
  load_data_with_pattern(pattern_name_3)

# ---------------------------------------------------------
# 5. FDR Power (Either & Shared)
# ---------------------------------------------------------

# --- Either ---
FDR_Power_either <- FDR_Power_either %>%
  dplyr::filter(Method != "SuSiE_weighted")

FDR_mfd_either <- FDR_Power_either %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

FDR_Power_either_3 <- bind_rows(FDR_Power_either, FDR_mfd_either) %>%
  load_data_with_pattern(pattern_name_3)

# --- Shared ---
FDR_Power_shared <- FDR_Power_shared %>%
  dplyr::filter(Method != "SuSiE_weighted")

FDR_mfd_shared <- FDR_Power_shared %>%
  filter(Method == "MESuSiE") %>%
  mutate(Method = "MFD")

FDR_Power_shared_3 <- bind_rows(FDR_Power_shared, FDR_mfd_shared) %>%
  load_data_with_pattern(pattern_name_3)

# ---------------------------------------------------------
# 6. FDR Power (Ancestry)
# ---------------------------------------------------------

FDR_Power_ancestry <- FDR_Power_ancestry %>%
  filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))

FDR_mfd_ancestry <- FDR_Power_ancestry %>%
  filter(Method %in% c("MESuSiE WB", "MESuSiE BB")) %>%
  mutate(Method = case_when(
    Method == "MESuSiE WB" ~ "MFD WB",
    Method == "MESuSiE BB" ~ "MFD BB",
    TRUE ~ Method
  ))

FDR_Power_ancestry_3 <- bind_rows(FDR_Power_ancestry, FDR_mfd_ancestry) %>%
  load_data_with_pattern(pattern_name_3)



# Combine Data Together

# ---------------------------------------------------------
# 1. Define Ordered Levels
# ---------------------------------------------------------
# Standard Method Order (MFD first, MultiSuSiE included)
method_levels_std <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged", 
                       "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")


# Ancestry Method Order (Automatically generates: "MFD WB", "MFD BB", "MESuSiE WB", etc.)
method_levels_ancestry <- unlist(lapply(method_levels_std, function(m) c(paste(m, "WB"), paste(m, "BB"))))

# Pattern Order
pattern_levels <- c(pattern_name_1, pattern_name_2, pattern_name_3)

# ---------------------------------------------------------
# 2. Combine & Factorize: Standard Datasets
# ---------------------------------------------------------
# Apply 'method_levels_std' to these
Set_Data <- bind_rows(all_Set_data_dataframe_1, all_Set_data_dataframe_2, all_Set_data_dataframe_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

Set_Power_Data <- bind_rows(set_power_summary_1, set_power_summary_2, set_power_summary_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

ROC_either <- bind_rows(either_all_ROC_data_dataframe_1, either_all_ROC_data_dataframe_2, either_all_ROC_data_dataframe_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

ROC_shared <- bind_rows(shared_all_ROC_data_dataframe_1, shared_all_ROC_data_dataframe_2, shared_all_ROC_data_dataframe_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

FDR_Power_either <- bind_rows(FDR_Power_either_1, FDR_Power_either_2, FDR_Power_either_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

FDR_Power_shared <- bind_rows(FDR_Power_shared_1, FDR_Power_shared_2, FDR_Power_shared_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_std))

# ---------------------------------------------------------
# 3. Combine & Factorize: Ancestry Datasets
# ---------------------------------------------------------
# Apply 'method_levels_ancestry' to these so "MFD WB" etc. are ordered correctly
ROC_ancestry <- bind_rows(ancestry_all_ROC_data_dataframe_1, ancestry_all_ROC_data_dataframe_2, ancestry_all_ROC_data_dataframe_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_ancestry))

FDR_Power_ancestry <- bind_rows(FDR_Power_ancestry_1, FDR_Power_ancestry_2, FDR_Power_ancestry_3) %>%
  mutate(Pattern = factor(Pattern, levels = pattern_levels),
         Method  = factor(Method, levels = method_levels_ancestry))

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
#p_size_power <- (p_size_box / p_power_bar) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#p_size_power<-p_size_power+ plot_layout(guides = "collect",heights=c(2,1))&theme(legend.position = 'bottom')


# Function to generate FDR power comparison plot (defined elsewhere)
FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method,"MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

FDR_Power_Combined<-rbind( FDR_Power_either%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Either"),
                           FDR_Power_shared%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%
                             filter(FDR==0.05)%>%mutate(Cat = "Shared"),
                           FDR_Power_ancestry%>%group_by(FDR,Method,Pattern)%>%summarise(Power= round(mean(Power),2))%>%filter(FDR==0.05)%>%mutate(Cat = "Ancestry"))



FDR_Power_Combined$Method <- factor(
  FDR_Power_Combined$Method,
  levels = c(
    "MFD", "MESuSiE", "SuSiE", "SuSiE_merged", 
    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX"
  )
)

FDR_Power_Combined$Cat<-factor(	 FDR_Power_Combined$Cat,levels = c("Either","Shared","Ancestry"))

p_FDR_Power<-FDR_Power_Compare_fun_noannot(FDR_Power_Combined,pattern_input_Power)











##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
#for (External_index in 1:2) {

External_index = 1
External_index_name <- c("", "_External")[External_index]

# Create a list to store the column plots for each iteration
plot_list <- list()

for (causal_index in 1:2) {
  
  causal_index_name <- c("Both", "One")[causal_index]
  
  # 1. Define File Prefix
  file_prefix <- paste0(
    "Missing_causal_",
    ifelse(External_index_name == "_External", "External_", ""), 
    causal_index_name
  )
  
  # 2. Conditional Data Loading (Logic adapted from source)
  if (causal_index_name == "One"){
    
    # Load MFD specific data
    env_mfd <- new.env(parent = emptyenv())
    load(paste0(res_dir, file_prefix, "_MFD.RData"), envir = env_mfd)
    
    # Extract only MFD rows
    Signal_ld_mfd       <- env_mfd$Signal_ld_data       %>% filter(Method == "MFD")
    Signal_distance_mfd <- env_mfd$Signal_distance_data %>% filter(Method == "MFD")
    
    # Load Meta data
    env_meta <- new.env(parent = emptyenv())
    load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
    
    # Combine LD Data
    Signal_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
    
    # Combine Distance Data
    Signal_distance_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
    
    # Clean up
    rm(env_mfd, env_meta)
    
  } else {
    # Logic for "Both" (Mock MFD using MESuSiE)
    env_meta <- new.env(parent = emptyenv())
    load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
    
    Signal_ld_mfd       <- env_meta$Signal_ld_data       %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
    Signal_distance_mfd <- env_meta$Signal_distance_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
    
    # Combine LD Data
    Signal_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
    
    # Combine Distance Data
    Signal_distance_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
    
    # Clean up
    rm(env_meta)
  }
  
  # 3. Factor Standardization (Crucial for consistent colors/legends)
  Signal_ld_data <- Signal_ld_data %>%
    mutate(Method = factor(Method, levels = method_levels_std))
  
  Signal_distance_data <- Signal_distance_data %>%
    mutate(Method = factor(Method, levels = method_levels_std))
  
  
  # 4. Plot Generation (Original formatting preferences)
  p_LD_signal_box <- Signal_LD_Fun(Signal_ld_data) + 
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()  # Removes the little tick marks
    ) 
  
  p_distance_signal_box <- Signal_Distance_Fun(Signal_distance_data) + 
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank() 
    ) 
  
  # Add a title to the top plot to distinguish the scenario
  p_LD_signal_box <- p_LD_signal_box + 
    ggtitle(paste("Scenario: Missing Causal SNP in", causal_index_name, "Ancestry group"))
  
  # Combine vertically and store in list
  plot_list[[causal_index]] <- p_LD_signal_box / p_distance_signal_box
}

# Combine the two columns side-by-side
# index 1 is "Both", index 2 is "One"
final_figure <- (plot_list[[2]] | plot_list[[1]]) +
  plot_layout(guides = "collect") &    # Collects legends into one
  theme(legend.position = "bottom")

# View the final combined object
final_figure

####################################################################################
#
#
#			Combine all the missing SNP figures
#
#
#
####################################################################################



p_a = p_size_box+custom_theme()
p_b = p_power_bar+custom_theme()
p_c = p_FDR_Power
p_d = final_figure+custom_theme()

blank_x <- theme(
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

drop_legends <- theme(legend.position = "none")

# Global theme: DO NOT touch axis.text / axis.text.x here
base_theme <- theme(
  axis.text.y  = element_text(size = 9),
  axis.title.y = element_text(size = 9, face = "bold"),
  strip.text.x = element_text(size = 9),
  strip.text.y = element_text(size = 9)
)

p_a <- p_a + blank_x + drop_legends
p_c <- p_c + drop_legends
p_d <- p_d + drop_legends

# p_b is the ONLY legend source; force single-row legends HERE
p_b <- p_b + blank_x +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    fill  = guide_legend(nrow = 1, byrow = TRUE, order = 2),
    shape = guide_legend(nrow = 1, byrow = TRUE, order = 3)
    # Disable the ggpattern legend if present:
    # pattern = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

row_top    <- p_a + p_b
row_bottom <- p_c + p_d + plot_layout(widths = c(1, 1))

final <- row_top / row_bottom +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    )
  ) &
  base_theme

ggsave(
  paste0(plot_dir, "Missing_simulation_summary_updated.pdf"),
  final,
  height = 300, width = 500, units = "mm", dpi = 600
)

