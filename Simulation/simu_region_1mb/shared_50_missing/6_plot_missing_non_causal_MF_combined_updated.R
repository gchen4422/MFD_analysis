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




res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_MFD_combined"
res_out_1<-"Missing_non_causal_One_MFD.RData"
res_out_2<-"Missing_non_causal_Both_MFD.RData"
res_out_1_external = "Missing_non_causal_One_meta_updated.RData"
res_out_2_external = "Missing_non_causal_Both_meta_updated.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 




Compare_Plot_Missing <- function(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) {
  
  
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
  
  
  library(openxlsx)
  library(stringr)
  
  res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
  res_out_xlsx <- paste0(plot_dir_name,".xlsx")
  full_path <- paste0(res_dir_table, res_out_xlsx)
  sheet_name_raw <- str_remove(res_out_xlsx, ".xlsx")
  sheet_name <- substr(sheet_name_raw, 1, 31) 
  
  # 1. Define ALL Orders (Method, Pattern, Signal Type)
  method_order_base <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged", 
                         "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
  
  # Generate Ancestry-Specific Levels (e.g., "MESuSiE WB", "MESuSiE BB")
  method_order_ancestry <- unlist(lapply(method_order_base, function(m) paste0(m, c(" WB", " BB"))))
  master_method_levels <- c(method_order_base, method_order_ancestry)
  
  pattern_order <- c("None Missing","Missing in One Ancestry", "Missing in Both Ancestry")
  signal_type_order <- c("Either", "Shared", "Ancestry")
  
  
  # 2. Process and Sort the Data
  combined_data <- bind_rows(
    "Either"   = FDR_Power_either,
    "Shared"   = FDR_Power_shared,
    "Ancestry" = FDR_Power_ancestry,
    .id = "Signal_Type" 
  ) %>%
    select(Method, h2,causal_num, Signal_Type, Pattern, FDR, Power) %>%
    
    # Apply Factors for ALL sorting columns
    mutate(
      Method = factor(Method, levels = master_method_levels),
      Pattern = factor(Pattern, levels = pattern_order),
      Signal_Type = factor(Signal_Type, levels = signal_type_order) # Enforces Either -> Shared -> Ancestry
    ) %>%
    
    # Order methods for reporting
    arrange(Pattern,Signal_Type, Method)
  
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
  p_power_bar<-Set_Power_Compare_fun_supp(Set_Power_Data,pattern_input_Power)
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
  p_ROC_either <- ROC_ancestry_Compare_fun(ROC_either%>%mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")), pattern_input_line)
  
  # Create and save the combined plot for either PIP
  p_ROC_FDR_Power_either <- (p_ROC_either + p_FDR_Power_either) + 
    plot_annotation(tag_levels = 'a') & 
    theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_either <- p_ROC_FDR_Power_either + plot_layout(heights = c(2, 1))
  
  ggsave(paste0(plot_dir, "PIP_ROC_FDR_Power_Either_Compare.pdf"), p_ROC_FDR_Power_either, height = 180, width = 210, units = "mm", dpi = 600)
  
  ###Plot for Shared PIP
  
  p_FDR_Power_shared<-FDR_Power_Compare_fun_supp(FDR_Power_shared,pattern_input_Power)
  p_ROC_shared<-ROC_ancestry_Compare_fun(ROC_shared%>%mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")),pattern_input_line )
  
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
                                               mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")),pattern_input_Power)
  
  

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
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")),pattern_input_line)

    
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
                                               mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")),pattern_input_Power)
  
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
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MFD","MESuSiE","SuSiE","SuSiE_merged","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX")),pattern_input_line)
  
  p_ROC_FDR_Power_BB<-(p_ROC_BB + p_FDR_Power_BB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_BB<-p_ROC_FDR_Power_BB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_BB_Compare.pdf"),p_ROC_FDR_Power_BB,height=180, width=210, units = "mm",dpi=600)
  
}



res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_MFD_combined"
res_out_1<-"Missing_non_causal_One_MFD.RData"
res_out_2<-"Missing_non_causal_Both_MFD.RData"
res_out_1_external = "Missing_non_causal_One_meta_updated.RData"
res_out_2_external = "Missing_non_causal_Both_meta_updated.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
#Compare_Plot_Missing(res_dir, plot_dir_name, res_out_1, res_out_2,pattern_name_1, pattern_name_2,pattern_name_3) 


res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_Ref_MFD_combined"
res_out_1<-"Missing_non_causal_External_One_MFD.RData"
res_out_2<-"Missing_non_causal_External_Both_MFD.RData"
res_out_1_external = "Missing_non_causal_External_One_meta_updated.RData"
res_out_2_external = "Missing_non_causal_External_Both_meta_updated.RData"
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
for (External_index in 1:2) {
  External_index_name <- c("", "_External")[External_index]
  
  plot_dir <- paste0(
    res_dir,
    "Figure/Missing_Causal", External_index_name, "_MFD_combined/"
  )
  system(paste0("mkdir -p ", plot_dir))
  
  for (causal_index in 1:2) {
    
    
    
    causal_index_name <- c("Both", "One")[causal_index]
    file_prefix <- paste0(
      "Missing_causal_",
      ifelse(External_index_name == "_External", "External_", ""), 
      causal_index_name
    )
    if (causal_index_name == "One"){

      env_mfd <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_MFD.RData"), envir = env_mfd)
      
      # Extract only MFD rows
      Signal_ld_mfd       <- env_mfd$Signal_ld_data       %>% filter(Method == "MFD")
      Signal_distance_mfd <- env_mfd$Signal_distance_data %>% filter(Method == "MFD")
      
      env_meta <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
      
      # Combine LD Data
      Signal_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
  
      # Combine Distance Data
      Signal_distance_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
      
      # Clean up temporary environments to save memory
      rm(env_mfd, env_meta)
    } else {
      env_meta <- new.env(parent = emptyenv())
      
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
      
      Signal_ld_mfd       <- env_meta$Signal_ld_data       %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
      Signal_distance_mfd <- env_meta$Signal_distance_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
      
      # Combine LD Data
      Signal_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
      
      # Combine Distance Data
      Signal_distance_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
      
      # Clean up temporary environments to save memory
      rm(env_meta)
      
    }
    
    # 2. Apply to Signal_ld_data
    Signal_ld_data <- Signal_ld_data %>%
      mutate(Method = factor(Method, levels = method_levels_std))
    
    # 3. Apply to Signal_distance_data
    Signal_distance_data <- Signal_distance_data %>%
      mutate(Method = factor(Method, levels = method_levels_std))
    
    
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




res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
# 1. Initialize lists to store combined data
all_ld_data_list   <- list()
all_dist_data_list <- list()

# Define levels once (assuming these variables are defined in the environment)
# master_method_levels <- ... 

for (External_index in 1:2) {
  External_index_name <- c("", "_External")[External_index]
  
  # Optional: Keep directory creation when plots will be saved later
  plot_dir <- paste0(res_dir, "Figure/Missing_Causal", External_index_name, "_MFD_combined/")
  system(paste0("mkdir -p ", plot_dir))
  
  for (causal_index in 1:2) {
    
    causal_index_name <- c("One", "Both")[causal_index]
    
    file_prefix <- paste0(
      "Missing_causal_",
      ifelse(External_index_name == "_External", "External_", ""), 
      causal_index_name
    )
    
    # --- DATA LOADING & PROCESSING LOGIC ---
    if (causal_index_name == "One") {
      
      env_mfd <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_MFD.RData"), envir = env_mfd)
      
      # Extract only MFD rows
      Signal_ld_mfd       <- env_mfd$Signal_ld_data       %>% filter(Method == "MFD")
      Signal_distance_mfd <- env_mfd$Signal_distance_data %>% filter(Method == "MFD")
      
      env_meta <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
      
      # Combine LD Data
      current_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
      
      # Combine Distance Data
      current_dist_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
      
      rm(env_mfd, env_meta)
      
    } else {
      # For "Both" case
      env_meta <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)
      
      Signal_ld_mfd       <- env_meta$Signal_ld_data       %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
      Signal_distance_mfd <- env_meta$Signal_distance_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD")
      
      current_ld_data <- bind_rows(env_meta$Signal_ld_data, Signal_ld_mfd) 
      current_dist_data <- bind_rows(env_meta$Signal_distance_data, Signal_distance_mfd) 
      
      rm(env_meta)
    }
    
    # --- TAGGING AND ACCUMULATION ---
    
    # Add columns to identify which loop iteration this data belongs to
    current_ld_data <- current_ld_data %>%
      mutate(
        External_Type = ifelse(External_index_name == "", "Standard", "External"),
        Causal_Type   = causal_index_name # "One" or "Both"
      )
    
    current_dist_data <- current_dist_data %>%
      mutate(
        External_Type = ifelse(External_index_name == "", "Standard", "External"),
        Causal_Type   = causal_index_name
      )
    
    # Append to lists
    all_ld_data_list[[paste(External_index, causal_index, sep="_")]]   <- current_ld_data
    all_dist_data_list[[paste(External_index, causal_index, sep="_")]] <- current_dist_data
  }
}

causal_order   <- c("One", "Both")        
external_order <- c("Standard", "External")
  
  
# 2. Final Combination

Global_Signal_ld_data <- bind_rows(all_ld_data_list) %>%
  mutate(
    Method        = factor(Method, levels = method_levels_std),
    Causal_Type   = factor(Causal_Type, levels = causal_order),
    External_Type = factor(External_Type, levels = external_order)
  ) %>%
  select(Method, h2, External_Type, Causal_Type, Cor) %>%
  
  # --- NEW: Group and Average ---
  group_by(Causal_Type, External_Type, Method, h2) %>%
  summarise(
    Cor = mean(Cor, na.rm = TRUE),
    .groups = "drop" # Ungroups data automatically after summarizing
  ) %>%
  # ------------------------------

arrange(Causal_Type, External_Type, Method)

# Convert to tibble (optional, as summarize usually returns a tibble)
Global_Signal_ld_data <- as_tibble(Global_Signal_ld_data) %>% select(Method, h2, External_Type, Causal_Type, Cor)

# 3. Refine Distance Data (Sorted by: Method -> Causal -> External)
Global_Signal_distance_data <- bind_rows(all_dist_data_list) %>%
  mutate(
    Method        = factor(Method, levels = method_levels_std),
    Causal_Type   = factor(Causal_Type, levels = causal_order),
    External_Type = factor(External_Type, levels = external_order)
  ) %>%
  select(Method, h2, External_Type, Causal_Type, Distance) %>%
  # This will now respect the factor levels defined above
  arrange(Causal_Type, External_Type,Method)

# 3. Check the result
print(table(Global_Signal_ld_data$External_Type, Global_Signal_ld_data$Causal_Type))

library(openxlsx)
library(stringr)

res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
res_out_xlsx <- "Missing_Causal_LD_Distance_Combined.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

# Initialize Workbook
wb <- createWorkbook()

# Define Styles
title_style <- createStyle(textDecoration = "bold", fontSize = 12)
header_style <- createStyle(textDecoration = "bold", border = "Bottom")

# --- Sheet 1: LD Analysis ---
sheet1 <- "LD Analysis"
addWorksheet(wb, sheet1)
writeData(wb, sheet1, "Supplementary Table: Signal LD Analysis (Combined)", startRow = 1)
addStyle(wb, sheet1, title_style, rows = 1, cols = 1)
writeData(wb, sheet1, Global_Signal_ld_data, startRow = 2, headerStyle = header_style)

# --- Sheet 2: Distance Analysis ---
sheet2 <- "Distance Analysis"
addWorksheet(wb, sheet2)
writeData(wb, sheet2, "Supplementary Table: Signal Distance Analysis (Combined)", startRow = 1)
addStyle(wb, sheet2, title_style, rows = 1, cols = 1)
writeData(wb, sheet2, Global_Signal_distance_data, startRow = 2, headerStyle = header_style)

# Save
saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined Excel file to:", full_path))


########################################
test =Set_Data%>%group_by(Method,Pattern)%>%summarise(round(median(Size)),round(mean(Power),2))


Signal_ld_data%>%group_by(Method,h2)%>%summarise(round(mean(Cor),3))
Signal_distance_data%>%group_by(Method,h2)%>%summarise(round(mean(Distance )))


all_Set_data_dataframe%>%group_by(Method,h2,causal_num)%>%summarise(round(median(Size)),round(mean(Power),2))
all_Set_data_dataframe%>%group_by(Method)%>%summarise(round(median(Size)),round(mean(Power),2))


test = all_Set_data_all%>%group_by(Method,Pattern,flip,rho)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2))
test = all_Set_data_all%>%group_by(Method,Pattern,rho)%>%summarise(round(median(Size)),round(mean(Power),2),round(mean(FDR),2)) 


FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2))
FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2)) %>% filter(Method %in% c("MESuSiE","SuSiE","MFD"))
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2)) %>% filter(Pattern == "Missing in One Ancestry")
FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2))%>% filter(Method %in% c("MESuSiE WB","MESuSiE BB","MFD WB","MFD BB","SuSiE BB","SuSiE WB"))
