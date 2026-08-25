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




Compare_Plot_Missing_Revision <- function(res_dir, plot_dir_name, res_out_1, res_out_2, pattern_name_1, pattern_name_2, pattern_name_3) {

  # MFD methods to keep
  mfd_methods <- c("MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")
  mfd_methods_ancestry <- c("MFD_MESuSiE WB", "MFD_MESuSiE BB",
                             "MFD_SuSiEx WB", "MFD_SuSiEx BB",
                             "MFD_XMAP WB", "MFD_XMAP BB")

  pattern_name_1 = "Missing in One Ancestry"

  # ---------------------------------------------------------
  # 1. LOAD Pattern 1 DATA (Missing in One) from endpoint-sensitivity file
  # ---------------------------------------------------------
  load(paste0(res_dir, res_out_1))

  all_Set_data_dataframe_1 <- load_data_with_pattern(all_Set_data_dataframe, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  set_power_summary_1 <- load_data_with_pattern(set_power_summary, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  # ROC Data
  either_all_ROC_data_dataframe_1 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  shared_all_ROC_data_dataframe_1 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  ancestry_all_ROC_data_dataframe_1 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_1) %>%
    filter(Method %in% mfd_methods_ancestry)

  # FDR Power Data
  FDR_Power_either_1 <- load_data_with_pattern(FDR_Power_either, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  FDR_Power_shared_1 <- load_data_with_pattern(FDR_Power_shared, pattern_name_1) %>%
    filter(Method %in% mfd_methods)

  FDR_Power_ancestry_1 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_1) %>%
    filter(Method %in% mfd_methods_ancestry)


  pattern_name_2 = "Missing in Both Ancestry"

  # ---------------------------------------------------------
  # 2. LOAD Pattern 2 DATA (Missing in Both) from endpoint-sensitivity file
  # ---------------------------------------------------------
  load(paste0(res_dir, res_out_2))

  all_Set_data_dataframe_2 <- load_data_with_pattern(all_Set_data_dataframe, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  set_power_summary_2 <- load_data_with_pattern(set_power_summary, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  # ROC Data
  either_all_ROC_data_dataframe_2 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  shared_all_ROC_data_dataframe_2 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  ancestry_all_ROC_data_dataframe_2 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_2) %>%
    filter(Method %in% mfd_methods_ancestry) %>% arrange(h2, causal_num, Pattern, Method)

  # FDR Power Data
  FDR_Power_either_2 <- load_data_with_pattern(FDR_Power_either, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  FDR_Power_shared_2 <- load_data_with_pattern(FDR_Power_shared, pattern_name_2) %>%
    filter(Method %in% mfd_methods) %>% arrange(h2, causal_num, Pattern, Method)

  FDR_Power_ancestry_2 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_2) %>%
    filter(Method %in% mfd_methods_ancestry) %>% arrange(h2, causal_num, Pattern, Method)



  # ---------------------------------------------------------
  # 3. LOAD Pattern 3 DATA (None Missing)
  #    Clone from baseline: MESuSiE -> MFD_MESuSiE, SuSiEx -> MFD_SuSiEx, XMAP -> MFD_XMAP
  # ---------------------------------------------------------

  pattern_name_3 = "None Missing"

  load(paste0(res_dir, "shared_50_baseline_updated_meta_updated_xmap.RData"))

  # Helper: clone standard methods to MFD variants
  clone_mfd_std <- function(df) {
    bind_rows(
      df %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD_MESuSiE"),
      df %>% filter(Method == "SuSiEx") %>% mutate(Method = "MFD_SuSiEx"),
      df %>% filter(Method == "XMAP") %>% mutate(Method = "MFD_XMAP")
    )
  }

  # Helper: clone ancestry methods to MFD ancestry variants
  clone_mfd_ancestry <- function(df) {
    bind_rows(
      df %>% filter(Method == "MESuSiE WB") %>% mutate(Method = "MFD_MESuSiE WB"),
      df %>% filter(Method == "MESuSiE BB") %>% mutate(Method = "MFD_MESuSiE BB"),
      df %>% filter(Method == "SuSiEx WB") %>% mutate(Method = "MFD_SuSiEx WB"),
      df %>% filter(Method == "SuSiEx BB") %>% mutate(Method = "MFD_SuSiEx BB"),
      df %>% filter(Method == "XMAP WB") %>% mutate(Method = "MFD_XMAP WB"),
      df %>% filter(Method == "XMAP BB") %>% mutate(Method = "MFD_XMAP BB")
    )
  }

  # --- Set Power Summary ---
  set_power_summary <- set_power_summary %>%
    filter(Method != "SuSiE_meta_weighted") %>%
    mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))

  set_power_summary_3 <- clone_mfd_std(set_power_summary) %>%
    load_data_with_pattern(pattern_name_3) %>%
    arrange(h2, causal_num, Method)

  # --- All Set Data ---
  all_Set_data_dataframe <- all_Set_data_dataframe %>%
    filter(Method != "SuSiE_meta_weighted") %>%
    mutate(Method = if_else(Method == "SuSiE_meta_merged", "SuSiE_merged", Method))

  all_Set_data_dataframe_3 <- clone_mfd_std(all_Set_data_dataframe) %>%
    load_data_with_pattern(pattern_name_3) %>%
    dplyr::select(-any_of("FDR")) %>%
    arrange(h2, causal_num, Method)

  # --- Either ROC ---
  either_all_ROC_data_dataframe <- either_all_ROC_data_dataframe %>%
    filter(Method != "SuSiE_weighted")

  either_all_ROC_data_dataframe_3 <- clone_mfd_std(either_all_ROC_data_dataframe) %>%
    load_data_with_pattern(pattern_name_3)

  # --- Shared ROC ---
  shared_all_ROC_data_dataframe <- shared_all_ROC_data_dataframe %>%
    filter(Method != "SuSiE_weighted")

  shared_all_ROC_data_dataframe_3 <- clone_mfd_std(shared_all_ROC_data_dataframe) %>%
    load_data_with_pattern(pattern_name_3)

  # --- Ancestry ROC ---
  ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%
    filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))

  ancestry_all_ROC_data_dataframe_3 <- clone_mfd_ancestry(ancestry_all_ROC_data_dataframe) %>%
    load_data_with_pattern(pattern_name_3)

  # --- FDR Power Either ---
  FDR_Power_either <- FDR_Power_either %>%
    dplyr::filter(Method != "SuSiE_weighted")

  FDR_Power_either_3 <- clone_mfd_std(FDR_Power_either) %>%
    load_data_with_pattern(pattern_name_3)

  # --- FDR Power Shared ---
  FDR_Power_shared <- FDR_Power_shared %>%
    dplyr::filter(Method != "SuSiE_weighted")

  FDR_Power_shared_3 <- clone_mfd_std(FDR_Power_shared) %>%
    load_data_with_pattern(pattern_name_3)

  # --- FDR Power Ancestry ---
  FDR_Power_ancestry <- FDR_Power_ancestry %>%
    filter(!Method %in% c("SuSiE_weighted WB", "SuSiE_weighted BB"))

  FDR_Power_ancestry_3 <- clone_mfd_ancestry(FDR_Power_ancestry) %>%
    load_data_with_pattern(pattern_name_3)




  # Combine Data Together

  # ---------------------------------------------------------
  # 4. Define Ordered Levels
  # ---------------------------------------------------------
  method_levels_std <- mfd_methods

  # Ancestry Method Order
  method_levels_ancestry <- unlist(lapply(method_levels_std, function(m) c(paste(m, "WB"), paste(m, "BB"))))

  # Pattern Order
  pattern_levels <- c(pattern_name_1, pattern_name_2, pattern_name_3)

  # ---------------------------------------------------------
  # 5. Combine & Factorize: Standard Datasets
  # ---------------------------------------------------------
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
  # 6. Combine & Factorize: Ancestry Datasets
  # ---------------------------------------------------------
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
  method_order_base <- mfd_methods

  # Generate Ancestry-Specific Levels
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
      Signal_Type = factor(Signal_Type, levels = signal_type_order)
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
  writeData(wb, sheet_name, "Supplementary Table: FDR Controlled Power Analysis (MFD Revision Sensitivity)", startRow = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)

  # Write the Combined Data Table starting at Row 2
  writeData(wb, sheet_name, combined_data, startRow = 2)

  # 5. Save
  saveWorkbook(wb, full_path, overwrite = TRUE)

  print(paste("Saved combined table to:", full_path))

  # Define a directory to save the plots
  plot_dir <- paste0(res_dir, "Figure/",plot_dir_name,"/")
  system(paste0("mkdir -p ",plot_dir))
  # Define patterns for plotting
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

  # Function to generate FDR power comparison plot
  p_FDR_Power_either <- FDR_Power_Compare_fun_supp(FDR_Power_either, pattern_input_Power)

  # Function to generate ROC comparison plot
  p_ROC_either <- ROC_ancestry_Compare_fun(ROC_either%>%mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")), pattern_input_line)

  # Create and save the combined plot for either PIP
  p_ROC_FDR_Power_either <- (p_ROC_either + p_FDR_Power_either) +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_either <- p_ROC_FDR_Power_either + plot_layout(heights = c(2, 1))

  ggsave(paste0(plot_dir, "PIP_ROC_FDR_Power_Either_Compare.pdf"), p_ROC_FDR_Power_either, height = 180, width = 210, units = "mm", dpi = 600)

  ###Plot for Shared PIP

  p_FDR_Power_shared<-FDR_Power_Compare_fun_supp(FDR_Power_shared,pattern_input_Power)
  p_ROC_shared<-ROC_ancestry_Compare_fun(ROC_shared%>%mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")),pattern_input_line )

  p_ROC_FDR_Power_shared<-(p_ROC_shared + p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_shared<-p_ROC_FDR_Power_shared+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared_Compare.pdf"),p_ROC_FDR_Power_shared,height=180, width=210, units = "mm",dpi=600)

  ###Plot for WB specific

  # Helper to strip ancestry suffix from MFD method names
  strip_ancestry_suffix <- function(method) {
    case_when(
      str_detect(method, "^MFD_MESuSiE(\\b|_)") ~ "MFD_MESuSiE",
      str_detect(method, "^MFD_SuSiEx(\\b|_)") ~ "MFD_SuSiEx",
      str_detect(method, "^MFD_XMAP(\\b|_)") ~ "MFD_XMAP",
      TRUE ~ method
    )
  }

  p_FDR_Power_WB<-FDR_Power_Compare_fun_supp(FDR_Power_ancestry[grep("WB",FDR_Power_ancestry$Method),] %>%
                                               mutate(Method = strip_ancestry_suffix(Method))%>%  ungroup() %>%
                                               mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")),pattern_input_Power)


  p_ROC_WB<-ROC_ancestry_Compare_fun(ROC_ancestry[grep("WB",ROC_ancestry$Method),]%>%
                                       mutate(Method = strip_ancestry_suffix(Method))%>% ungroup() %>%
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")),pattern_input_line)


  p_ROC_FDR_Power_WB<-(p_ROC_WB + p_FDR_Power_WB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_WB<-p_ROC_FDR_Power_WB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_WB_Compare.pdf"),p_ROC_FDR_Power_WB,height=180, width=210, units = "mm",dpi=600)


  ###Plot for BB specific
  p_FDR_Power_BB<-FDR_Power_Compare_fun_supp(FDR_Power_ancestry[grep("BB",FDR_Power_ancestry$Method),] %>%
                                               mutate(Method = strip_ancestry_suffix(Method))%>% ungroup() %>%
                                               mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")),pattern_input_Power)

  p_ROC_BB<-ROC_ancestry_Compare_fun(ROC_ancestry[grep("BB",ROC_ancestry$Method),]%>%
                                       mutate(Method = strip_ancestry_suffix(Method))%>% ungroup() %>%
                                       mutate(Method = as_factor(Method) %>%fct_relevel("MFD_MESuSiE","MFD_SuSiEx","MFD_XMAP")),pattern_input_line)

  p_ROC_FDR_Power_BB<-(p_ROC_BB + p_FDR_Power_BB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_BB<-p_ROC_FDR_Power_BB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_BB_Compare.pdf"),p_ROC_FDR_Power_BB,height=180, width=210, units = "mm",dpi=600)

}



res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_MFD_endpoint_sensitivity"
res_out_1<-"Missing_non_causal_One_MFD_endpoint_sensitivity.RData"
res_out_2<-"Missing_non_causal_Both_MFD_endpoint_sensitivity.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
Compare_Plot_Missing_Revision(res_dir, plot_dir_name, res_out_1, res_out_2, pattern_name_1, pattern_name_2, pattern_name_3)


res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name<-"Missing_Non_Causal_Ref_MFD_endpoint_sensitivity"
res_out_1<-"Missing_non_causal_External_One_MFD_endpoint_sensitivity.RData"
res_out_2<-"Missing_non_causal_External_Both_MFD_endpoint_sensitivity.RData"
pattern_name_1 = "Missing in One Ancestry"
pattern_name_2 = "Missing in Both Ancestry"
pattern_name_3 = "None Missing"
Compare_Plot_Missing_Revision(res_dir, plot_dir_name, res_out_1, res_out_2, pattern_name_1, pattern_name_2, pattern_name_3)


########################################
# Summary statistics
########################################
#test =Set_Data%>%group_by(Method,Pattern)%>%summarise(round(median(Size)),round(mean(Power),2))
#FDR_Power_either%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2))
#FDR_Power_shared%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2))
#FDR_Power_ancestry%>%filter(FDR==0.05)%>%group_by(Method,Pattern)%>%summarise(round(mean(Power),2))


##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

mfd_methods_causal <- c("MFD_MESuSiE", "MFD_SuSiEx", "MFD_XMAP")

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
for (External_index in 1:2) {
  External_index_name <- c("", "_External")[External_index]

  plot_dir <- paste0(
    res_dir,
    "Figure/Missing_Causal", External_index_name, "_MFD_endpoint_sensitivity/"
  )
  system(paste0("mkdir -p ", plot_dir))

  for (causal_index in 1:2) {

    causal_index_name <- c("Both", "One")[causal_index]
    file_prefix <- paste0(
      "Missing_causal_",
      ifelse(External_index_name == "_External", "External_", ""),
      causal_index_name
    )

    if (causal_index_name == "One") {
      # Load the endpoint sensitivity file directly — already has MFD variants
      env_rev <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_MFD_endpoint_sensitivity.RData"), envir = env_rev)

      Signal_ld_data <- env_rev$Signal_ld_data
      Signal_distance_data <- env_rev$Signal_distance_data

      rm(env_rev)

    } else {
      # For "Both": clone from _updated_meta.RData
      # MFD always picks multi-ancestry for "Both", so MFD_X ≈ X
      env_meta <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)

      Signal_ld_data <- bind_rows(
        env_meta$Signal_ld_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD_MESuSiE"),
        env_meta$Signal_ld_data %>% filter(Method == "SuSiEx") %>% mutate(Method = "MFD_SuSiEx"),
        env_meta$Signal_ld_data %>% filter(Method == "XMAP") %>% mutate(Method = "MFD_XMAP")
      )

      Signal_distance_data <- bind_rows(
        env_meta$Signal_distance_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD_MESuSiE"),
        env_meta$Signal_distance_data %>% filter(Method == "SuSiEx") %>% mutate(Method = "MFD_SuSiEx"),
        env_meta$Signal_distance_data %>% filter(Method == "XMAP") %>% mutate(Method = "MFD_XMAP")
      )

      rm(env_meta)
    }

    # Apply factor levels
    Signal_ld_data <- Signal_ld_data %>%
      mutate(Method = factor(Method, levels = mfd_methods_causal))

    Signal_distance_data <- Signal_distance_data %>%
      mutate(Method = factor(Method, levels = mfd_methods_causal))

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



# Combined supplementary table for Missing Causal
res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
all_ld_data_list   <- list()
all_dist_data_list <- list()

for (External_index in 1:2) {
  External_index_name <- c("", "_External")[External_index]

  for (causal_index in 1:2) {

    causal_index_name <- c("One", "Both")[causal_index]

    file_prefix <- paste0(
      "Missing_causal_",
      ifelse(External_index_name == "_External", "External_", ""),
      causal_index_name
    )

    if (causal_index_name == "One") {
      env_rev <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_MFD_endpoint_sensitivity.RData"), envir = env_rev)

      current_ld_data <- env_rev$Signal_ld_data
      current_dist_data <- env_rev$Signal_distance_data

      rm(env_rev)

    } else {
      env_meta <- new.env(parent = emptyenv())
      load(paste0(res_dir, file_prefix, "_updated_meta.RData"), envir = env_meta)

      current_ld_data <- bind_rows(
        env_meta$Signal_ld_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD_MESuSiE"),
        env_meta$Signal_ld_data %>% filter(Method == "SuSiEx") %>% mutate(Method = "MFD_SuSiEx"),
        env_meta$Signal_ld_data %>% filter(Method == "XMAP") %>% mutate(Method = "MFD_XMAP")
      )

      current_dist_data <- bind_rows(
        env_meta$Signal_distance_data %>% filter(Method == "MESuSiE") %>% mutate(Method = "MFD_MESuSiE"),
        env_meta$Signal_distance_data %>% filter(Method == "SuSiEx") %>% mutate(Method = "MFD_SuSiEx"),
        env_meta$Signal_distance_data %>% filter(Method == "XMAP") %>% mutate(Method = "MFD_XMAP")
      )

      rm(env_meta)
    }

    # Tag with External_Type and Causal_Type
    current_ld_data <- current_ld_data %>%
      mutate(
        External_Type = ifelse(External_index_name == "", "Standard", "External"),
        Causal_Type   = causal_index_name
      )

    current_dist_data <- current_dist_data %>%
      mutate(
        External_Type = ifelse(External_index_name == "", "Standard", "External"),
        Causal_Type   = causal_index_name
      )

    all_ld_data_list[[paste(External_index, causal_index, sep="_")]]   <- current_ld_data
    all_dist_data_list[[paste(External_index, causal_index, sep="_")]] <- current_dist_data
  }
}

causal_order   <- c("One", "Both")
external_order <- c("Standard", "External")

Global_Signal_ld_data <- bind_rows(all_ld_data_list) %>%
  mutate(
    Method        = factor(Method, levels = mfd_methods_causal),
    Causal_Type   = factor(Causal_Type, levels = causal_order),
    External_Type = factor(External_Type, levels = external_order)
  ) %>%
  select(Method, h2, External_Type, Causal_Type, Cor) %>%
  group_by(Causal_Type, External_Type, Method, h2) %>%
  summarise(
    Cor = mean(Cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Causal_Type, External_Type, Method)

Global_Signal_ld_data <- as_tibble(Global_Signal_ld_data) %>% select(Method, h2, External_Type, Causal_Type, Cor)

Global_Signal_distance_data <- bind_rows(all_dist_data_list) %>%
  mutate(
    Method        = factor(Method, levels = mfd_methods_causal),
    Causal_Type   = factor(Causal_Type, levels = causal_order),
    External_Type = factor(External_Type, levels = external_order)
  ) %>%
  select(Method, h2, External_Type, Causal_Type, Distance) %>%
  arrange(Causal_Type, External_Type, Method)

print(table(Global_Signal_ld_data$External_Type, Global_Signal_ld_data$Causal_Type))

library(openxlsx)

res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
res_out_xlsx <- "Missing_Causal_LD_Distance_MFD_endpoint_sensitivity.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

wb <- createWorkbook()

title_style <- createStyle(textDecoration = "bold", fontSize = 12)
header_style <- createStyle(textDecoration = "bold", border = "Bottom")

sheet1 <- "LD Analysis"
addWorksheet(wb, sheet1)
writeData(wb, sheet1, "Supplementary Table: Signal LD Analysis (MFD Revision Sensitivity)", startRow = 1)
addStyle(wb, sheet1, title_style, rows = 1, cols = 1)
writeData(wb, sheet1, Global_Signal_ld_data, startRow = 2, headerStyle = header_style)

sheet2 <- "Distance Analysis"
addWorksheet(wb, sheet2)
writeData(wb, sheet2, "Supplementary Table: Signal Distance Analysis (MFD Revision Sensitivity)", startRow = 1)
addStyle(wb, sheet2, title_style, rows = 1, cols = 1)
writeData(wb, sheet2, Global_Signal_distance_data, startRow = 2, headerStyle = header_style)

saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined Excel file to:", full_path))
