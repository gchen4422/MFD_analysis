library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)
library(ggbreak)
library(DescTools)
library(coin)
library(susieR)
library(ggrepel)
library(stringr)

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summary_res/"
load(paste0(res_dir,"res_all_mf_updated.RData"))   

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region+ 204, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region+ 243, res_all$Region)
res_all$CS <- ifelse(res_all$CS != 0, 1, res_all$CS)
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
res_all$ASV <- ifelse(is.na(res_all$MAF_eur) | is.na(res_all$MAF_afr), 1, 0)



classify_signal <- function(shared, p1, p2) {
  case_when(
    shared > 0.5 ~ 1,
    p1 > 0.5 & p2 < 0.5 ~ 2,
    p1 < 0.5 & p2 > 0.5 ~ 3,
    .default = 0
  )
}

res_all <- res_all %>%
  mutate(
    MFD_Signal     = classify_signal(PIP_Shared, PIP_Ancestry_1, PIP_Ancestry_2),
    SuSiE_Signal   = classify_signal(SuSiE_PIP_Shared, SuSiE_PIP_EU, SuSiE_PIP_BB),
    MESuSiE_Signal = classify_signal(MESuSiE_PIP_Shared, MESuSiE_PIP_WB, MESuSiE_PIP_BB)
  )

signal_cols <- c("MFD_Signal", "SuSiE_Signal", "MESuSiE_Signal")

# Print tables for ALL data
print("--- Counts for ALL Data ---")
res_all %>% 
  select(all_of(signal_cols)) %>% 
  apply(2, table)

# Print tables for ASV == 1 subset
print("--- Counts for ASV == 1 Subset ---")
res_all %>% 
  filter(ASV == 1) %>% 
  select(all_of(signal_cols)) %>% 
  apply(2, table)

#########################################################


####Combined Shared vs. Unique plot PIP 0.5 to 0.9


##########################################################

res_all = res_all %>% mutate(MFD_PIP_Unique = pmax(PIP_Ancestry_1,PIP_Ancestry_2),
                             SuSiE_PIP_Unique = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),
                             MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_WB,MESuSiE_PIP_BB),
                             CARMAX_PIP_Unique = pmax(CARMAX_PIP_WB,CARMAX_PIP_BB))

res_all %>% filter(MFD_PIP_Unique > 0.5) %>% filter(!(MFD_Signal %in% c(2,3))) %>% mutate(MFD_PIP_Unique = 0)
res_all %>% filter(SuSiE_PIP_Unique > 0.5) %>% filter(!(SuSiE_Signal %in% c(2,3)))%>% mutate(SuSiE_PIP_Unique = 0)


res_all <- res_all %>%
  mutate(
    MFD_PIP_Unique = case_when(
      MFD_PIP_Unique > 0.5 & !(MFD_Signal %in% c(2, 3)) ~ 0,
      TRUE ~ MFD_PIP_Unique # Keep original value otherwise
    ),
    
    SuSiE_PIP_Unique = case_when(
      SuSiE_PIP_Unique > 0.5 & !(SuSiE_Signal %in% c(2, 3)) ~ 0,
      TRUE ~ SuSiE_PIP_Unique # Keep original value otherwise
    )
  )

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    strip.background = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7, face="bold"),
    plot.title = element_text(size=7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

# Define annotation names and thresholds
ann_col_name<-c( "non_synonymous_common","synonymous", "promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")


#counts_df <- res_all %>%
#  select(contains("PIP_Shared")) %>%
#  summarise(across(everything(), ~ sum(. > 0.5, na.rm = TRUE))) %>%
#  pivot_longer(everything(), names_to = "Method_Column", values_to = "Count_gt_0.5")

#print(counts_df)



#counts_df <- res_all %>%
#  select(contains("PIP_Unique")) %>%
#  summarise(across(everything(), ~ sum(. > 0.5, na.rm = TRUE))) %>%
#  pivot_longer(everything(), names_to = "Method_Column", values_to = "Count_gt_0.5")

#print(counts_df)

#counts_df <- res_all %>%filter(ASV==1) %>% 
#  select(contains("PIP_Unique"))  %>% 
#  summarise(across(everything(), ~ sum(. > 0.5, na.rm = TRUE))) %>%
#  pivot_longer(everything(), names_to = "Method_Column", values_to = "Count_gt_0.5")

#print(counts_df)

top_n_steps <- seq(10, 200, by = 10)
binary_ann  <- ann_col_name
cont_ann    <- c("Z_eur")

# Pre-calculate TOTAL sums for binary annotations (do this once globally)
total_n <- nrow(res_all)
total_sums_bin <- res_all %>%
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))

# Define the Column Map
pip_config <- list(
  Shared = c(
    MFD     = "PIP_Shared",
    MESuSiE = "MESuSiE_PIP_Shared",
    SuSiE   = "SuSiE_PIP_Shared"
  ),
  Unique = c(
    MFD     = "MFD_PIP_Unique",
    MESuSiE = "MESuSiE_PIP_Unique",
    SuSiE   = "SuSiE_PIP_Unique"
  )
)

# 2. Compute Enrichment (Fast Method) ------------------------------------------

enrichment_results_combined <- map_dfr(names(pip_config), function(variant_type) {
  
  # Get the map for this variant type (e.g., Shared methods)
  current_method_map <- pip_config[[variant_type]]
  
  map_dfr(names(current_method_map), function(method_name) {
    
    # Identify the correct column
    pip_col <- current_method_map[[method_name]]
    
    # OPTIMIZATION: Sort the dataframe ONCE per method/variant combo
    # Select only necessary columns to keep memory usage low
    sorted_data <- res_all %>%
      select(all_of(c(pip_col, binary_ann, cont_ann))) %>%
      arrange(desc(.data[[pip_col]]))
    
    # Loop through the Top N steps
    map_dfr(top_n_steps, function(k) {
      
      # Slice the top k signals
      sig   <- sorted_data %>% slice(1:k)
      n_sig <- k
      n_bg  <- total_n - k
      
      # --- FAST Signal Calculation ---
      sig_sums_bin <- sig %>%
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))
      
      # --- FAST Background Calculation (Subtraction) ---
      # Background Sum = Total Sum - Signal Sum
      bg_sums_bin <- total_sums_bin - sig_sums_bin
      
      # Calculate Frequencies
      sig_freq <- sig_sums_bin / n_sig
      bg_freq  <- bg_sums_bin / n_bg
      
      # Enrichment Ratio
      enrich_bin <- sig_freq / bg_freq
      
      # --- Continuous Annotations ---
      # Just mean() the signal (no background math needed for this part)
      enrich_cont <- sig %>%
        summarise(across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE)))
      
      # --- Formatting & Combining ---
      bind_rows(
        enrich_bin %>% pivot_longer(everything(), names_to="annotation", values_to="enrichment"),
        enrich_cont %>% pivot_longer(everything(), names_to="annotation", values_to="enrichment")
      ) %>%
        mutate(
          method       = method_name,
          variant      = variant_type,
          top_n        = k
        )
      
    }) # End Top N Loop
  }) # End Method Loop
}) # End Variant Type Loop
# 3. Factor Recoding and Formatting --------------------------------------------

enrichment_final <- enrichment_results_combined %>%
  # 1. Recode Annotation Names
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common", "synonymous", "promotor", "utr_comb", 
                                        "CRE", "heart_ind_eQTL", "artery_ind_eQTL", 
                                        "brain_ind_eQTL", "Z_eur"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common"        = "Non-synonymous",
                             "synonymous"      = "Synonymous",
                             "promotor"        = "Promotor",
                             "utr_comb"        = "UTR",
                             "CRE"             = "CRE",
                             "heart_ind_eQTL"  = "Heart_ind_eQTL",
                             "artery_ind_eQTL" = "Artery_ind_eQTL",
                             "brain_ind_eQTL"  = "Brain_ind_eQTL",
                             "Z_eur"           = "Z_eur")) %>%
  
  # 2. Set Variant Factors
  mutate(variant = factor(variant, levels = c("Shared", "Unique"))) %>%
  
  # 3. Set Method Factors and Labels
  # Includes only the methods calculated above (MFD, MESuSiE, SuSiE, CARMAX)
  mutate(method = factor(method, 
                         levels = c("MFD", "MESuSiE", "SuSiE"),
                         labels = c("MFD", "MESuSiE", "SuSiE")))

# Preview
head(enrichment_final)



### Define x-axis labeling for selected thresholds only: 0.5, 0.6, 0.7, 0.8, and 0.9.
#selected_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
#top_signals_shared <- enrichment_results_combined %>%
#  filter(threshold %in% selected_thresholds) %>%
#  group_by(variant,threshold) %>%
#  summarise(top_signal = unique(top_N_signal)) %>%
#  arrange(threshold) %>% filter(variant == "Shared")

#top_signals_unique <- enrichment_results_combined %>%
#  filter(threshold %in% selected_thresholds) %>%
#  group_by(variant,threshold) %>%
#  summarise(top_signal = unique(top_N_signal)) %>%
#  arrange(threshold) %>% filter(variant == "Unique")

#x_breaks <- selected_thresholds
#x_labels <- paste0(x_breaks, " \n(n=", top_signals_shared$top_signal, ")","\n(n=", top_signals_unique$top_signal, ")")

### Split the data into two panels:
panelA <- enrichment_final %>% 
  filter(!(annotation %in% c("Z_eur"))) %>% dplyr::rename(PIP = variant)

panelB <- enrichment_results_combined %>% 
  filter(annotation %in% c("CADD", "brain_ind_eQTL")) %>% dplyr::rename(PIP = variant)

### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations
plotA <- ggplot(panelA, aes(x = top_n, y = enrichment, color = method, linetype = PIP)) +
  geom_line() +
  geom_point(size = 1.5)+
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel A: Enrichment Results (Signal/Background)",
       x = "Top number based on PIP",
       y = "Enrichment (Signal / Background)",
       color    = "Method",   # legend title for color
       linetype = "Variant"   # legend title for linetype
  ) +
  scale_color_manual(name   = "Method",values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
  )) +
  scale_linetype_manual(name   = "Variant",values = c("Shared" = "solid", "Unique" = "dashed")) +
  theme_bw()+custom_theme()+  theme(plot.title = element_blank(),legend.position  = "bottom")


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline//Figure/"
ggsave(paste0(plot_dir, "Supp_MFD_shared_unique_pip_5_to_9.pdf"),
      plotA, width = 210, height = 180, dpi = 600, units = 'mm')


################################################################
library(openxlsx)
library(stringr) # For str_replace



# Replace .RData with .xlsx for the excel output
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
system(paste0("mkdir -p ",res_dir_table))
res_out_xlsx <-  "Real_data_enrichment_shared_unique.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

sheet_name_raw <- str_remove(res_out_xlsx, ".xlsx")
sheet_name <- substr(sheet_name_raw, 1, 31) 


panelA_to_supp = panelA %>% select(method,annotation,PIP, top_n,enrichment)

wb <- createWorkbook()
addWorksheet(wb, sheet_name)

# Create a bold style for the Title
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# Write the Main Title at Row 1
writeData(wb, sheet_name, "Supplementary Table: Functional enrichment for priotized Shared and Unique Signal", startRow = 1)
addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)

# Write the Combined Data Table starting at Row 2
writeData(wb, sheet_name, panelA_to_supp, startRow = 2)

# 5. Save
saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined file to:", full_path))
print(paste("Saved combined table to:", full_path))





##########################################################


####Combined Shared vs. Unique plot at PIP 0.5 and 0.9


##########################################################




### Split the data into two panels:
panelA <- enrichment_final %>% 
  filter(!(annotation %in% c("Z_eur"))) %>% dplyr::rename(PIP = variant) %>% filter(top_n == 50)

panelA_upper <- enrichment_final %>% 
  filter(!(annotation %in% c("Z_eur", "GTEx_eQTL_MaxCPP"))) %>% dplyr::rename(PIP = variant) %>% filter(top_n == 100)


### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.5
plotA <- ggplot(panelA, 
                aes(x = top_n, 
                    y = enrichment,
                    linetype = PIP,
                    fill = method 
                )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at top 50 Signals",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
    
  )) +
  scale_linetype_manual(values = c(
    "Shared" = "solid",
    "Unique" = "dashed"
  )) +
  theme_bw() +
  custom_theme() +
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
  )


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir, "Enrichment_pip_Shared_Unique_top50_combined.pdf"),
       plotA, width = 210, height = 180, dpi = 600, units = 'mm')

ggsave(paste0(plot_dir, "Enrichment_pip_Shared_Unique_top50_combined.pdf"),
       plotA, width = 10.7, height = 5.92, dpi = 600, units = 'in')


### Create Panel B plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.9
plotA_upper <- ggplot(panelA_upper, 
                      aes(x = top_n, 
                          y = enrichment,
                          linetype = PIP,
                          fill = method 
                      )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at top 100 Signals",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
    
  )) +
  scale_linetype_manual(values = c(
    "Shared" = "solid",
    "Unique" = "dashed"
  )) +
  theme_bw() +
  custom_theme() +
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
  )


### Combine the two panels vertically using patchwork
combined_plot <- (plotA / plotA_upper)+plot_annotation(tag_levels = "a") +  custom_theme()+plot_layout(guides = 'collect') &
  theme(
    axis.title.x  = element_blank(),  # remove x-axis title
    axis.text.x   = element_blank(),  # remove x-axis tick labels
    axis.ticks.x  = element_blank(),   # remove x-axis ticks
    legend.position  = "bottom"
  )


# Display the combined plot
combined_plot


# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"

ggsave(paste0(plot_dir, "combined_enrichment_plot_Shared_Unique_figure5.pdf"),
       combined_plot, width = 7, height = 6.5, dpi = 500, units = 'in')









########################Classfication summary for MFD##################


bmi_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_BMI/result/Finemapping_Method_Decisions_BMI_20251208_162748.csv")
dbp_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_DBP/result/Finemapping_Method_Decisions_DBP_20251208_162748.csv")
sbp_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_SBP/result/Finemapping_Method_Decisions_SBP_20251208.csv")


# 1. Create the reference table used for classification
reference_table <- data.table(
  Method = c(
    "MESuSiE",
    "SuSiE post-hoc",
    "SuSiE post-hoc",
    "MESuSiE",
    "SuSiE post-hoc"
  ),
  Reason = c(
    "AS-V not significant",
    "AS-V significant and in high LD with a Shared SNP significant in ONLY specific ancestry",
    "AS-V significant, but not in high LD with shared signals",
    "AS-V significant and in high LD with a Shared SNP significant in ALL ancestries",
    "AS-V significant, but no significant shared SNPs found"
  )
)

# 2. Perform the Merge (Lookup)
# Map 'V2' from the data to 'Reason' in the reference table.
# all.x = TRUE retains all rows from the original data, even if a match fails.

dbp_result <- merge(dbp_class, reference_table, by.x = "V2", by.y = "Reason", all.x = TRUE)
sbp_result <- merge(sbp_class, reference_table, by.x = "V2", by.y = "Reason", all.x = TRUE)

# 3. (Optional) Restore original sort order based on ID (V1)
setorder(dbp_result, V1)
setorder(sbp_result, V1)

setnames(dbp_result, old = c("V1", "V2"), new = c("Locus_ID", "Reason"))
setnames(sbp_result, old = c("V1", "V2"), new = c("Locus_ID", "Reason"))

# 2. Add a 'Trait' column to distinguish the data sources
bmi_class[, Trait := "BMI"]
dbp_result[, Trait := "DBP"]
sbp_result[, Trait := "SBP"]

# 3. Merge (Row Bind) the datasets
# Select only the common columns (Locus_ID, Method, Reason, Trait)
cols_to_keep <- c("Locus_ID", "Method", "Reason", "Trait")
all_data <- rbind(
  bmi_class[, ..cols_to_keep], 
  dbp_result[, ..cols_to_keep], 
  sbp_result[, ..cols_to_keep]
)


all_data$Locus_ID <- seq_len(nrow(all_data))


library(openxlsx)

res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
res_out_xlsx <- "Real_data_MFD_summary.xlsx"
full_path_eqtl <- paste0(res_dir_table, res_out_xlsx)

write.xlsx(all_data, file = full_path_eqtl, overwrite = TRUE)


summary_stats <- all_data %>%
  count(Trait, Method) %>%          # Get counts (creates column 'n')
  group_by(Trait) %>%               # Group by Trait to calculate sum within each Trait
  mutate(Proportion = n / sum(n))   # Calculate proportion

print(summary_stats)

# 1. Summary by Trait and Method
summary_stats <- all_data %>%
  count(Trait, Method,Reason) %>%          # Get counts (creates column 'n')
  group_by(Trait) %>%               # Group by Trait to calculate sum within each Trait
  mutate(Proportion = n / sum(n))   # Calculate proportion

print(summary_stats)

overall_stats <- all_data %>%
  count(Method) %>%                 # Get counts by Method
  mutate(Proportion = n / sum(n))   # Calculate proportion of total

print(overall_stats)
