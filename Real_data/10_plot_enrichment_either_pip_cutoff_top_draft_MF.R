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
library(purrr)
library(ggpattern)

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summary_res/"
load(paste0(res_dir,"res_all_mf_updated.RData"))   

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region+ 204, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region+ 243, res_all$Region)
res_all$CS <- ifelse(res_all$CS != 0, 1, res_all$CS)
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))




ann_col_name<-c( "non_synonymous_common","synonymous", "promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")




custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),  
    axis.title.x = element_text(size = 10, face="bold"),
    axis.title.y = element_text(size = 10, face="bold"),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    strip.background = element_blank(),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12, face="bold"),
    plot.title = element_text(size=11, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}


# Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
binary_ann <- ann_col_name

cont_ann <- c("Z_eur")

# Select columns containing "PIP_Either" and apply the sum condition
col_sums <- res_all %>%
  select(contains("PIP_Either")) %>%
  summarise(across(everything(), ~ sum(. > 0.5, na.rm = TRUE)))

# Print the result (transposed for easier reading)
print(t(col_sums))

# 1) Define the fixed numbers of top signals to retrieve
top_n_steps <- seq(10, 200, by = 10)


# 2) PIP columns
pip_cols <- c(
  MFD          = "PIP_Either",
  MESuSiE      = "MESuSiE_PIP_Either",
  SuSiE        = "SuSiE_PIP_Either",
  Paintor      = "Paintor_PIP_Either",
  SuSiEx       = "SuSiEx_PIP_Either",
  XMAP         = "XMAP_PIP_Either",
  CARMAX       = "CARMAX_PIP_Either",
  SuSiE_merged = "SuSiE_merged_PIP_Either"
)

total_n <- nrow(res_all)
total_sums_bin <- res_all %>%
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))

# 2. Loop by METHOD first (Outer Loop) to minimize sorting
enrichment_results <- map_dfr(names(pip_cols), function(method) {
  
  pip_col <- pip_cols[[method]]
  
  # OPTIMIZATION: Sort the ENTIRE dataframe just ONCE per method
  # Select only necessary columns to save memory
  sorted_data <- res_all %>%
    select(all_of(c(pip_col, binary_ann, cont_ann))) %>%
    arrange(desc(.data[[pip_col]]))
  
  # 3. Loop through K steps (Inner Loop)
  # Slice only the small top_k subset without loading the full background
  map_dfr(top_n_steps, function(k) {
    
    # Slice the top k signals
    sig <- sorted_data %>% slice(1:k)
    n_sig <- k
    n_bg  <- total_n - k
    
    # --- FAST Signal Calculation ---
    # Sum of binary annotations in Top K
    sig_sums_bin <- sig %>%
      summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))
    
    # --- FAST Background Calculation (Subtraction Method) ---
    # Background Sum = Total Sum - Signal Sum
    bg_sums_bin <- total_sums_bin - sig_sums_bin
    
    # Calculate Frequencies (Mean = Sum / N)
    sig_freq <- sig_sums_bin / n_sig
    bg_freq  <- bg_sums_bin / n_bg
    
    # Enrichment
    enrich_bin <- sig_freq / bg_freq
    
    # --- Continuous Annotations ---
    # Calculate the signal mean without loading the background
    enrich_cont <- sig %>%
      summarise(across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE)))
    
    # --- Formatting ---
    bin_long  <- enrich_bin %>%
      pivot_longer(everything(), names_to="annotation", values_to="enrichment")
    
    cont_long <- enrich_cont %>%
      pivot_longer(everything(), names_to="annotation", values_to="enrichment")
    
    bind_rows(bin_long, cont_long) %>%
      mutate(
        method = method,
        top_n  = k
      )
  })
})

# Reorder the 'annotation' factor in enrichment_results
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, 
                             levels = c("non_synonymous_common",
                                        "synonymous", 
                                        "promotor","utr_comb",
                                        "CRE",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL","brain_ind_eQTL","Z_eur"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "CRE" = "CRE",
                             "heart_ind_eQTL" = "Heart ind. eQTL",
                             "artery_ind_eQTL" = "Artery ind. eQTL",
                             "brain_ind_eQTL" = "Brain ind. eQTL",
                             "Z_eur" = "Z_eur"))

enrichment_results$method <- factor(
  enrichment_results$method, 
  # 1. 'levels' must match the CURRENT strings in the data frame
  levels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged","Paintor", "SuSiEx", "XMAP", "CARMAX"),
  # 2. 'labels' are the names assigned to the levels
  labels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta","Paintor", "SuSiEx", "XMAP", "CARMAx")
)

# Create the plot with colorblind-friendly colors
enrichment_plot <- ggplot(enrichment_results%>% 
                            filter(!(annotation %in% c("Z_eur"))), aes(x = top_n, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Enrichment Results Across  Cutoffs",
    x = "Top number based on PIP",
    y = "Enrichment (Signal / Background)"
  ) +
  # Manually set colors for each method
  scale_color_manual(values = 
                       
                       c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")
                     
  ) +
  theme_bw()+  custom_theme()+theme(plot.title = element_blank(),legend.position  = "bottom")

enrichment_plot
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir, "Supp_MFD_Either_pip_5_to_9.pdf"),
       enrichment_plot, width = 210, height = 180, dpi = 600, units = 'mm')





######################################


#Barplot at PIP threshold 0.5


######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

fill_colors =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")

# Create Panel A plot: Enrichment (Signal / Background)
enrichment_panelA = enrichment_results %>% filter(top_n==100) %>% filter(annotation != "Z_eur")
# Create the base plot with the 2x4 layout
plotA <- ggplot(enrichment_panelA, aes(x = method, y = enrichment, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  
  # 1. Facet Layout: 2 rows, 4 columns
  facet_wrap(~ annotation, scales = "free_y", nrow = 2, ncol = 4) +
  
  labs(
    title = "Panel A: Enrichment Results",
    x = "Method",
    y = "Enrichment (Signal / Background)",
    fill = "Method"
  ) +
  
  scale_fill_manual(name = "Method", values = fill_colors) +
  
  # 2. Force Legend to 1 Row
  guides(fill = guide_legend(nrow = 1)) +
  
  theme_bw() + 
  custom_theme() + 
  theme(
    legend.position = "bottom",    
    plot.title = element_blank(),
    
    # 3. Remove X-axis text and ticks (cleaner for small multiples)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Add annotations
plotA_annot <- plotA + 
  geom_text(
    aes(label = round(enrichment, 2)), 
    position = position_dodge(width = 1), 
    vjust = -0.5, 
    size = 5 * 5 / 10
  )

print(plotA_annot)

#####Data for summary table 1#####
processed_enrichment <- enrichment_panelA %>%
  # Create a grouping column
  mutate(group = case_when(
    annotation %in% c("Non-synonymous", "Synonymous") ~ "Coding",
    annotation %in% c("Promotor", "UTR", "CRE") ~ "Regulatory",
    str_detect(annotation, "eQTL") ~ "eQTL"
  )) %>%
  # Group by Method + Category and Calculate Average
  group_by(method, group) %>%
  summarise(avg_enrichment = mean(enrichment), .groups = "drop") %>%
  # Reshape to Wide Format (Columns for each category)
  pivot_wider(
    names_from = group, 
    values_from = avg_enrichment
  ) %>%
  # Rename columns to match the FunkyHeatmap code used below
  rename(
    val_coding = Coding,
    val_regulatory = Regulatory, # Regulatory annotation
    val_eqtl = eQTL,
    id_map = method # Join key for method summaries
  )


#####Data for rank based evaluation#####

# --- 1. Setup Paths and Levels ---
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
system(paste0("mkdir -p ", res_dir_table))
res_out_xlsx <- "Real_data_enrichment_either_rank.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

# Define the order from the factor levels
desired_levels <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP", "CARMAx")

# --- 2. Process Enrichment Data ---
enrichment_results_to_supp <- enrichment_results %>% 
  filter(!(annotation %in% c("Z_eur"))) %>% 
  select(method, annotation, top_n, enrichment) %>% 
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(method)

# --- 3. Process Rank Based Evaluation ---

# A. Summarize Data
method_summary <- enrichment_panelA %>%
  filter(!grepl("eQTL", annotation), annotation != "Z_eur") %>%
  group_by(method, annotation) %>%
  summarise(
    mean_value = mean(enrichment, na.rm = TRUE), 
    .groups = "drop"
  )

# B. Calculate Rank by Annotation (Long Format)
# Apply desired_levels here
ranked_by_annotation <- method_summary %>%
  group_by(annotation) %>%
  mutate(
    rank = min_rank(desc(mean_value)) 
  ) %>%
  ungroup() %>%
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(annotation, rank)

# C. Create Pivot Table for Ranks (Wide Format)
# Apply desired_levels here as well so the rows are sorted by method order
rank_table <- ranked_by_annotation %>% 
  select(annotation, rank, method) %>% 
  pivot_wider(names_from = annotation, values_from = rank) %>%
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(method)

# D. Overall Average Rank (Optional Check)
overall_method_rank <- ranked_by_annotation %>%
  group_by(method) %>%
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    med_rank = median(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(avg_rank)

# --- 4. Write to Excel (Multiple Sheets) ---

wb <- createWorkbook()
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# Sheet 1: Enrichment Results
addWorksheet(wb, "Enrichment_Results")
writeData(wb, "Enrichment_Results", "Supplementary Table: Functional enrichment for prioritized Either Signal with Rank", startRow = 1)
addStyle(wb, "Enrichment_Results", title_style, rows = 1, cols = 1)
writeData(wb, "Enrichment_Results", enrichment_results_to_supp, startRow = 2)

# Sheet 2: Ranked Long Format
addWorksheet(wb, "Ranked_Annotation_Long")
writeData(wb, "Ranked_Annotation_Long", "Supplementary Table: Method Ranks by Annotation (Long Format)", startRow = 1)
addStyle(wb, "Ranked_Annotation_Long", title_style, rows = 1, cols = 1)
writeData(wb, "Ranked_Annotation_Long", ranked_by_annotation, startRow = 2)

# Sheet 3: Ranked Matrix (Wide)
addWorksheet(wb, "Ranked_Matrix_Wide")
writeData(wb, "Ranked_Matrix_Wide", "Supplementary Table: Method Ranks Matrix", startRow = 1)
addStyle(wb, "Ranked_Matrix_Wide", title_style, rows = 1, cols = 1)
writeData(wb, "Ranked_Matrix_Wide", rank_table, startRow = 2)

# --- 5. Save ---
saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined file to:", full_path))
print("Results Check:")
print("Overall Average Rank:")
print(overall_method_rank)



################################################

# Number of CSs by traits
region_cs<-res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE)))

cs_number_by_trait = region_cs%>%group_by(PHENONAME) %>% summarise(sum(CS!=0),sum(MESuSiE_cs!=0),sum(SuSiE_cs!=0),sum(Paintor_cs!=0),sum(SuSiEx_cs!=0),sum(XMAP_cs !=0),sum(CARMAX_cs !=0),sum(SuSiE_metal_cs !=0))

data_a <- cs_number_by_trait %>%
  # Rename grouping column
  rename(Trait = PHENONAME) %>%
  
  # Pivot all columns except Trait
  pivot_longer(
    cols = -Trait,
    names_to = "raw_col",
    values_to = "CS_Count"
  ) %>%
  
  # Clean up the method names
  mutate(Method = case_when(
    # 1. Map base CS column to MFD
    raw_col == "sum(CS != 0)" ~ "MFD",
    
    # Standardize method names
    grepl("SuSiE_metal", raw_col) ~ "SuSiE_Meta",
    grepl("CARMAX", raw_col)      ~ "CARMAx",
    
    # 3. Default rule: clean the wrapper text for everything else
    TRUE ~ str_remove_all(raw_col, "^sum\\(|_cs != 0\\)$")
  )) %>%
  
  # Reorder columns
  select(Method, Trait, CS_Count)





# 1. Force the stacking order: BMI (left) -> DBP -> SBP (right)
data_a$Trait <- factor(data_a$Trait, levels = c("SBP", "DBP", "BMI"))

# 2. Define the colors and order (as previously set)
# Ensure Method levels are reversed so MFD is at the top of the Y-axis
data_a$Method <- factor(data_a$Method, levels = rev(c("MFD", "MESuSiE", "SuSiE","SuSiE_Meta","Paintor", "SuSiEx", 
                                                      "MultiSuSiE", "XMAP", "CARMAx")))

# 3. Create the plot

pattern_map <- c("BMI" = "none", "DBP" = "stripe", "SBP" = "circle")


# 4. Create the Plot
p_a <- ggplot(data_a, aes(x = CS_Count, y = Method, fill = Method)) +
  geom_col_pattern(
    aes(pattern = Trait),
    position = "stack",
    
    # --- Aesthetic Optimizations ---
    color = "white",                # White borders
    size = 0.2,                     # Thin border
    
    # Pattern Styling
    pattern_colour = "white",       # Pattern lines are white
    pattern_fill   = "white",       # Pattern fill is white
    pattern_alpha  = 0.6,           # Transparency
    pattern_density = 0.2,          # Density
    pattern_spacing = 0.02,         # Spacing
    pattern_angle = 45,             # Angle
    pattern_key_scale_factor = 0.5  # Legend key scale
  ) +
  
  # Method Colors (Hidden Legend)
  scale_fill_manual(values = fill_colors, guide = "none") + 
  
  # Pattern mapping
  scale_pattern_manual(
    values = pattern_map,
    name = "Trait",
    
    # Set the legend order to BMI -> DBP -> SBP
    # This changes the legend display without changing the plot stacking order.
    breaks = c("BMI", "DBP", "SBP"), 
    
    guide = guide_legend(override.aes = list(fill = "grey70")) 
  ) +
  
  # Labels
  labs(x = "Number of credible sets", y = NULL) +
  
  # Count Labels
  geom_text(aes(label = after_stat(x), group = Method), 
            stat = 'summary', fun = sum, hjust = -0.2, size = 2.8, fontface = "plain") +
  
  # Scales & Spacing
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Theme
  custom_theme() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Legend Position (inside plot, bottom right or custom coords)
    legend.position = c(0.95, 0.2), 
    legend.background = element_rect(fill = "white", color = NA)
  )

print(p_a)

################################################




################################################

# Size and Z-score plot


################################################
################################################
#
#		Set SiZe Part
#
###############################################		
###Credible set size by PHENONAME


all_sets_info<-data.frame(res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE))))#%>%filter(CS!= 0,MESuSiE_cs!=0,SuSiE_cs!=0,Paintor_cs!=0,SuSiEx_cs !=0,XMAP_cs!=0,CARMAX_cs!=0) #Median Set Size across all locus

# 1. Prepare the data (assuming all_sets_info is already loaded)
all_sets_info_long <- all_sets_info %>%
  # Select only the required columns (excluding the indices if they exist as row names)
  select(PHENONAME, Region, CS, MESuSiE_cs, SuSiE_cs, Paintor_cs, 
         SuSiEx_cs, XMAP_cs, CARMAX_cs, SuSiE_metal_cs) %>%
  # Pivot to long format
  pivot_longer(
    cols = c("CS", "MESuSiE_cs", "SuSiE_cs", "Paintor_cs", 
             "SuSiEx_cs", "XMAP_cs", "CARMAX_cs", "SuSiE_metal_cs"),
    names_to = "Method", 
    values_to = "Count"
  )

# 2. Filter Zeros, Re-level, and Log-Transform
plot_data <- all_sets_info_long %>%
  # CRITICAL STEP: Remove 0s before plotting/transforming
  filter(Count > 0) %>% 
  # Log transform (using Count instead of Count+1 is usually cleaner if 0s are gone, 
  # Retain the log2 transformation used below.
  mutate(LogCount = log(Count+1, base = 2)) %>% 
  # Re-factor and Rename levels
  mutate(Method = factor(Method, 
                         levels = c("CS", "MESuSiE_cs", "SuSiE_cs", "SuSiE_metal_cs","Paintor_cs", 
                                    "SuSiEx_cs", "XMAP_cs", "CARMAX_cs"),
                         labels = c("MFD", "MESuSiE", "SuSiE","SuSiE_Meta", "Paintor", 
                                    "SuSiEx", "XMAP", "CARMAx")))


# all_sets_info_long_save <- all_sets_info_long %>%
#   mutate(causal = num_causal,
#          h2 = h2_num)


p_set = ggplot(data =plot_data,aes(x = PHENONAME, y=LogCount,fill=Method))+geom_boxplot(aes(x = PHENONAME,fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray")+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"),guide=FALSE)
p_set =p_set + theme_bw() + xlab("") +ylab("log2(Set Size + 1)")+coord_cartesian(ylim=c(0,10))
p_set= p_set+custom_theme()

median_stats <- plot_data %>%
  group_by(PHENONAME, Method) %>%
  summarise(
    Median_CS_Size = median(Count, na.rm = TRUE), 
    .groups = 'drop'
  ) %>%
  arrange(PHENONAME, Median_CS_Size)

# Display the result
print(median_stats)

median_stats %>% group_by(Method) %>% summarise(mean(Median_CS_Size))
################################################
#
#		Z-score Part
#
###############################################		
MFD_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CS==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))
MESuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)
Paintor_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiEx_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiEx_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
XMAP_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(XMAP_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
CARMAX_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CARMAX_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_metal_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_metal_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)

set_size_z_info<-data.frame(cbind(MFD_cs_Z,MESuSiE_cs_Z,SuSiE_cs_Z,Paintor_cs_Z,SuSiEx_cs_Z,XMAP_cs_Z,CARMAX_cs_Z,SuSiE_metal_cs_Z))
colnames(set_size_z_info)<-c("PHENONAME",c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta"))
set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(PHENONAME), names_to = "Method", values_to = "Z")%>%mutate(Method = factor(Method, levels=c("MFD","MESuSiE","SuSiE", "SuSiE_Meta","Paintor","SuSiEx","XMAP","CARMAx")))

set_size_z_info_long %>% group_by(Method) %>% summarise(mean(Z))
#p_z = ggplot(data = set_size_z_info_long,aes(x = PHENONAME, y=Z,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"))
#p_z = p_z + theme_bw() + xlab("") +ylab("Median |Z|")+ ylim(0,max(round(set_size_z_info_long$Z,2)+1)) +custom_theme()
#p_z = p_z + geom_text(label = round(set_size_z_info_long$Z,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)

p_z = ggplot(data = set_size_z_info_long, aes(x = PHENONAME, y = Z, fill = Method)) +
  # 1. Be explicit about the width in geom_bar to match it later
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values = c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8",
                               "SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3",
                               "SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"))

# 2. Move label INSIDE aes() and match the position_dodge width
p_z = p_z + geom_text(
  aes(label = round(Z, 2)),  # Display rounded Z values
  position = position_dodge(width = 0.9), # Match the geom_bar width
  vjust = -0.5,
  size = 5 * 5 /10
)

p_z = p_z + theme_bw() + xlab("") + ylab("Median |Z|") + 
  ylim(0, max(round(set_size_z_info_long$Z, 2) + 1))
p_z = p_z + custom_theme()


set_size_z_info_long %>% group_by( Method) %>%   summarise(
  Mean_Z = mean(Z, na.rm = TRUE), 
  .groups = 'drop'
)

####################################################################################
#
#
#			Combine the first real data plots
#
#
#
####################################################################################


p_a = p_a
p_b = p_set
p_c = p_z
p_d = plotA_annot


p_a_ready <- p_a 

# p_b, p_c: Remove legends entirely
p_b_ready <- p_b + guides(fill = "none", color = "none", pattern = "none")
p_c_ready <- p_c + guides(fill = "none", color = "none")

# p_d: Keep legend (it will be collected for the bottom row)
p_d_ready <- p_d + theme(legend.position = "bottom")

# 2. Define Rows Separately (The "Nesting" Strategy)
# Top Row: Just p_a and p_b. No collecting here, so p_a keeps its internal legend.
row1 <- (p_a_ready | p_b_ready)

# Bottom Row: p_c and p_d. Collect legends HERE.
# This takes p_d's legend and centers it at the bottom of this row.
row2 <- (p_c_ready | p_d_ready) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.box.margin = margin(t = 10))

# 3. Combine Rows
combined_plot <- (row1 / row2) +
  plot_annotation(
    #title = "Multi-ancestry fine-mapping on All of Us data",
    tag_levels = 'a',
    theme = theme(
     # plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = "sans"),
      plot.tag = element_text(size = 12, face = "bold", family = "sans")
    )
  ) &
  # Global font styling
  theme(text = element_text(family = "sans"))

# 4. Print/Save
print(combined_plot)


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir, "set_size_enrichment_pip_top100_combined_updated_test.pdf"),
       combined_plot, width = 550, height = 300, dpi = 600, units = 'mm')

ggsave(paste0(plot_dir, "set_size_enrichment_pip_top100_combined_updated_test.pdf"),
       combined_plot, width = 210, height = 180, dpi = 600, units = 'mm')






####################################################################################
#
#
#			Enrichment of annotations group by traits
#
#
#
####################################################################################

library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(stringr)

# --- 1. Setup and Data Loading ---
res_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summary_res/"
load(paste0(res_dir, "res_all_mf_updated.RData"))

# Pre-processing
res_all <- res_all %>%
  mutate(
    n_regions = n_distinct(Region), # Calculate distinct regions if needed later
    Region = case_when(
      PHENONAME == "DBP" ~ Region + 204,
      PHENONAME == "SBP" ~ Region + 243,
      TRUE ~ Region
    ),
    CS = ifelse(CS != 0, 1, CS),
    utr_comb = ifelse((utr_3 + utr_5) > 0, 1, 0)
  )

# Define Annotations
ann_col_name <- c("non_synonymous_common", "synonymous", "promotor", 
                  "utr_comb", "CRE", "heart_ind_eQTL", 
                  "artery_ind_eQTL", "brain_ind_eQTL")
binary_ann <- ann_col_name
cont_ann <- c("Z_eur")

# Define Methods Map
pip_cols <- c(
  MFD          = "PIP_Either",
  MESuSiE      = "MESuSiE_PIP_Either",
  SuSiE        = "SuSiE_PIP_Either",
  Paintor      = "Paintor_PIP_Either",
  SuSiEx       = "SuSiEx_PIP_Either",
  XMAP         = "XMAP_PIP_Either",
  CARMAX       = "CARMAX_PIP_Either",
  SuSiE_merged = "SuSiE_merged_PIP_Either"
)

# --- 2. Pre-Calculate Totals per Phenotype ---
# Calculate the total sum of annotations and total N per phenotype *once*
# to serve as the "Universe" for background calculations.
pheno_stats <- res_all %>%
  group_by(PHENONAME) %>%
  summarise(
    total_n = n(),
    across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE), .names = "total_{.col}"),
    .groups = "drop"
  )

# --- 3. Main Enrichment Loop ---
target_k <- 100 # Use the prespecified value

enrichment_results <- map_dfr(names(pip_cols), function(method) {
  
  pip_col <- pip_cols[[method]]
  
  # Select relevant data
  method_data <- res_all %>%
    select(PHENONAME, all_of(c(pip_col, binary_ann, cont_ann))) %>%
    rename(PIP = !!sym(pip_col))
  
  # Group by Phenotype, Sort by PIP, Slice Top 100
  top_hits <- method_data %>%
    group_by(PHENONAME) %>%
    arrange(desc(PIP)) %>%
    slice_head(n = target_k) %>% 
    ungroup()
  
  # Summarize the Signals (Top 100) per Phenotype
  sig_stats <- top_hits %>%
    group_by(PHENONAME) %>%
    summarise(
      across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE), .names = "sig_{.col}"),
      across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )
  
  # Merge Signal stats with Total stats to calculate Background and Enrichment
  merged <- inner_join(sig_stats, pheno_stats, by = "PHENONAME")
  
  # Calculate Enrichment
  results <- merged %>%
    mutate(
      method = method,
      top_n = target_k
    )
  
  # Loop through annotations to calculate enrichment ratio for this method/pheno
  for (ann in binary_ann) {
    sig_col <- paste0("sig_", ann)
    tot_col <- paste0("total_", ann)
    
    # Define variables for calculation
    n_sig <- target_k
    # Background N = Total rows for that pheno - 100
    n_bg <- results$total_n - target_k 
    
    # Signal Sum
    s_sum <- results[[sig_col]]
    # Background Sum = Total Sum - Signal Sum
    b_sum <- results[[tot_col]] - s_sum
    
    # Frequencies
    sig_freq <- s_sum / n_sig
    bg_freq <- b_sum / n_bg
    
    # Calculate Ratio (Handle division by zero if bg_freq is 0)
    enrichment_val <- ifelse(bg_freq > 0, sig_freq / bg_freq, NA)
    
    # Assign back to a clean column name
    results[[ann]] <- enrichment_val
  }
  
  # Clean up and pivot long
  long_res <- results %>%
    select(PHENONAME, method, top_n, all_of(binary_ann), starts_with("mean_")) %>%
    rename(Z_eur = mean_Z_eur) %>% # Rename continuous back
    pivot_longer(
      cols = c(all_of(binary_ann), "Z_eur"),
      names_to = "annotation",
      values_to = "enrichment"
    )
  
  return(long_res)
})

# --- 4. Formatting for Plotting ---
enrichment_results <- enrichment_results %>%
  mutate(
    annotation = factor(annotation, 
                        levels = c("non_synonymous_common", "synonymous", "promotor", "utr_comb", 
                                   "CRE", "heart_ind_eQTL", "artery_ind_eQTL", "brain_ind_eQTL", "Z_eur")),
    annotation = recode(annotation,
                        "non_synonymous_common" = "Non-synonymous",
                        "synonymous" = "Synonymous",
                        "promotor" = "Promotor",
                        "utr_comb" = "UTR",
                        "CRE" = "CRE",
                        "heart_ind_eQTL" = "Heart ind. eQTL",
                        "artery_ind_eQTL" = "Artery ind. eQTL",
                        "brain_ind_eQTL" = "Brain ind. eQTL",
                        "Z_eur" = "Z_eur"),
    method = factor(method, 
                    levels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged", "Paintor", "SuSiEx", "XMAP", "CARMAX"),
                    labels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP", "CARMAx"))
  )


# Define Colors
method_colors <- c(
  "MESuSiE" = "#8da0cb", "SuSiE" = "#66c2a5", "MFD" = "#377eb8", 
  "SuSiE_Meta" = "#9FBA95", "Paintor" = "#fc8d62", "MultiSuSiE" = "#e78ac3", 
  "SuSiEx" = "#E89DA0", "XMAP" = "#ffd92f", "CARMAx" = "#f2b56e"
)

enrichment_plot <- ggplot(enrichment_results %>% filter(!(annotation %in% c("Z_eur"))), aes(x = PHENONAME, y = enrichment, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ annotation, scales = "free_y", ncol = 3) +
  labs(
    title = paste0("Enrichment Results (Top ", target_k, " signals) by Phenotype"),
    x = "Phenotype",
    y = "Enrichment (Signal / Background)"
  ) +
  scale_fill_manual(values = method_colors) +
  theme_bw() +  custom_theme()+theme(plot.title = element_blank(),legend.position  = "bottom")

print(enrichment_plot)

plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir, "Supp_MFD_Either_pip_top_100_by_trait.pdf"),
       enrichment_plot, width = 210, height = 180, dpi = 600, units = 'mm')



#####Data for rank based evaluation#####

# --- 1. Setup Paths and Levels ---
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
# No need to mkdir again if it exists, but safe to keep
system(paste0("mkdir -p ", res_dir_table))

res_out_xlsx_eqtl <- "Real_data_enrichment_either_eQTL.xlsx"
full_path_eqtl <- paste0(res_dir_table, res_out_xlsx_eqtl)

# Define the order from the factor levels
desired_levels <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP", "CARMAx")

# --- 2. Process eQTL Enrichment Data ---

# Restrict to eQTL-related annotations
enrichment_results_eqtl <- enrichment_results %>% 
  filter(str_detect(annotation, "eQTL"))

# Summarize enrichment by Method x Annotation x Trait
method_summary <- enrichment_results_eqtl %>% 
  group_by(PHENONAME, method, annotation) %>% 
  summarise(
    mean_value = mean(enrichment, na.rm = TRUE), 
    .groups = "drop"
  )

# --- 3. Rank Methods within each Annotation x Trait (Long Format) ---
ranked_by_annotation <- method_summary %>% 
  group_by(PHENONAME, annotation) %>% 
  mutate(
    rank = min_rank(desc(mean_value)) 
  ) %>% 
  ungroup() %>%
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(PHENONAME, annotation, rank)

# --- 4. Create Ranking Table (Wide Format) ---
rank_table <- ranked_by_annotation %>% 
  select(PHENONAME, annotation, method, rank) %>% 
  pivot_wider(
    names_from = annotation, 
    values_from = rank
  ) %>%
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(PHENONAME, method)

# --- 5. Compute Average and Median Rank per Method within each Trait ---
overall_method_rank <- ranked_by_annotation %>% 
  # keep only traits of interest
  filter(PHENONAME %in% c("BMI", "DBP", "SBP")) %>% 
  # apply trait-specific annotation inclusion logic
  filter(
    (PHENONAME == "BMI" & str_detect(annotation, "Brain\\s*ind\\.\\s*eQTL")) | 
      (PHENONAME %in% c("DBP", "SBP") & str_detect(annotation, "(Heart|Artery)\\s*ind\\.\\s*eQTL"))
  ) %>% 
  group_by(PHENONAME, method) %>% 
  summarise(
    avg_rank = mean(rank, na.rm = TRUE), 
    med_rank = median(rank, na.rm = TRUE), 
    n_ann    = n_distinct(annotation), 
    .groups = "drop"
  ) %>% 
  mutate(
    method = factor(method, levels = desired_levels)
  ) %>%
  arrange(PHENONAME, avg_rank)


combined_table <- rank_table %>%
  left_join(overall_method_rank, by = c("PHENONAME", "method"))

# --- 6. Write to New Excel File (Multiple Sheets) ---

wb <- createWorkbook()
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# Sheet 1: Ranked Long Format
addWorksheet(wb, "eQTL_Ranked_Long")
writeData(wb, "eQTL_Ranked_Long", "Supplementary Table: eQTL Method Ranks by Annotation (Long Format)", startRow = 1)
addStyle(wb, "eQTL_Ranked_Long", title_style, rows = 1, cols = 1)
writeData(wb, "eQTL_Ranked_Long", ranked_by_annotation, startRow = 2)

# Sheet 2: Ranked Matrix (Wide)
addWorksheet(wb, "eQTL_Ranked_Wide")
writeData(wb, "eQTL_Ranked_Wide", "Supplementary Table: eQTL Method Ranks Matrix", startRow = 1)
addStyle(wb, "eQTL_Ranked_Wide", title_style, rows = 1, cols = 1)
writeData(wb, "eQTL_Ranked_Wide", rank_table, startRow = 2)

# Sheet 3: Overall Summary
addWorksheet(wb, "eQTL_Overall_Summary")
writeData(wb, "eQTL_Overall_Summary", "Supplementary Table: Overall eQTL Rank Summary by Trait", startRow = 1)
addStyle(wb, "eQTL_Overall_Summary", title_style, rows = 1, cols = 1)
writeData(wb, "eQTL_Overall_Summary", overall_method_rank, startRow = 2)

# --- 7. Save ---
saveWorkbook(wb, full_path_eqtl, overwrite = TRUE)

print(paste("Saved eQTL results file to:", full_path_eqtl))
print("Final Check - Overall Ranks:")
print(overall_method_rank)
