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
load(paste0(res_dir,"res_all_mf.RData"))   

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region+ 204, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region+ 243, res_all$Region)
res_all$CS <- ifelse(res_all$CS != 0, 1, res_all$CS)
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))

#########################################################


####Combined Shared vs. Unique plot PIP 0.5 to 0.9


##########################################################

res_all = res_all %>% mutate(MFD_PIP_Unique = pmax(PIP_Ancestry_1,PIP_Ancestry_2),
                             SuSiE_PIP_Unique = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),
                             MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_WB,MESuSiE_PIP_BB),
                             CARMAX_PIP_Unique = pmax(CARMAX_PIP_WB,CARMAX_PIP_BB))


custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
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
ann_col_name<-c( "missense","synonymous", "promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")




thresholds <- seq(0.5, 0.9, by = 0.05)
binary_ann <- ann_col_name
cont_ann   <- c("Z_eur")

# Define the Column Map
# Select the configured column for each method and variant category.
pip_config <- list(
  Shared = c(
    MFD     = "PIP_Shared",           # Shared-variant PIP column
    MESuSiE = "MESuSiE_PIP_Shared",
    SuSiE   = "SuSiE_PIP_Shared",
    CARMAX  = "CARMAX_PIP_Shared"
  ),
  Unique = c(
    MFD     = "MFD_PIP_Unique",       # Ancestry-specific PIP column
    MESuSiE = "MESuSiE_PIP_Unique",
    SuSiE   = "SuSiE_PIP_Unique",
    CARMAX  = "CARMAX_PIP_Unique"
  )
)

# 2. Compute Enrichment --------------------------------------------------------

enrichment_results_combined <- map_dfr(names(pip_config), function(variant_type) {
  
  # Extract the specific method-to-column map for this variant type
  current_method_map <- pip_config[[variant_type]]
  
  map_dfr(thresholds, function(thresh) {
    
    map_dfr(names(current_method_map), function(method_name) {
      
      # 1. Identify the correct column
      pip_col <- current_method_map[[method_name]]
      
      # 2. Filter Signal (n_sig determined by specific column > thresh)
      sig_df <- res_all %>% filter(.data[[pip_col]] > thresh)
      n_sig  <- nrow(sig_df)
      
      # Handle edge case if n_sig is 0 to avoid division by zero errors
      if(n_sig == 0) return(NULL)
      
      # 3. Calculate Background Frequency (Non-Signal)
      # Usage: (Sum of annotation in non-signal rows) / (Total rows - n_sig)
      bg_bin <- res_all %>%
        summarise(across(all_of(binary_ann), 
                         ~ sum(.x, na.rm = TRUE) / (n() - n_sig)))
      
      # 4. Calculate Signal Frequency
      sig_bin <- sig_df %>%
        summarise(across(all_of(binary_ann), 
                         ~ sum(.x, na.rm = TRUE) / n_sig))
      
      # 5. Compute Enrichment Ratio (Signal / Background)
      enrich_bin <- sig_bin / bg_bin
      
      # 6. Compute Continuous Means (e.g., Z_eur, CADD, etc.)
      enrich_cont <- sig_df %>%
        summarise(across(all_of(cont_ann), 
                         ~ mean(.x, na.rm = TRUE)))
      
      # 7. Pivot to Long Format and Combine
      bind_rows(
        enrich_bin %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment"),
        enrich_cont %>% pivot_longer(everything(), names_to = "annotation", values_to = "enrichment")
      ) %>%
        mutate(
          method = method_name,
          threshold = thresh,
          variant = variant_type,
          top_N_signal = n_sig
        )
      
    }) # End Method Loop
  }) # End Threshold Loop
}) # End Variant Type Loop


# 3. Factor Recoding and Formatting --------------------------------------------

enrichment_final <- enrichment_results_combined %>%
  # 1. Recode Annotation Names
  mutate(annotation = factor(annotation, 
                             levels = c("missense", "synonymous", "promotor", "utr_comb", 
                                        "CRE", "heart_ind_eQTL", "artery_ind_eQTL", 
                                        "brain_ind_eQTL", "Z_eur"))) %>%
  mutate(annotation = recode(annotation,
                             "missense"        = "Missense",
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
                         levels = c("MFD", "MESuSiE", "SuSiE", "CARMAX"),
                         labels = c("MFD", "MESuSiE", "SuSiE", "CARMAx")))

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
plotA <- ggplot(panelA, aes(x = threshold, y = enrichment, color = method, linetype = PIP)) +
  geom_line() +
  geom_point(size = 1.5)+
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel A: Enrichment Results (Signal/Background)",
       x = "PIP Threshold",
       y = "Enrichment (Signal / Background)",
       color    = "Method",   # legend title for color
       linetype = "Variant"   # legend title for linetype
  ) +
  scale_color_manual(name   = "Method",values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
  )) +
  scale_linetype_manual(name   = "Variant",values = c("Shared" = "solid", "Unique" = "dashed")) +
  theme_bw()+custom_theme()+  theme(plot.title = element_blank(),legend.position  = "bottom")


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_enrichment_plot_bmi_updated_tidy_mesusiemlk_cadd_Shared_Unique.pdf"),
       plotA, width = 6, height = 3.9, dpi = 500, units = 'in')


### Create Panel B plot: Average Scores for continuous annotations (CADD and GTEx_eQTL_MaxCPP)
plotB <- ggplot(panelB, aes(x = threshold, y = enrichment, color = method, linetype = PIP)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(title = "Panel B: Average Scores for BMI",
       x = "PIP Threshold",
       y = "Average Score") +
  scale_color_manual(values = c(
    "MESuSiE" = "#009E73",   # Green
    "MLK"     = "#0072B2"    # Blue
  )) +
  scale_linetype_manual(values = c("Shared" = "solid", "Unique" = "dashed")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw()+custom_theme()+  theme(plot.title = element_blank(),legend.position  = "bottom")

### Combine the two panels vertically using patchwork
combined_plot <- plotA / plotB

# Display the combined plot
combined_plot


# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_enrichment_and_score_plot_bmi_updated_tidy_mesusiemlk_cadd_Shared_Unique.png"),
       combined_plot, width = 275, height = 300, dpi = 500, units = 'mm')

plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "combined_eqtl_cadd_annotations_plot_bmi_mesusiemlk_Shared_Unique.pdf"),
       plotB, width = 6, height = 3.9, dpi = 500, units = 'in')






##########################################################


####Combined Shared vs. Unique plot at PIP 0.5 and 0.9


##########################################################



selected_thresholds <- c(0.5)

### Split the data into two panels:
panelA <- enrichment_final %>% 
  filter(!(annotation %in% c("Z_eur"))) %>% dplyr::rename(PIP = variant) %>% filter(threshold == 0.5)

panelA_upper <- enrichment_final %>% 
  filter(!(annotation %in% c("Z_eur", "GTEx_eQTL_MaxCPP"))) %>% dplyr::rename(PIP = variant) %>% filter(threshold == 0.9)


### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.5
plotA <- ggplot(panelA, 
                aes(x = threshold, 
                    y = enrichment,
                    linetype = PIP,
                    fill = method 
                )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at PIP = 0.5",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
    
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
ggsave(paste0(plot_dir, "Figure_5_combined_enrichment_plot_bp_updated_tidy_mesusiemlk_cadd_Shared_Unique.png"),
       plotA, width = 10.7, height = 5.92, dpi = 500, units = 'in')


### Create Panel B plot: Enrichment (Signal / Background) for non-continuous annotations at PIIP 0.9
plotA_upper <- ggplot(panelA_upper, 
                      aes(x = threshold, 
                          y = enrichment,
                          linetype = PIP,
                          fill = method 
                      )) +
  geom_col(position = position_dodge(width = 0.8), 
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title    = "Shared and Ancestry-specific Enrichment at PIP = 0.9",
    y        = "Enrichment (Signal / Background)",
    fill     = "Method",       # legend title for fill
    linetype = "Variant"       # legend title for linetype
  ) +
  scale_fill_manual(values = c(
    "MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"
    
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
