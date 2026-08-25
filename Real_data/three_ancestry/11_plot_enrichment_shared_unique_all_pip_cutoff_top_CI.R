# ============================================================================
# 11_plot_enrichment_shared_unique_all_pip_cutoff_top_CI.R
#
# Copy of 11_plot_enrichment_shared_unique_all_pip_cutoff_top.R, EXCEPT the
# Shared vs. Unique functional enrichment is now computed with 95% confidence
# intervals, backed by chi-square / Fisher's exact testing of the enrichment
# ratio (relative risk). The CI logic mirrors compute_enrichment_stats() in
# 10_enrichment_chisq_CI_MF.R, and the bar plots gain matching error bars.
# ============================================================================

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

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summary_res/"
load(paste0(res_dir,"res_all_mf.RData"))

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region + 213, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region + 273, res_all$Region)
res_all$CS <- ifelse(res_all$CS != 0, 1, res_all$CS)
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
res_all$ASV <- ifelse(is.na(res_all$MAF_eur) | is.na(res_all$MAF_afr) | is.na(res_all$MAF_eas), 1, 0)



classify_signal <- function(shared, p1, p2, p3) {
    case_when(
      shared > 0.5 ~ 1,
      p1 > 0.5 & p2 < 0.5 & p3 < 0.5 ~ 2,  # EUR-specific
      p1 < 0.5 & p2 > 0.5 & p3 < 0.5 ~ 3,  # AFR-specific
      p1 < 0.5 & p2 < 0.5 & p3 > 0.5 ~ 4,  # EAS-specific
      .default = 0
    )
}

res_all <- res_all %>%
    mutate(
      MFD_Signal     = classify_signal(PIP_Shared, PIP_Ancestry_1, PIP_Ancestry_2, PIP_Ancestry_3),
      SuSiE_Signal   = classify_signal(SuSiE_PIP_Shared, SuSiE_PIP_EU, SuSiE_PIP_BB, SuSiE_PIP_EA),
      MESuSiE_Signal = classify_signal(MESuSiE_PIP_Shared, MESuSiE_PIP_EUR, MESuSiE_PIP_AFR, MESuSiE_PIP_EAS)
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

res_all = res_all %>% mutate(MFD_PIP_Unique = pmax(PIP_Ancestry_1, PIP_Ancestry_2, PIP_Ancestry_3),
                             SuSiE_PIP_Unique = pmax(SuSiE_PIP_EU, SuSiE_PIP_BB, SuSiE_PIP_EA),
                             MESuSiE_PIP_Unique = pmax(MESuSiE_PIP_EUR, MESuSiE_PIP_AFR, MESuSiE_PIP_EAS))

res_all <- res_all %>%
  mutate(
    MFD_PIP_Unique = case_when(
      MFD_PIP_Unique > 0.5 & !(MFD_Signal %in% c(2, 3, 4)) ~ 0,
      TRUE ~ MFD_PIP_Unique
    ),
    SuSiE_PIP_Unique = case_when(
      SuSiE_PIP_Unique > 0.5 & !(SuSiE_Signal %in% c(2, 3, 4)) ~ 0,
      TRUE ~ SuSiE_PIP_Unique
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


top_n_steps <- seq(10, 200, by = 10)
binary_ann  <- ann_col_name
cont_ann    <- c("Z_eur")

# Pre-calculate TOTAL sums for binary annotations (do this once globally)
total_n <- nrow(res_all)
total_sums_bin <- res_all %>%
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))

# ============================================================================
# HELPER: ENRICHMENT TEST WITH CI
# ============================================================================
#
# Constructs 2x2 table:
#              Annotated   Not-annotated
#   Signal     a           K - a
#   Background c_bg        N_bg - c_bg
#
# Returns: enrichment (RR = (a/K)/(c_bg/N_bg)); an APPROXIMATE 95% CI obtained by
#          taking a Wilson score interval for the top-K annotation proportion
#          a/K and scaling it by the background rate p_bg = c_bg/N_bg, which is
#          treated as fixed (valid here because N_bg is very large, so the
#          uncertainty in p_bg is much smaller than in a/K); a p-value from the
#          2x2 table (chi-square, or Fisher's exact when expected counts < 5);
#          and which test was used.
#
# NOTE: the CI is for the enrichment ratio (RR) and the p-value/test are computed
#       separately from the 2x2 table -- this is NOT Fisher's CI (which would be
#       for the odds ratio). If the background were small, a two-proportion score
#       CI for the RR (Koopman / Miettinen-Nurminen) or a bootstrap would be the
#       better, fully-unconditional choice.

compute_enrichment_stats <- function(a, K, c_bg, N_bg) {
  b <- K - a
  d <- N_bg - c_bg

  # POINT ESTIMATE: raw enrichment ratio = sig_freq / bg_freq, computed from the
  # uncorrected counts. This is IDENTICAL to the original script's enrichment, so
  # the plotted bar heights / lines are unchanged (e.g. an all-zero signal cell
  # still plots at 0).
  p_top <- a / K
  p_bg  <- c_bg / N_bg
  rr    <- p_top / p_bg

  # 95% CI: Wilson score interval for the top-K annotation proportion (a/K),
  # then scaled by the (treated-as-fixed) background annotation rate p_bg. This
  # is an APPROXIMATE 95% CI for the enrichment ratio, conditional on p_bg being
  # known -- reasonable here because N_bg is very large, so the uncertainty in
  # p_bg is much smaller than in a/K. It behaves correctly at a = 0 (lower bound
  # 0, interval contains the point estimate) and avoids the Haldane +0.5 fudge.
  z      <- qnorm(0.975)   # exact 97.5% normal quantile (= 1.959964), matches prop.test
  denom  <- 1 + z^2 / K
  centre <- (p_top + z^2 / (2 * K)) / denom
  halfw  <- (z * sqrt(p_top * (1 - p_top) / K + z^2 / (4 * K^2))) / denom
  ci_lo  <- max(0, (centre - halfw) / p_bg)   # clamp tiny FP noise to 0
  ci_hi  <- (centre + halfw) / p_bg

  # 2x2 table (rows = group, cols = annotation status)
  mat <- matrix(c(a, c_bg, b, d), nrow = 2)

  # Test selection: chi-square if all expected >= 5, else Fisher's exact
  expected  <- suppressWarnings(chisq.test(mat)$expected)
  if (any(expected < 5)) {
    pval      <- fisher.test(mat)$p.value
    test_type <- "Fisher"
  } else {
    pval      <- chisq.test(mat, correct = FALSE)$p.value
    test_type <- "Chi-square"
  }

  data.frame(
    enrichment = rr,
    ci_lower   = ci_lo,
    ci_upper   = ci_hi,
    p_value    = pval,
    test_used  = test_type,
    sig_count  = a,
    bg_count   = c_bg,
    stringsAsFactors = FALSE
  )
}

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

# 2. Compute Enrichment WITH CI / Statistical Tests --------------------------

enrichment_results_combined <- map_dfr(names(pip_config), function(variant_type) {

  # Get the map for this variant type (e.g., Shared methods)
  current_method_map <- pip_config[[variant_type]]

  map_dfr(names(current_method_map), function(method_name) {

    # Identify the correct column
    pip_col <- current_method_map[[method_name]]

    # Sort the dataframe ONCE per method/variant combo
    sorted_data <- res_all %>%
      select(all_of(c(pip_col, binary_ann, cont_ann))) %>%
      arrange(desc(.data[[pip_col]]))

    # Loop through the Top N steps
    map_dfr(top_n_steps, function(k) {

      # Slice the top k signals
      sig   <- sorted_data %>% slice(1:k)
      n_sig <- k
      n_bg  <- total_n - k

      # --- Signal / Background counts for binary annotations ---
      sig_sums_bin <- sig %>%
        summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))
      # Background Sum = Total Sum - Signal Sum
      bg_sums_bin <- total_sums_bin - sig_sums_bin

      # --- Binary annotations: enrichment ratio + 95% CI + test ---
      bin_stats <- map_dfr(binary_ann, function(ann) {
        a    <- sig_sums_bin[[ann]]
        c_bg <- bg_sums_bin[[ann]]
        stats <- compute_enrichment_stats(a, n_sig, c_bg, n_bg)
        stats$annotation <- ann
        stats
      })

      # --- Continuous Annotations: mean only (no test) ---
      cont_long <- sig %>%
        summarise(across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "annotation", values_to = "enrichment") %>%
        mutate(ci_lower = NA_real_, ci_upper = NA_real_,
               p_value = NA_real_, test_used = NA_character_,
               sig_count = NA_integer_, bg_count = NA_integer_)

      # --- Formatting & Combining ---
      bind_rows(bin_stats, cont_long) %>%
        mutate(
          method       = method_name,
          variant      = variant_type,
          top_n        = k
        )

    }) # End Top N Loop
  }) # End Method Loop
}) # End Variant Type Loop

# Per-annotation Bonferroni correction across methods within each top_n x variant
enrichment_results_combined <- enrichment_results_combined %>%
  group_by(top_n, variant, annotation) %>%
  mutate(
    n_tests      = sum(!is.na(p_value)),
    p_bonferroni = pmin(p_value * n_tests, 1)
  ) %>%
  ungroup() %>%
  mutate(
    sig_label = case_when(
      is.na(p_bonferroni)    ~ "",
      p_bonferroni < 0.001   ~ "***",
      p_bonferroni < 0.01    ~ "**",
      p_bonferroni < 0.05    ~ "*",
      TRUE                   ~ "ns"
    )
  )

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
  mutate(method = factor(method,
                         levels = c("MFD", "MESuSiE", "SuSiE"),
                         labels = c("MFD", "MESuSiE", "SuSiE")))

# Preview
head(enrichment_final)


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


plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"
system(paste0("mkdir -p ", plot_dir))
ggsave(paste0(plot_dir, "Supp_MFD_shared_unique_pip_5_to_9_CI.pdf"),
      plotA, width = 210, height = 180, dpi = 600, units = 'mm')


################################################################
library(openxlsx)
library(stringr) # For str_replace



# Replace .RData with .xlsx for the excel output
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
system(paste0("mkdir -p ",res_dir_table))
res_out_xlsx <-  "Real_data_enrichment_shared_unique_CI_3anc.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

sheet_name_raw <- str_remove(res_out_xlsx, ".xlsx")
sheet_name <- substr(sheet_name_raw, 1, 31)


panelA_to_supp = panelA %>% select(method,annotation,PIP, top_n,enrichment,
                                   ci_lower, ci_upper, p_value, p_bonferroni,
                                   test_used, sig_count, bg_count, sig_label)

wb <- createWorkbook()
addWorksheet(wb, sheet_name)

# Create a bold style for the Title
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# Write the Main Title at Row 1
writeData(wb, sheet_name, "Supplementary Table: Functional enrichment (with 95% CI) for priotized Shared and Unique Signal", startRow = 1)
addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)

# Write the Combined Data Table starting at Row 2
writeData(wb, sheet_name, panelA_to_supp, startRow = 2)

# 5. Save
saveWorkbook(wb, full_path, overwrite = TRUE)

print(paste("Saved combined file to:", full_path))
print(paste("Saved combined table to:", full_path))





##########################################################


####Combined Shared vs. Unique plot at PIP 0.5 and 0.9 (now with 95% CIs)


##########################################################




### Split the data into two panels:
panelA <- enrichment_final %>%
  filter(!(annotation %in% c("Z_eur"))) %>% dplyr::rename(PIP = variant) %>% filter(top_n == 50)

panelA_upper <- enrichment_final %>%
  filter(!(annotation %in% c("Z_eur", "GTEx_eQTL_MaxCPP"))) %>% dplyr::rename(PIP = variant) %>% filter(top_n == 100)


### Create Panel A plot: Enrichment (Signal / Background) for non-continuous annotations at top 50, with CIs
plotA <- ggplot(panelA,
                aes(x = top_n,
                    y = enrichment,
                    linetype = PIP,
                    fill = method,
                    group = interaction(method, PIP)
                )) +
  geom_col(position = position_dodge(width = 0.8),
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.3, linetype = "solid") +
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


plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"
ggsave(paste0(plot_dir, "Enrichment_pip_Shared_Unique_top50_combined_CI.pdf"),
       plotA, width = 210, height = 180, dpi = 600, units = 'mm')

ggsave(paste0(plot_dir, "Enrichment_pip_Shared_Unique_top50_combined_CI.pdf"),
       plotA, width = 10.7, height = 5.92, dpi = 600, units = 'in')


### Create Panel B plot: Enrichment (Signal / Background) for non-continuous annotations at top 100, with CIs
plotA_upper <- ggplot(panelA_upper,
                      aes(x = top_n,
                          y = enrichment,
                          linetype = PIP,
                          fill = method,
                          group = interaction(method, PIP)
                      )) +
  geom_col(position = position_dodge(width = 0.8),
           color    = "black",       # solid black border
           width    = 0.7) +         # narrower bars so dodge shows both
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.3, linetype = "solid") +
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
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"

ggsave(paste0(plot_dir, "combined_enrichment_plot_Shared_Unique_figure5_CI.pdf"),
       combined_plot, width = 7, height = 6.5, dpi = 500, units = 'in')


# ---- Summary table of test results at top 50 and top 100 ----
test_summary <- enrichment_final %>%
  filter(!(annotation %in% c("Z_eur")), top_n %in% c(50, 100)) %>%
  select(method, variant, annotation, top_n, enrichment, ci_lower, ci_upper,
         p_value, p_bonferroni, test_used, sig_count, bg_count, sig_label) %>%
  arrange(top_n, variant, annotation, method)

cat("\n=== Enrichment Test Results at top 50 / top 100 (Bonferroni-corrected) ===\n")
print(as.data.frame(test_summary), row.names = FALSE)

cat("\n=== Test type breakdown ===\n")
print(table(test_summary$test_used))

write.csv(test_summary,
          paste0(res_dir, "enrichment_test_results_shared_unique_top50_top100.csv"),
          row.names = FALSE)
cat("\nTest results saved to:",
    paste0(res_dir, "enrichment_test_results_shared_unique_top50_top100.csv"), "\n")




########################Classfication summary for MFD##################


bmi_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_BMI/result/Finemapping_Method_Decisions_BMI_3ancestry.csv")
dbp_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_DBP/result/Finemapping_Method_Decisions_DBP_3ancestry.csv")
sbp_class = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_SBP/result/Finemapping_Method_Decisions_SBP_3ancestry.csv")


# The CSV files already have proper headers (Locus_ID, Method, Reason, Status)
# and Method is already assigned — no merge with reference_table needed.

# Add a 'Trait' column to distinguish the data sources
bmi_class[, Trait := "BMI"]
dbp_class[, Trait := "DBP"]
sbp_class[, Trait := "SBP"]

# Row Bind the datasets
cols_to_keep <- c("Locus_ID", "Method", "Reason", "Trait")
all_data <- rbind(
  bmi_class[, ..cols_to_keep],
  dbp_class[, ..cols_to_keep],
  sbp_class[, ..cols_to_keep]
)


all_data$Locus_ID <- seq_len(nrow(all_data))


library(openxlsx)

res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
res_out_xlsx <- "Real_data_MFD_summary_3anc.xlsx"
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
