# ============================================================================
# 10_enrichment_chisq_CI_MF.R
#
# Full real-data plotting pipeline (line plot, CS-count plot, set-size boxplots,
# Z-score plot, combined 4-panel figure, by-trait enrichment, supplementary rank
# tables) -- identical to 10_plot_enrichment_either_pip_cutoff_top_draft_MF.R,
# EXCEPT the plain functional enrichment bar plot is replaced by the updated
# bar plot with 95% confidence intervals, backed by chi-square / Fisher's exact
# testing of the enrichment ratio (relative risk).
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
library(openxlsx)

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summary_res/"
load(paste0(res_dir,"res_all_mf.RData"))

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region + 213, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region + 273, res_all$Region)
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
  #CARMAX       = "CARMAX_PIP_Either",
  SuSiE_merged = "SuSiE_merged_PIP_Either"
)

method_colors <- c(
  "MESuSiE" = "#8da0cb", "SuSiE" = "#66c2a5", "MFD" = "#377eb8",
  "SuSiE_Meta" = "#9FBA95", "Paintor" = "#fc8d62", "MultiSuSiE" = "#e78ac3",
  "SuSiEx" = "#E89DA0", "XMAP" = "#ffd92f", "CARMAx" = "#f2b56e"
)

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

  # POINT ESTIMATE: raw enrichment ratio = sig_freq / bg_freq, from the
  # uncorrected counts. This is IDENTICAL to the original draft script's
  # enrichment, so the plotted bar heights / lines are unchanged (an all-zero
  # signal cell still plots at 0).
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

# ============================================================================
# MAIN ENRICHMENT LOOP WITH STATISTICAL TESTS
# ============================================================================

total_n <- nrow(res_all)
total_sums_bin <- res_all %>%
  summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))

enrichment_results <- map_dfr(names(pip_cols), function(method) {

  pip_col <- pip_cols[[method]]

  sorted_data <- res_all %>%
    select(all_of(c(pip_col, binary_ann, cont_ann))) %>%
    arrange(desc(.data[[pip_col]]))

  map_dfr(top_n_steps, function(k) {

    sig    <- sorted_data %>% slice(1:k)
    n_sig  <- k
    n_bg   <- total_n - k

    sig_sums_bin <- sig %>%
      summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))
    bg_sums_bin  <- total_sums_bin - sig_sums_bin

    # Binary annotations: enrichment + statistical test
    bin_stats <- map_dfr(binary_ann, function(ann) {
      a    <- sig_sums_bin[[ann]]
      c_bg <- bg_sums_bin[[ann]]
      stats <- compute_enrichment_stats(a, n_sig, c_bg, n_bg)
      stats$annotation <- ann
      stats
    })

    # Continuous annotations: mean only (no test)
    cont_long <- sig %>%
      summarise(across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "annotation", values_to = "enrichment") %>%
      mutate(ci_lower = NA_real_, ci_upper = NA_real_,
             p_value = NA_real_, test_used = NA_character_,
             sig_count = NA_integer_, bg_count = NA_integer_)

    bind_rows(bin_stats, cont_long) %>%
      mutate(method = method, top_n = k)
  })
})

# ============================================================================
# BONFERRONI CORRECTION & FORMATTING
# ============================================================================

# Per-annotation Bonferroni correction across 7 methods within each top_n step
enrichment_results <- enrichment_results %>%
  group_by(top_n, annotation) %>%
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

# Reorder the 'annotation' factor in enrichment_results
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation,
                             levels = c("non_synonymous_common",
                                        "synonymous",
                                        "promotor","utr_comb",
                                        "CRE",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL","brain_ind_eQTL","Z_eur","Z_afr","Z_eas"))) %>%
  mutate(annotation = recode(annotation,
                             "non_synonymous_common" = "Non-synonymous",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "CRE" = "CRE",
                             "heart_ind_eQTL" = "Heart ind. eQTL",
                             "artery_ind_eQTL" = "Artery ind. eQTL",
                             "brain_ind_eQTL" = "Brain ind. eQTL",
                             "Z_eur" = "Z_eur",
                             "Z_afr" = "Z_afr",
                             "Z_eas" = "Z_eas"))

enrichment_results$method <- factor(
  enrichment_results$method,
  # 1. 'levels' must match the CURRENT strings in the data frame
  levels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged","Paintor", "SuSiEx", "XMAP"),
  # 2. 'labels' are the names assigned to the levels
  labels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta","Paintor", "SuSiEx", "XMAP")
)

# Create the plot with colorblind-friendly colors
enrichment_plot <- ggplot(enrichment_results%>%
                            filter(!(annotation %in% c("Z_eur","Z_afr","Z_eas"))), aes(x = top_n, y = enrichment, color = method)) +
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
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"
ggsave(paste0(plot_dir, "Supp_MFD_Either_pip_5_to_9.pdf"),
       enrichment_plot, width = 210, height = 180, dpi = 600, units = 'mm')





######################################


#Barplot at PIP threshold 0.5 (bar plot with 95% confidence intervals,
# backed by chi-square / Fisher's exact testing of the enrichment ratio)


######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

fill_colors =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#377eb8","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")

# Create Panel A plot: Enrichment (Signal / Background) at K=100, with CIs
enrichment_panelA = enrichment_results %>% filter(top_n==100) %>% filter(!(annotation %in% c("Z_eur","Z_afr","Z_eas")))

# Bar plot: error bars are 95% confidence intervals on the enrichment ratio
plotA <- ggplot(enrichment_panelA, aes(x = method, y = enrichment, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, linewidth = 0.4) +

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

print(plotA)

ggsave(paste0(plot_dir, "Supp_MFD_Either_pip_top100_barplot_CI.pdf"),
       plotA, width = 210, height = 180, dpi = 600, units = 'mm')

# ---- Summary table of test results at K=100 ----
test_summary <- enrichment_panelA %>%
  select(method, annotation, enrichment, ci_lower, ci_upper,
         p_value, p_bonferroni, test_used, sig_count, bg_count, sig_label) %>%
  arrange(annotation, method)

cat("\n=== Enrichment Test Results at K=100 (Bonferroni-corrected) ===\n")
print(as.data.frame(test_summary), row.names = FALSE)

cat("\n=== Test type breakdown ===\n")
print(table(test_summary$test_used))

cat("\n=== Significant results (Bonferroni p < 0.05) ===\n")
sig_results <- test_summary %>% filter(p_bonferroni < 0.05)
cat(nrow(sig_results), "out of", nrow(test_summary), "tests significant\n")
print(as.data.frame(sig_results), row.names = FALSE)

write.csv(test_summary,
          paste0(res_dir, "enrichment_test_results_k100.csv"),
          row.names = FALSE)
cat("\nTest results saved to:", paste0(res_dir, "enrichment_test_results_k100.csv"), "\n")

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
  dplyr::rename(
    val_coding = Coding,
    val_regulatory = Regulatory, # Regulatory annotation
    val_eqtl = eQTL,
    id_map = method # Join key for method summaries
  )


#####Data for rank based evaluation#####

# --- 1. Setup Paths and Levels ---
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
system(paste0("mkdir -p ", res_dir_table))
res_out_xlsx <- "Real_data_enrichment_either_rank_revision_3anc.xlsx"
full_path <- paste0(res_dir_table, res_out_xlsx)

# Define the order from the factor levels
desired_levels <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP")

# --- 2. Process Enrichment Data (Supplementary Table 13) ---
# Report fold enrichment AND the 95% CI + p-values for every
# method x annotation x top-ranked threshold (top_n = 10..200).
enrichment_results_to_supp <- enrichment_results %>%
  filter(!(annotation %in% c("Z_eur"))) %>%
  mutate(method = factor(method, levels = desired_levels)) %>%
  select(method, annotation, top_n,
         fold_enrichment = enrichment,
         ci_lower, ci_upper,
         p_value, p_bonferroni, test_used,
         sig_count, bg_count, sig_label) %>%
  arrange(method, annotation, top_n)

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

# Sheet 1: Enrichment Results (Supplementary Table 13)
addWorksheet(wb, "Enrichment_Results")
writeData(wb, "Enrichment_Results",
          "Supplementary Table 13: Functional enrichment of fine-mapped variants across varying top-ranked thresholds (PIP_either)",
          startRow = 1)
addStyle(wb, "Enrichment_Results", title_style, rows = 1, cols = 1)
writeData(wb, "Enrichment_Results",
          paste0("Functional enrichment metrics for variants identified by evaluated fine-mapping methods, ",
                 "aggregated across three quantitative traits (BMI, DBP, SBP). Eight functional annotation categories: ",
                 "non-synonymous, synonymous, promoter, UTR, cis-regulatory elements (CRE), and independent eQTLs ",
                 "(heart, artery, brain). Fold enrichment is the ratio of the proportion of annotated variants within ",
                 "the selected top set to the proportion of that annotation in the background set, reported as a function ",
                 "of the number of top-ranked variants (top 10 to 200 based on PIP_either). ci_lower/ci_upper give the ",
                 "95% Wilson-scaled confidence interval for the enrichment ratio; p_value and p_bonferroni are from the ",
                 "corresponding 2x2 chi-square / Fisher's exact test (test_used)."),
          startRow = 2)
writeData(wb, "Enrichment_Results", enrichment_results_to_supp, startRow = 4)

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
region_cs<-res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE)))

cs_number_by_trait = region_cs%>%group_by(PHENONAME) %>% summarise(sum(CS!=0),sum(MESuSiE_cs!=0),sum(SuSiE_cs!=0),sum(Paintor_cs!=0),sum(SuSiEx_cs!=0),sum(XMAP_cs !=0),sum(SuSiE_metal_cs !=0))

data_a <- cs_number_by_trait %>%
  # Rename grouping column
  dplyr::rename(Trait = PHENONAME) %>%

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
   # grepl("CARMAX", raw_col)      ~ "CARMAx",

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
                                                      "MultiSuSiE", "XMAP")))

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


all_sets_info<-data.frame(res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE))))#%>%filter(CS!= 0,MESuSiE_cs!=0,SuSiE_cs!=0,Paintor_cs!=0,SuSiEx_cs !=0,XMAP_cs!=0,CARMAX_cs!=0) #Median Set Size across all locus

# 1. Prepare the data (assuming all_sets_info is already loaded)
all_sets_info_long <- all_sets_info %>%
  # Select only the required columns (excluding the indices if they exist as row names)
  select(PHENONAME, Region, CS, MESuSiE_cs, SuSiE_cs, Paintor_cs,
         SuSiEx_cs, XMAP_cs, SuSiE_metal_cs) %>%
  # Pivot to long format
  pivot_longer(
    cols = c("CS", "MESuSiE_cs", "SuSiE_cs", "Paintor_cs",
             "SuSiEx_cs", "XMAP_cs", "SuSiE_metal_cs"),
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
                                    "SuSiEx_cs", "XMAP_cs"),
                         labels = c("MFD", "MESuSiE", "SuSiE","SuSiE_Meta", "Paintor",
                                    "SuSiEx", "XMAP")))


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
MFD_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CS==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))
MESuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
SuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
Paintor_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
SuSiEx_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiEx_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
XMAP_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(XMAP_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
#CARMAX_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CARMAX_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)
SuSiE_metal_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_metal_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),abs(Z_eas),na.rm=T)))%>%pull(zmax)

set_size_z_info<-data.frame(cbind(MFD_cs_Z,MESuSiE_cs_Z,SuSiE_cs_Z,Paintor_cs_Z,SuSiEx_cs_Z,XMAP_cs_Z,SuSiE_metal_cs_Z))
colnames(set_size_z_info)<-c("PHENONAME",c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","SuSiE_Meta"))
set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(PHENONAME), names_to = "Method", values_to = "Z")%>%mutate(Method = factor(Method, levels=c("MFD","MESuSiE","SuSiE", "SuSiE_Meta","Paintor","SuSiEx","XMAP")))

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
p_d = plotA   # enrichment bar plot with 95% CIs


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


plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"
ggsave(paste0(plot_dir, "set_size_enrichment_pip_top100_combined_updated_test.pdf"),
       combined_plot, width = 550, height = 300, dpi = 600, units = 'mm')

#ggsave(paste0(plot_dir, "set_size_enrichment_pip_top100_combined_updated_test.pdf"),
#       combined_plot, width = 210, height = 180, dpi = 600, units = 'mm')






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
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summary_res/"
load(paste0(res_dir,"res_all_mf.RData"))

# Pre-processing
res_all <- res_all %>%
  mutate(
    n_regions = n_distinct(Region), # Calculate distinct regions if needed later
    Region = case_when(
      PHENONAME == "DBP" ~ Region + 213,
      PHENONAME == "SBP" ~ Region + 273,
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
  #CARMAX       = "CARMAX_PIP_Either",
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
    dplyr::rename(PIP = !!sym(pip_col))

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

  # Calculate enrichment with CI per phenotype x annotation
  bin_results <- map_dfr(1:nrow(merged), function(i) {
    row <- merged[i, ]
    n_sig <- target_k
    n_bg  <- row$total_n - target_k

    map_dfr(binary_ann, function(ann) {
      sig_col <- paste0("sig_", ann)
      tot_col <- paste0("total_", ann)
      a    <- row[[sig_col]]
      c_bg <- row[[tot_col]] - a

      stats <- compute_enrichment_stats(a, n_sig, c_bg, n_bg)
      stats$annotation <- ann
      stats$PHENONAME  <- row$PHENONAME
      stats
    })
  })

  # Continuous annotations: mean only (no test)
  cont_results <- top_hits %>%
    group_by(PHENONAME) %>%
    summarise(
      across(all_of(cont_ann), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = all_of(cont_ann), names_to = "annotation", values_to = "enrichment") %>%
    mutate(ci_lower = NA_real_, ci_upper = NA_real_,
           p_value = NA_real_, test_used = NA_character_,
           sig_count = NA_integer_, bg_count = NA_integer_)

  long_res <- bind_rows(bin_results, cont_results) %>%
    mutate(method = method, top_n = target_k)

  return(long_res)
})

# --- 4. Bonferroni Correction ---
# Per-annotation x PHENONAME Bonferroni correction across methods
enrichment_results <- enrichment_results %>%
  group_by(PHENONAME, annotation) %>%
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

# --- 5. Formatting for Plotting ---
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
                    levels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_merged", "Paintor", "SuSiEx", "XMAP" ),
                    labels = c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP"))
  )


# Define Colors
method_colors <- c(
  "MESuSiE" = "#8da0cb", "SuSiE" = "#66c2a5", "MFD" = "#377eb8",
  "SuSiE_Meta" = "#9FBA95", "Paintor" = "#fc8d62", "MultiSuSiE" = "#e78ac3",
  "SuSiEx" = "#E89DA0", "XMAP" = "#ffd92f"
)

enrichment_plot <- ggplot(enrichment_results %>% filter(!(annotation %in% c("Z_eur"))), aes(x = PHENONAME, y = enrichment, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8), width = 0.25) +
  facet_wrap(~ annotation, scales = "free_y", ncol = 3) +
  labs(
    title = paste0("Enrichment Results (Top ", target_k, " signals) by Phenotype"),
    x = "Phenotype",
    y = "Enrichment (Signal / Background)"
  ) +
  scale_fill_manual(values = method_colors) +
  theme_bw() +  custom_theme()+theme(plot.title = element_blank(),legend.position  = "bottom")

print(enrichment_plot)

plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Figure/CI_updated/"
ggsave(paste0(plot_dir, "Supp_MFD_Either_pip_top_100_by_trait.pdf"),
       enrichment_plot, width = 210, height = 180, dpi = 600, units = 'mm')



#####Data for rank based evaluation#####

# --- 1. Setup Paths and Levels ---
res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
# No need to mkdir again if it exists, but safe to keep
system(paste0("mkdir -p ", res_dir_table))

res_out_xlsx_eqtl <- "Real_data_enrichment_either_eQTL_revision_3anc.xlsx"
full_path_eqtl <- paste0(res_dir_table, res_out_xlsx_eqtl)

# Define the order from the factor levels
desired_levels <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_Meta", "Paintor", "SuSiEx", "XMAP")

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
