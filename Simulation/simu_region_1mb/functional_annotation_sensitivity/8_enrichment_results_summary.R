###############################################################################
#
#  8_enrichment_results_summary.R   (MFD revision / baseline_functional_annotation)
#
#  Adapted from shared_50_mesusie_inf/7_enrichment_results_summary.R
#
#  Differences vs INFUSE:
#   * Method set follows the new pipeline (no MESuSiE_mlk / MESuSiE_glm):
#       MESuSiE, SuSiE, SuSiE_weighted, SuSiE_merged, Paintor, Paintor_fun,
#       Paintor_all_fun, SuSiEx, XMAP, MultiSuSiE, CARMAX
#   * Per-SNP columns come from data_all (see 7_enrichment_preparation.R).
#   * Enrichment is evaluated for the 5 baseline BINARY annotations
#     (fold enrichment = freq_in_set / freq_out_of_set) PLUS evo2 as a
#     CONTINUOUS mean-ratio (mean_in_set / mean_out_of_set). Because evo2 is
#     continuous, the same sum/n formula yields the mean ratio automatically.
###############################################################################

suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/"
plot_dir <- paste0(base_dir, "plot/Enrichment/")
system(paste0("mkdir -p ", plot_dir))

# ---- method / column maps -------------------------------------------------
method_levels <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged",
                   "Paintor", "Paintor_fun", "Paintor_all_fun",
                   "SuSiEx", "XMAP", "MultiSuSiE", "CARMAX")

cs_cols <- c(MESuSiE = "MESuSiE_cs", SuSiE = "SuSiE_cs",
             SuSiE_weighted = "SuSiE_weighted_cs", SuSiE_merged = "SuSiE_merged_cs",
             Paintor = "Paintor_cs", Paintor_fun = "Paintor_fun_cs",
             Paintor_all_fun = "Paintor_all_fun_cs", SuSiEx = "SuSiEx_cs",
             XMAP = "XMAP_cs", MultiSuSiE = "MultiSuSIE_cs", CARMAX = "CARMAX_cs")

# "Either"-style overall PIP for each method (used for top-signal ranking)
pip_either <- c(MESuSiE = "MESuSiE_Either", SuSiE = "SuSiE_Either",
                SuSiE_weighted = "SuSiE_weighted_PIP", SuSiE_merged = "SuSiE_merged_PIP",
                Paintor = "Paintor_PIP", Paintor_fun = "Paintor_fun_PIP",
                Paintor_all_fun = "Paintor_all_fun_PIP", SuSiEx = "SuSiEx_PIP",
                XMAP = "XMAP_PIP", MultiSuSiE = "MultiSuSiE_PIP",
                CARMAX = "CARMAX_Either")

# Method colors aligned with simu_region_1mb/utility.R
method_colors <- c(MESuSiE = "#8da0cb", SuSiE = "#66c2a5",
                   SuSiE_weighted = "#B2D3A4", SuSiE_merged = "#9FBA95",
                   Paintor = "#fc8d62", Paintor_fun = "#e5703e",
                   Paintor_all_fun = "#d4533a", SuSiEx = "#E89DA0",
                   XMAP = "#ffd92f", MultiSuSiE = "#e78ac3", CARMAX = "#f2b56e")

binary_ann <- c("H3K27ac_Hnisz_common", "DHS_Trynka_common", "H3K4me3_Trynka_common",
                "Conserved_LindbladToh_common", "non_synonymous_common")
eval_ann <- c(binary_ann, "evo2_score")   # evo2 last -> continuous mean ratio

ann_levels <- c("non_synonymous_common", "Conserved_LindbladToh_common",
                "DHS_Trynka_common", "H3K4me3_Trynka_common",
                "H3K27ac_Hnisz_common", "evo2_score")
ann_recode <- c("non_synonymous_common" = "Non-synonymous",
                "Conserved_LindbladToh_common" = "Conserved",
                "H3K4me3_Trynka_common" = "H3K4me3",
                "DHS_Trynka_common" = "DHS",
                "H3K27ac_Hnisz_common" = "H3K27ac",
                "evo2_score" = "evo2")

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
    strip.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, face = "bold"),
    plot.title = element_text(size = 7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )
}

# ---- enrichment helpers ---------------------------------------------------
# 95% credible-set fold enrichment (marginal): among regions with >=1 set SNP,
# freq(annotation | in set) / freq(annotation | not in set).
calc_cs_enrichment <- function(df, cs_col, ann) {
  df %>%
    group_by(Region) %>%
    filter(sum(.data[[cs_col]], na.rm = TRUE) != 0) %>%
    group_by(.data[[cs_col]]) %>%
    summarise(across(all_of(ann), ~ sum(.x, na.rm = TRUE) / n()), .groups = "drop") %>%
    summarise(across(all_of(ann),
                     ~ .x[.data[[cs_col]] == 1] / .x[.data[[cs_col]] == 0]))
}

all_sets <- list()
all_CS_enrichment <- list()
all_Signal_enrichment <- list()

for (num_causal in 1:3) {
  for (h2_num in 1:2) {

    out_dir <- paste0(base_dir, "causal_num_", num_causal, "/out/")
    load(paste0(out_dir, "res_simulation_h2_", h2_num, ".RData"))   # res_all

    ############################################################
    #  Credible-set fold enrichment (per method)
    ############################################################
    CS_enrichment <- bind_rows(lapply(method_levels, function(m) {
      calc_cs_enrichment(res_all, cs_cols[[m]], eval_ann) %>% mutate(Method = m)
    }))

    CS_enrichment_long <- CS_enrichment %>%
      pivot_longer(cols = -Method, names_to = "Cat", values_to = "Prop") %>%
      mutate(Method = factor(Method, levels = method_levels),
             Cat = recode(factor(Cat, levels = ann_levels), !!!ann_recode),
             Prop = round(Prop, 2))

    all_CS_enrichment[[paste0("causal", num_causal, "_h2", h2_num)]] <-
      CS_enrichment_long %>% mutate(causal = num_causal, h2 = h2_num,
                                    enrichment_type = "CS")

    ############################################################
    #  Top-signal fold enrichment (top-N by overall/Either PIP)
    ############################################################
    top_N_signal <- ifelse(num_causal == 3, 150, 100)
    bg_an <- res_all %>% summarise(across(all_of(eval_ann),
                                          ~ sum(.x, na.rm = TRUE) / n()))

    Signal_enrichment <- bind_rows(lapply(method_levels, function(m) {
      pcol <- pip_either[[m]]
      sig <- res_all %>%
        filter(.data[[pcol]] > 0.5) %>%
        arrange(desc(.data[[pcol]])) %>%
        slice_head(n = top_N_signal)
      if (nrow(sig) == 0) {
        enr <- bg_an; enr[] <- NA_real_
      } else {
        enr <- (sig %>% summarise(across(all_of(eval_ann),
                                         ~ sum(.x, na.rm = TRUE) / n()))) / bg_an
      }
      enr %>% mutate(Method = m)
    }))

    Signal_enrichment_long <- Signal_enrichment %>%
      pivot_longer(cols = -Method, names_to = "Cat", values_to = "Prop") %>%
      mutate(Method = factor(Method, levels = method_levels),
             Cat = recode(factor(Cat, levels = ann_levels), !!!ann_recode),
             Prop = round(Prop, 2))

    all_Signal_enrichment[[paste0("causal", num_causal, "_h2", h2_num)]] <-
      Signal_enrichment_long %>% mutate(causal = num_causal, h2 = h2_num,
                                        enrichment_type = "Signal")

    ############################################################
    #  Set size (per region, per method)
    ############################################################
    set_info <- res_all %>%
      group_by(Region) %>%
      summarise(across(all_of(unname(cs_cols)), ~ sum(.x, na.rm = TRUE)),
                .groups = "drop")
    set_info_long <- set_info %>%
      pivot_longer(-Region, names_to = "cs_col", values_to = "Count") %>%
      mutate(Method = factor(names(cs_cols)[match(cs_col, cs_cols)],
                             levels = method_levels)) %>%
      filter(Count > 0)
    all_sets[[paste0("causal", num_causal, "_h2", h2_num)]] <-
      set_info_long %>% mutate(causal = num_causal, h2 = h2_num)

    ############################################################
    #  Per-iteration plots
    ############################################################
    p_set <- ggplot(set_info_long, aes(x = Method, y = Count, fill = Method)) +
      geom_boxplot(outlier.size = 0.1, fatten = 0.5, color = "darkgray") +
      scale_fill_manual(values = method_colors, guide = "none") +
      theme_bw() + xlab("") + ylab("Set Size") +
      coord_cartesian(ylim = c(0, 100)) + custom_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    p_cs <- ggplot(CS_enrichment_long, aes(x = Cat, y = Prop, fill = Method)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = method_colors) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      xlab("") + ylab("Fold Enrichment Credible Set") +
      theme_bw() + custom_theme() + guides(fill = guide_legend(nrow = 2))
    p_out <- p_set / p_cs +
      plot_annotation(tag_levels = 'a') +
      plot_layout(guides = "collect", heights = c(1, 1)) &
      theme(legend.position = 'bottom', plot.tag = element_text(size = 6, face = "bold"))
    ggsave(paste0(plot_dir, "Causal_num_", num_causal, "_h2_", h2_num,
                  "_Enrichment_Comparison_Set_CS.pdf"),
           p_out, dpi = 500, width = 200, height = 130, units = 'mm')

    p_sig <- ggplot(Signal_enrichment_long, aes(x = Cat, y = Prop, fill = Method)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = method_colors) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      xlab("") + ylab("Fold Enrichment Top Signal") +
      theme_bw() + custom_theme() + guides(fill = guide_legend(nrow = 2))
    p_out <- p_set / p_sig +
      plot_annotation(tag_levels = 'a') +
      plot_layout(guides = "collect", heights = c(1, 1)) &
      theme(legend.position = 'bottom', plot.tag = element_text(size = 7, face = "bold"))
    ggsave(paste0(plot_dir, "Causal_num_", num_causal, "_h2_", h2_num,
                  "_Enrichment_Comparison_Set_Top.pdf"),
           p_out, dpi = 500, width = 200, height = 130, units = 'mm')

    cat("done causal", num_causal, "h2", h2_num, "\n")
  }
}

###########################################################################
#  Summary plots averaged across all causal_num / h2
###########################################################################
combined_CS <- bind_rows(all_CS_enrichment)
CS_combined <- combined_CS %>%
  group_by(Method, Cat) %>%
  summarise(Prop = round(mean(Prop, na.rm = TRUE), 2), .groups = "drop") %>%
  mutate(Method = factor(Method, levels = method_levels),
         Cat = factor(Cat, levels = unname(ann_recode)))

combined_Signal <- bind_rows(all_Signal_enrichment)
Signal_combined <- combined_Signal %>%
  group_by(Method, Cat) %>%
  summarise(Prop = round(mean(Prop, na.rm = TRUE), 2), .groups = "drop") %>%
  mutate(Method = factor(Method, levels = method_levels),
         Cat = factor(Cat, levels = unname(ann_recode)))

combined_sets <- bind_rows(all_sets)
sets_combined <- combined_sets %>%
  group_by(Method, Region) %>%
  summarise(Count = round(mean(Count), 2), .groups = "drop") %>%
  mutate(Method = factor(Method, levels = method_levels))

p_set <- ggplot(sets_combined, aes(x = Method, y = Count, fill = Method)) +
  geom_boxplot(outlier.size = 0.1, fatten = 0.5, color = "darkgray") +
  scale_fill_manual(values = method_colors, guide = "none") +
  theme_bw() + xlab("") + ylab("Set Size") +
  coord_cartesian(ylim = c(0, 100)) + custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cs <- ggplot(CS_combined, aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colors) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("") + ylab("Fold Enrichment Credible Set") +
  theme_bw() + custom_theme() + guides(fill = guide_legend(nrow = 2))
p_out <- p_set / p_cs + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = 'bottom', plot.tag = element_text(size = 6, face = "bold"))
ggsave(paste0(plot_dir, "Summary_Enrichment_Comparison_Set_CS.pdf"),
       p_out, dpi = 500, width = 200, height = 130, units = 'mm')

p_sig <- ggplot(Signal_combined, aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colors) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("") + ylab("Fold Enrichment Top Signal") +
  theme_bw() + custom_theme() + guides(fill = guide_legend(nrow = 2))
p_out <- p_set / p_sig + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = 'bottom', plot.tag = element_text(size = 7, face = "bold"))
ggsave(paste0(plot_dir, "Summary_Enrichment_Comparison_Set_Top.pdf"),
       p_out, dpi = 500, width = 200, height = 130, units = 'mm')

# save tables
save(combined_CS, combined_Signal, combined_sets,
     CS_combined, Signal_combined, sets_combined,
     file = paste0(plot_dir, "enrichment_summary_tables.RData"))

cat("\n--- Combined CS fold enrichment ---\n"); print(as.data.frame(CS_combined))
cat("\n--- Combined Top-signal fold enrichment ---\n"); print(as.data.frame(Signal_combined))
cat("\nDone. Plots + tables written to ", plot_dir, "\n")
