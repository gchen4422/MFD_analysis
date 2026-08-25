###############################################################################
#
#  9_2_enrichment_results_pip_cutoff_shared.R
#  (MFD revision / baseline_functional_annotation)
#
#  Adapted from shared_50_mesusie_inf/9_2_enrichmen_results_pip_cutoff_shared.R
#
#  SHARED-signal enrichment. For each threshold, top_N_signal is the number of
#  SNPs whose MESuSiE shared PIP exceeds the threshold; each method then takes
#  its own top_N_signal SNPs (ranked by its shared-PIP column, or its single
#  overall PIP for methods without an ancestry decomposition) and evaluate
#  enrichment vs background.
#   - 5 baseline binary annotations: fold enrichment
#   - evo2: continuous mean-ratio
###############################################################################

suppressMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(purrr); library(stringr)
})

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/"
plot_dir <- paste0(base_dir, "plot/Enrichment/")
system(paste0("mkdir -p ", plot_dir))

binary_ann <- c("H3K27ac_Hnisz_common", "DHS_Trynka_common", "H3K4me3_Trynka_common",
                "Conserved_LindbladToh_common", "non_synonymous_common")
eval_ann <- c(binary_ann, "evo2_score")

# shared-signal PIP column per method (single-PIP methods use overall PIP)
pip_cols <- c(MESuSiE = "MESuSiE_Shared", SuSiE = "SuSiE_Shared",
              SuSiE_weighted = "SuSiE_weighted_PIP", SuSiE_merged = "SuSiE_merged_PIP",
              Paintor = "Paintor_PIP", Paintor_fun = "Paintor_fun_PIP",
              Paintor_all_fun = "Paintor_all_fun_PIP", SuSiEx = "SuSiEx_PIP",
              XMAP = "XMAP_PIP", MultiSuSiE = "MultiSuSiE_PIP",
              CARMAX = "CARMAX_Shared")
method_levels <- names(pip_cols)
ref_pip <- "MESuSiE_Shared"   # defines top_N_signal

# Method colors aligned with simu_region_1mb/utility.R
method_colors <- c(MESuSiE = "#8da0cb", SuSiE = "#66c2a5",
                   SuSiE_weighted = "#B2D3A4", SuSiE_merged = "#9FBA95",
                   Paintor = "#fc8d62", Paintor_fun = "#e5703e",
                   Paintor_all_fun = "#d4533a", SuSiEx = "#E89DA0",
                   XMAP = "#ffd92f", MultiSuSiE = "#e78ac3", CARMAX = "#f2b56e")

ann_levels <- c("non_synonymous_common", "Conserved_LindbladToh_common",
                "DHS_Trynka_common", "H3K4me3_Trynka_common",
                "H3K27ac_Hnisz_common", "evo2_score")
ann_recode <- c("non_synonymous_common" = "Non-synonymous",
                "Conserved_LindbladToh_common" = "Conserved",
                "H3K4me3_Trynka_common" = "H3K4me3", "DHS_Trynka_common" = "DHS",
                "H3K27ac_Hnisz_common" = "H3K27ac", "evo2_score" = "evo2")

thresholds <- seq(0.5, 0.9, by = 0.05)
all_Signal_enrichment <- list()

for (num_causal in 1:3) {
  for (h2_num in 1:2) {

    out_dir <- paste0(base_dir, "causal_num_", num_causal, "/out/")
    load(paste0(out_dir, "res_simulation_h2_", h2_num, ".RData"))   # res_all

    enrichment_results <- map_df(thresholds, function(thresh) {
      top_N_signal <- res_all %>% filter(.data[[ref_pip]] > thresh) %>% nrow()

      if (top_N_signal == 0) {
        return(map_df(method_levels, function(method) {
          data.frame(annotation = eval_ann, enrichment = NA_real_,
                     method = method, threshold = thresh)
        }))
      }

      bg <- res_all %>%
        summarise(across(all_of(eval_ann), ~ sum(.x, na.rm = TRUE) / (n() - top_N_signal)))

      map_df(method_levels, function(method) {
        pcol <- pip_cols[[method]]
        sig <- res_all %>% arrange(desc(.data[[pcol]])) %>%
          slice_head(n = top_N_signal)
        sig_freq <- sig %>% summarise(across(all_of(eval_ann), ~ sum(.x, na.rm = TRUE) / n()))
        enrich <- sig_freq / bg
        enrich %>%
          pivot_longer(cols = everything(), names_to = "annotation", values_to = "enrichment") %>%
          mutate(method = method, threshold = thresh)
      })
    })

    enrichment_results <- enrichment_results %>%
      mutate(annotation = recode(factor(annotation, levels = ann_levels), !!!ann_recode),
             method = factor(method, levels = method_levels),
             causal = num_causal, h2 = h2_num)

    p <- ggplot(enrichment_results, aes(x = threshold, y = enrichment, color = method)) +
      geom_line() +
      facet_wrap(~ annotation, scales = "free_y") +
      scale_color_manual(values = method_colors) +
      labs(title = paste("Shared-signal enrichment | causal =", num_causal, "h2 =", h2_num),
           x = "MESuSiE shared-PIP threshold", y = "Enrichment (Signal / Background)") +
      theme_bw()
    ggsave(paste0(plot_dir, "Causal_num_", num_causal, "_h2_", h2_num,
                  "_Enrichment_pip_cutoff_shared.pdf"),
           p, dpi = 500, width = 200, height = 130, units = 'mm')

    all_Signal_enrichment[[paste0("num_causal_", num_causal, "_h2_", h2_num)]] <- enrichment_results
    cat("done causal", num_causal, "h2", h2_num, "\n")
  }
}

combined_enrichment <- bind_rows(all_Signal_enrichment)
avg_enrichment <- combined_enrichment %>%
  group_by(method, annotation, threshold) %>%
  summarise(avg_enrichment = mean(enrichment, na.rm = TRUE), .groups = "drop")

avg_plot <- ggplot(avg_enrichment, aes(x = threshold, y = avg_enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  scale_color_manual(values = method_colors) +
  labs(title = "Average shared-signal enrichment across all causal_num / h2",
       x = "MESuSiE shared-PIP threshold", y = "Average Enrichment (Signal / Background)") +
  theme_bw()
ggsave(paste0(plot_dir, "Summary_Enrichment_pip_cutoff_shared.pdf"),
       avg_plot, dpi = 500, width = 200, height = 130, units = 'mm')

save(combined_enrichment, avg_enrichment,
     file = paste0(plot_dir, "enrichment_pip_cutoff_shared.RData"))
cat("Done (shared). Output in ", plot_dir, "\n")
