###############################################################################
#
#  9_1_enrichment_results_pip_cutoff_either_pip_method.R
#  (MFD revision / baseline_functional_annotation)
#
#  Adapted from shared_50_mesusie_inf/9_1_enrichmen_results_pip_cutoff_either_pip_method.R
#
#  Enrichment of SNPs with EITHER-ancestry PIP > threshold, swept over
#  thresholds 0.5..0.9, for every method in the new pipeline.
#  - 5 baseline binary annotations: fold enrichment (freq_signal/freq_bg)
#  - evo2: continuous mean-ratio (same sum/n formula -> mean ratio)
#
#  NOTE: the INFUSE version also computed an MSE-vs-true-enrichment panel that
#  relied on a precomputed True_enrichment_w2_either.Rdata. No such truth file
#  exists for this simulation, so that block is intentionally omitted.
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

# overall / Either PIP column per method
pip_cols <- c(MESuSiE = "MESuSiE_Either", SuSiE = "SuSiE_Either",
              SuSiE_weighted = "SuSiE_weighted_PIP", SuSiE_merged = "SuSiE_merged_PIP",
              Paintor = "Paintor_PIP", Paintor_fun = "Paintor_fun_PIP",
              Paintor_all_fun = "Paintor_all_fun_PIP", SuSiEx = "SuSiEx_PIP",
              XMAP = "XMAP_PIP", MultiSuSiE = "MultiSuSiE_PIP",
              CARMAX = "CARMAX_Either")
method_levels <- names(pip_cols)

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

# ============================================================================
# HELPER: ENRICHMENT TEST WITH CI (from 10_enrichment_chisq_CI_MF.R)
#
# Constructs 2x2 table for binary annotations:
#              Annotated   Not-annotated
#   Signal     a           K - a
#   Background c_bg        N_bg - c_bg
#
# Returns: enrichment (RR); Wilson-scaled 95% CI; chi-square / Fisher p-value.
# ============================================================================
compute_enrichment_stats <- function(a, K, c_bg, N_bg) {
  p_top <- a / K
  p_bg  <- c_bg / N_bg
  rr    <- p_top / p_bg

  # Wilson score interval for top-K proportion, scaled by background rate
  z      <- qnorm(0.975)
  denom  <- 1 + z^2 / K
  centre <- (p_top + z^2 / (2 * K)) / denom
  halfw  <- (z * sqrt(p_top * (1 - p_top) / K + z^2 / (4 * K^2))) / denom
  ci_lo  <- max(0, (centre - halfw) / p_bg)
  ci_hi  <- (centre + halfw) / p_bg

  # 2x2 test: chi-square or Fisher's exact
  mat <- matrix(c(a, c_bg, K - a, N_bg - c_bg), nrow = 2)
  expected <- suppressWarnings(chisq.test(mat)$expected)
  if (any(expected < 5)) {
    pval      <- fisher.test(mat)$p.value
    test_type <- "Fisher"
  } else {
    pval      <- chisq.test(mat, correct = FALSE)$p.value
    test_type <- "Chi-square"
  }

  data.frame(enrichment = rr, ci_lower = ci_lo, ci_upper = ci_hi,
             p_value = pval, test_used = test_type,
             stringsAsFactors = FALSE)
}

thresholds <- seq(0.5, 0.9, by = 0.05)
all_Signal_enrichment <- list()

for (num_causal in 1:3) {
  for (h2_num in 1:2) {

    out_dir <- paste0(base_dir, "causal_num_", num_causal, "/out/")
    load(paste0(out_dir, "res_simulation_h2_", h2_num, ".RData"))   # res_all

    total_n <- nrow(res_all)
    total_sums_bin <- res_all %>%
      summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))

    enrichment_results <- map_df(thresholds, function(pip_thresh) {
      map_df(method_levels, function(method) {
        pip_col <- pip_cols[[method]]
        sig   <- res_all %>% filter(.data[[pip_col]] > pip_thresh)
        n_sig <- nrow(sig)
        n_bg  <- total_n - n_sig

        if (n_sig == 0) {
          # No signal: return NA for all annotations
          data.frame(annotation = eval_ann, enrichment = NA_real_,
                     ci_lower = NA_real_, ci_upper = NA_real_,
                     p_value = NA_real_, test_used = NA_character_,
                     method = method, threshold = pip_thresh,
                     stringsAsFactors = FALSE)
        } else {
          sig_sums_bin <- sig %>%
            summarise(across(all_of(binary_ann), ~ sum(.x, na.rm = TRUE)))
          bg_sums_bin  <- total_sums_bin - sig_sums_bin

          # Binary annotations: enrichment + CI + test
          bin_stats <- map_df(binary_ann, function(ann) {
            a    <- sig_sums_bin[[ann]]
            c_bg <- bg_sums_bin[[ann]]
            stats <- compute_enrichment_stats(a, n_sig, c_bg, n_bg)
            stats$annotation <- ann
            stats
          })

          # evo2 (continuous): mean-ratio enrichment with bootstrap CI
          evo2_sig_vals <- sig$evo2_score
          evo2_bg_mean  <- (sum(res_all$evo2_score, na.rm = TRUE) -
                              sum(evo2_sig_vals, na.rm = TRUE)) / n_bg
          evo2_rr <- mean(evo2_sig_vals, na.rm = TRUE) / evo2_bg_mean
          # Bootstrap 95% CI for the mean-ratio
          set.seed(42)
          boot_rr <- replicate(2000, {
            bs <- sample(evo2_sig_vals, n_sig, replace = TRUE)
            mean(bs, na.rm = TRUE) / evo2_bg_mean
          })
          evo2_ci <- quantile(boot_rr, c(0.025, 0.975))
          # Wilcoxon rank-sum test: signal vs background evo2 scores
          evo2_bg_vals <- res_all$evo2_score[!(seq_len(total_n) %in%
                            which(res_all[[pip_col]] > pip_thresh))]
          evo2_pval <- wilcox.test(evo2_sig_vals, evo2_bg_vals)$p.value
          evo2_row <- data.frame(
            enrichment = evo2_rr,
            ci_lower = evo2_ci[[1]], ci_upper = evo2_ci[[2]],
            p_value = evo2_pval, test_used = "Wilcoxon",
            annotation = "evo2_score", stringsAsFactors = FALSE)

          bind_rows(bin_stats, evo2_row) %>%
            mutate(method = method, threshold = pip_thresh)
        }
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
      labs(title = paste("Either-PIP enrichment | causal =", num_causal, "h2 =", h2_num),
           x = "PIP threshold", y = "Enrichment (Signal / Background)") +
      theme_bw()
    ggsave(paste0(plot_dir, "Causal_num_", num_causal, "_h2_", h2_num,
                  "_Enrichment_pip_cutoff_either.pdf"),
           p, dpi = 500, width = 200, height = 130, units = 'mm')

    all_Signal_enrichment[[paste0("num_causal_", num_causal, "_h2_", h2_num)]] <- enrichment_results
    cat("done causal", num_causal, "h2", h2_num, "\n")
  }
}

combined_enrichment <- bind_rows(all_Signal_enrichment)
avg_enrichment <- combined_enrichment %>%
  group_by(method, annotation, threshold) %>%
  summarise(avg_enrichment = mean(enrichment, na.rm = TRUE),
            avg_ci_lower   = mean(ci_lower, na.rm = TRUE),
            avg_ci_upper   = mean(ci_upper, na.rm = TRUE),
            .groups = "drop")

avg_plot <- ggplot(avg_enrichment, aes(x = threshold, y = avg_enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  scale_color_manual(values = method_colors) +
  labs(title = "Average Either-PIP enrichment across all causal_num / h2",
       x = "PIP threshold", y = "Average Enrichment (Signal / Background)") +
  theme_bw()
ggsave(paste0(plot_dir, "Summary_Enrichment_pip_cutoff_either.pdf"),
       avg_plot, dpi = 500, width = 200, height = 130, units = 'mm')

save(combined_enrichment, avg_enrichment,
     file = paste0(plot_dir, "enrichment_pip_cutoff_either.RData"))

###############################################################################
#  GROUND-TRUTH ENRICHMENT + MSE  (EITHER case only, matching INFUSE)
#
#  Estimate ground-truth enrichment from the simulation:
#  the fold enrichment of each annotation among the TRUE causal SNPs
#  (Signal != 0 -> the "either"-ancestry causal set), region-restricted,
#  exactly as INFUSE defined cs_flag_true = as.integer(Signal != 0).
#    - binary annotations  -> freq(causal) / freq(non-causal)
#    - evo2 (continuous)   -> mean(causal)  / mean(non-causal)
#  Then MSE = (estimated enrichment at PIP > 0.5  -  truth)^2, per
#  (method, causal, h2, annotation); boxplot across the 6 sim settings.
###############################################################################

calc_true_enrichment <- function(df, ann) {
  df <- df %>% mutate(cs_flag_true = as.integer(Signal != 0))
  df %>%
    group_by(Region) %>%
    filter(sum(cs_flag_true) != 0) %>%
    group_by(cs_flag_true) %>%
    summarise(across(all_of(ann), ~ sum(.x, na.rm = TRUE) / n()), .groups = "drop") %>%
    summarise(across(all_of(ann), ~ .x[cs_flag_true == 1] / .x[cs_flag_true == 0]))
}

True_enrichment_long <- bind_rows(lapply(1:3, function(num_causal) {
  bind_rows(lapply(1:2, function(h2_num) {
    out_dir <- paste0(base_dir, "causal_num_", num_causal, "/out/")
    load(paste0(out_dir, "res_simulation_h2_", h2_num, ".RData"))   # res_all
    calc_true_enrichment(res_all, eval_ann) %>%
      pivot_longer(cols = everything(), names_to = "annotation", values_to = "Prop") %>%
      mutate(causal = num_causal, h2 = h2_num)
  }))
})) %>%
  mutate(annotation = as.character(recode(factor(annotation, levels = ann_levels), !!!ann_recode)))

save(True_enrichment_long, file = paste0(plot_dir, "True_enrichment_either.RData"))

# estimated enrichment at PIP > 0.5 (per method, annotation, causal, h2)
Signal_enrichment_5 <- combined_enrichment %>%
  filter(threshold == 0.5) %>%
  mutate(annotation = as.character(annotation))

merged <- Signal_enrichment_5 %>%
  inner_join(True_enrichment_long, by = c("annotation", "causal", "h2"))

mse_by_annotation <- merged %>%
  group_by(method, causal, h2, annotation) %>%
  summarise(mse = mean((enrichment - Prop)^2, na.rm = TRUE), .groups = "drop") %>%
  mutate(method = factor(method, levels = method_levels),
         annotation = factor(annotation, levels = unname(ann_recode)))

mse_summary <- mse_by_annotation %>%
  group_by(method, annotation) %>%
  summarise(mse = round(mean(mse, na.rm = TRUE), 3), .groups = "drop")

p_mse_box <- ggplot(mse_by_annotation, aes(x = method, y = mse, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = method_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ annotation, scales = "free_y", nrow = 1) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "MSE of Enrichment (vs ground truth)",
       title = "MSE of estimated enrichment vs ground truth (PIP > 0.5, Either)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste0(plot_dir, "MSE_enrichment_either.pdf"),
       p_mse_box, dpi = 500, width = 280, height = 120, units = 'mm')
save(mse_by_annotation, mse_summary, p_mse_box,
     file = paste0(plot_dir, "MSE_enrichment_either.RData"))

###############################################################################
#  Bar plot at PIP = 0.5 with 95% CI + ground-truth enrichment line
###############################################################################
true_avg <- True_enrichment_long %>%
  group_by(annotation) %>%
  summarise(true_enrichment = mean(Prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(annotation = factor(annotation, levels = unname(ann_recode)))

Signal_05 <- avg_enrichment %>% filter(threshold == 0.5) %>%
  mutate(method = factor(method, levels = method_levels))
p_bar <- ggplot(Signal_05, aes(x = method, y = avg_enrichment, fill = method)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = avg_ci_lower, ymax = avg_ci_upper),
                width = 0.3, linewidth = 0.4) +
  geom_hline(data = true_avg, aes(yintercept = true_enrichment),
             linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ annotation, scales = "free_y", nrow = 1) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Fold Enrichment (PIP > 0.5, Either)") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(paste0(plot_dir, "Bar_Enrichment_pip0.5_either.pdf"),
       p_bar, dpi = 500, width = 280, height = 120, units = 'mm')

cat("\n--- Ground-truth enrichment (mean over causal/h2) ---\n")
print(as.data.frame(true_avg))
cat("\n--- Mean MSE by method x annotation ---\n")
print(as.data.frame(mse_summary))
cat("Done (either). Output in ", plot_dir, "\n")
