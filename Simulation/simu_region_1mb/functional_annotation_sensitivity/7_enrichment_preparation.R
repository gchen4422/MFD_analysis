###############################################################################
#
#  6_enrichment_preparation.R   (MFD revision / baseline_functional_annotation)
#
#  Adapted from:
#  shared_50_mesusie_inf/6_enrichment_prepration.R
#
#  ---------------------------------------------------------------------------
#  ANNOTATION REUSE NOTE
#  ---------------------------------------------------------------------------
#  This simulation uses the SAME genomic regions / SNPs as the INFUSE run
#  (shared_50_mesusie_inf). The 5 baseline BINARY annotations
#     H3K27ac_Hnisz_common, DHS_Trynka_common, H3K4me3_Trynka_common,
#     Conserved_LindbladToh_common, non_synonymous_common
#  come from the identical source files:
#     risk_loci_annot/loci_<LD_BLOCK>.annot
#  and join 1:1 onto every SNP (verified: 0 mismatches, 0 NA). They are
#  therefore IDENTICAL to INFUSE and do NOT need to be regenerated.
#
#  The one difference vs INFUSE is the evo2 score (evo2_score), which the new
#  pipeline uses (it feeds Paintor_all_fun) but INFUSE did not. evo2 is a
#  rank-normalized [0,1] continuous score:
#     evo2_score = rank(-evo2_delta_score) / n_snp
#  built per locus from
#     risk_loci_annot/CAUSAL_1_LOCI_<LD_BLOCK>_h2_1_eur_scored.csv
#  (this matches annot/Region_<LD_BLOCK>.annot$evo2_score to ~1e-16). evo2 is a
#  per-SNP genomic property, identical across causal_num / h2.
#
#  Rather than re-parse every method's raw output, reuse the already
#  summarized per-SNP result table `data_all` produced by
#  5_result_summarzing_functional.R and simply attach the annotations to it,
#  then split into per-(causal_num, h2) `res_all` files for the downstream
#  enrichment scripts (7, 9_1, 9_2, 9_3).
###############################################################################

suppressMessages({
  library(data.table)
  library(dplyr)
})

base_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/"
annot_src <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/"
res_summary_file <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/functional_annotation_sensitivity_w2.RData"

binary_ann <- c("H3K27ac_Hnisz_common", "DHS_Trynka_common", "H3K4me3_Trynka_common",
                "Conserved_LindbladToh_common", "non_synonymous_common")

##############################################################################
# 1) Build the per-SNP annotation table (RSID = CHR:BP) once for all 100 loci
#    - 5 binary annotations from loci_<i>.annot
#    - evo2_score from CAUSAL_1_LOCI_<i>_h2_1_eur_scored.csv (genomic, h2/causal
#      invariant), rank-normalized exactly as the fine-mapping pipeline does.
##############################################################################
cat("Building annotation table (5 baseline binary + evo2) ...\n")
annot_all <- rbindlist(lapply(1:100, function(i) {
  bin <- fread(paste0(annot_src, "loci_", i, ".annot"),
               select = c("CHR", "BP", binary_ann))
  bin[, RSID := paste0(CHR, ":", BP)]

  sc <- fread(paste0(annot_src, "CAUSAL_1_LOCI_", i, "_h2_1_eur_scored.csv"))
  # row-aligned with loci_<i>.annot; key by RSID to be safe
  evo2_dt <- data.table(RSID = sc$RSID,
                        evo2_score = rank(-sc$evo2_delta_score) / nrow(sc))

  out <- merge(bin[, c("RSID", binary_ann), with = FALSE], evo2_dt,
               by = "RSID", all.x = TRUE, sort = FALSE)
  out
}))
cat("  annotation rows:", nrow(annot_all),
    " unique RSID:", length(unique(annot_all$RSID)), "\n")

##############################################################################
# 2) Load the summarized per-SNP results (all methods) and attach annotations
##############################################################################
cat("Loading data_all from result summary ...\n")
load(res_summary_file)   # provides data_all (among others)

stopifnot("SNP" %in% colnames(data_all))

##############################################################################
# 3) Split into per-(causal_num, h2) res_all and save for downstream scripts
##############################################################################
for (num_causal in 1:3) {
  for (h2_num in 1:2) {

    out_dir <- paste0(base_dir, "causal_num_", num_causal, "/out/")
    system(paste0("mkdir -p ", out_dir))

    res_all <- data_all %>%
      filter(causal_num == num_causal, h2 == h2_num) %>%
      dplyr::rename(Region = locus) %>%
      left_join(annot_all, by = c("SNP" = "RSID"))

    # sanity check on annotation coverage
    n_na <- sum(is.na(res_all$H3K27ac_Hnisz_common))
    cat("causal", num_causal, "h2", h2_num,
        "| SNPs:", nrow(res_all), "| regions:", dplyr::n_distinct(res_all$Region),
        "| annotation NA:", n_na, "\n")

    res_name <- paste0(out_dir, "res_simulation_h2_", h2_num, ".RData")
    save(res_all, file = res_name)
  }
}

cat("Done. res_simulation_h2_<h2>.RData written to each causal_num_*/out/.\n")
