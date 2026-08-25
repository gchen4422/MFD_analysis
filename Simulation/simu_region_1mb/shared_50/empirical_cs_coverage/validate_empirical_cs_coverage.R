#!/usr/bin/env Rscript

# Truth-based empirical coverage of INDIVIDUAL reported credible sets.
#
# This script deliberately does not use all_Set_data_dataframe/CS_summary as
# the estimator. Those objects collapse all CSs from a method within a region
# into one SNP union and therefore cannot answer per-CS
# question.
#
# Input consistency:
#   * Method outputs are read from the same paths and filenames used by
#     5_result_summarzing_meta.R to create the baseline summary RData.
#   * The simulated causal variants are read from the same PAINTOR input/output
#     table (Signal != 0) and are cross-checked against data_all in that RData.
#   * The analysis includes the same 100 loci x 3 causal-architecture codes x
#     2 heritability codes (600 region-setting combinations).
#
# Estimand:
#   number of individual reported CSs containing >=1 true causal variant
#   --------------------------------------------------------------------
#                 total number of individual reported CSs
#
# The point estimates are raw pooled proportions. No bootstrap or other
# resampling is performed.

suppressPackageStartupMessages({
  library(data.table)
  library(XMAP)
})

args <- commandArgs(trailingOnly = TRUE)
paintor_only <- "--paintor-only" %in% args

simulation_dir <- normalizePath(
  "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50",
  mustWork = TRUE
)
baseline_script <- file.path(simulation_dir, "5_result_summarzing_meta.R")
baseline_rdata <- normalizePath(
  file.path(
    dirname(simulation_dir), "res_summary",
    "shared_50_baseline_updated_meta_updated_xmap.RData"
  ),
  mustWork = TRUE
)
out_dir <- file.path(simulation_dir, "empirical_cs_coverage_results_per_cs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c(
  "MESuSiE", "SuSiE post-hoc", "SuSiE-weighted", "SuSiE-merged",
  "PAINTOR", "SuSiEx", "XMAP", "MultiSuSiE", "CARMA-X"
)

existing_detail_file <- file.path(
  out_dir, "empirical_cs_coverage_per_set.csv"
)
existing_manifest_file <- file.path(
  out_dir, "empirical_cs_coverage_input_manifest.csv"
)

if (paintor_only && !file.exists(existing_detail_file)) {
  stop(
    "--paintor-only requires an existing per-set result at ",
    existing_detail_file
  )
}

as_indices <- function(x, n) {
  if (is.null(x) || length(x) == 0L) return(integer())
  x <- as.integer(round(unlist(x, use.names = FALSE)))
  unique(x[is.finite(x) & x >= 1L & x <= n])
}

collapse_sorted <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  paste(sort(x), collapse = ";")
}

coverage_value <- function(x, i) {
  if (is.null(x) || length(x) < i) return(NA_real_)
  as.numeric(x[[i]])
}

# PAINTOR reports marginal PIPs but no component/CS labels. Construct proxy
# signal-specific CSs using CARMA's credible.set.fun.improved algorithm, adapted
# to multiple ancestry LD matrices. For correlation thresholds 0.50,...,0.99,
# repeatedly examine lead SNPs in descending PIP order. A lead produces a CS
# only when the ORIGINAL (not renormalized) PIPs of its correlated group exceed
# `coverage`; the CS is the smallest descending raw-PIP prefix exceeding that
# value. Remove the selected CS members, then continue. Retain the threshold
# producing the most CSs, choosing the highest threshold on ties.
paintor_ld_credible_sets <- function(
    pip, ld_matrices, coverage = 0.95,
    correlation_thresholds = c(
      seq(0.50, 0.95, 0.05), seq(0.96, 0.99, 0.01)
    )) {
  pip <- as.numeric(pip)
  n <- length(pip)
  if (
    length(ld_matrices) == 0L ||
    any(vapply(
      ld_matrices,
      function(ld) !is.matrix(ld) || !identical(dim(ld), c(n, n)),
      logical(1)
    ))
  ) {
    stop("Every PAINTOR LD matrix must be square and aligned to the PIPs.")
  }
  if (
    !is.finite(coverage) || coverage <= 0 || coverage > 1 ||
    length(correlation_thresholds) == 0L ||
    any(!is.finite(correlation_thresholds)) ||
    any(correlation_thresholds < 0 | correlation_thresholds > 1)
  ) {
    stop("Invalid PAINTOR CS coverage or LD-correlation thresholds.")
  }

  pip[!is.finite(pip) | pip < 0] <- 0
  threshold_sets <- vector("list", length(correlation_thresholds))

  for (threshold_index in seq_along(correlation_thresholds)) {
    threshold <- correlation_thresholds[threshold_index]
    remaining <- which(pip > 0)
    pip_order <- remaining[order(pip[remaining], decreasing = TRUE)]
    lead_position <- 1L
    sets <- list()

    while (lead_position <= length(pip_order)) {
      tail <- pip_order[seq.int(lead_position, length(pip_order))]
      if (sum(pip[tail]) <= coverage) break

      lead <- pip_order[lead_position]
      correlated <- Reduce(
        `|`,
        lapply(
          ld_matrices,
          function(ld) {
            r <- abs(as.numeric(ld[lead, remaining]))
            is.finite(r) & r > threshold
          }
        )
      )
      # Retain the lead defensively if an LD diagonal is malformed.
      cluster <- unique(c(lead, remaining[correlated]))
      cluster_mass <- sum(pip[cluster])

      if (cluster_mass > coverage) {
        ord <- order(pip[cluster], decreasing = TRUE)
        endpoint <- which(cumsum(pip[cluster][ord]) > coverage)[1]
        cs_indices <- cluster[ord[seq_len(endpoint)]]
        sets[[length(sets) + 1L]] <- list(
          indices = cs_indices,
          lead_index = lead,
          cluster_indices = cluster,
          cluster_mass = cluster_mass,
          fitted_mass = sum(pip[cs_indices]),
          correlation_threshold = threshold
        )
        remaining <- setdiff(remaining, cs_indices)
        pip_order <- pip_order[!pip_order %in% cs_indices]
      } else {
        lead_position <- lead_position + 1L
      }
    }
    threshold_sets[[threshold_index]] <- sets
  }

  set_counts <- lengths(threshold_sets)
  if (sum(set_counts) == 0L) return(list())
  selected <- max(which(set_counts == max(set_counts)))
  threshold_sets[[selected]]
}

if (paintor_only) {
  old_detail <- fread(existing_detail_file)
  rows <- list(old_detail[Method != "PAINTOR"])
  row_counter <- 1L
  manifest <- if (file.exists(existing_manifest_file)) {
    list(fread(existing_manifest_file))
  } else {
    list()
  }
  manifest_counter <- length(manifest)
} else {
  rows <- list()
  row_counter <- 0L
  manifest <- list()
  manifest_counter <- 0L
}
truth_rows <- list()
truth_counter <- 0L

record_file <- function(region_id, input_type, path) {
  info <- file.info(path)
  manifest_counter <<- manifest_counter + 1L
  manifest[[manifest_counter]] <<- data.table(
    Region_ID = region_id,
    Input_type = input_type,
    Path = normalizePath(path, mustWork = TRUE),
    Bytes = as.numeric(info$size),
    Modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %z")
  )
}

add_cs <- function(
    method, causal_code, h2_code, locus, region_id, source, cs_id,
    indices, paintor, fitted_mass = NA_real_, lead_index = NA_integer_,
    signal_cluster_size = NA_integer_, signal_posterior_mass = NA_real_,
    ld_correlation_threshold = NA_real_) {
  indices <- as_indices(indices, nrow(paintor))
  if (length(indices) == 0L) return(invisible(NULL))
  causal_indices <- which(paintor$Signal != 0)
  hit <- intersect(indices, causal_indices)
  row_counter <<- row_counter + 1L
  rows[[row_counter]] <<- data.table(
    Method = method,
    Causal_code = causal_code,
    N_causal_total = length(causal_indices),
    H2_code = h2_code,
    Locus = locus,
    Region_ID = region_id,
    Source = source,
    CS_ID = as.character(cs_id),
    CS_size = length(indices),
    Covered = length(hit) > 0L,
    N_causal_hit = length(hit),
    Causal_SNPs = collapse_sorted(paintor$SNP[causal_indices]),
    Hit_causal_SNPs = collapse_sorted(paintor$SNP[hit]),
    CS_SNPs = collapse_sorted(paintor$SNP[indices]),
    Fitted_posterior_mass = as.numeric(fitted_mass),
    Lead_SNP = if (
      length(lead_index) == 1L && is.finite(lead_index) &&
      lead_index >= 1L && lead_index <= nrow(paintor)
    ) {
      as.character(paintor$SNP[lead_index])
    } else {
      NA_character_
    },
    Signal_cluster_size = as.integer(signal_cluster_size),
    Signal_posterior_mass = as.numeric(signal_posterior_mass),
    LD_correlation_threshold = as.numeric(ld_correlation_threshold)
  )
  invisible(NULL)
}

for (causal_code in 1:3) {
  result_dir <- file.path(
    simulation_dir, paste0("causal_num_", causal_code), "result"
  )
  summary_dir <- file.path(
    simulation_dir, paste0("causal_num_", causal_code), "summary_data"
  )

  for (h2_code in 1:2) {
    for (locus in 1:100) {
      suffix <- sprintf(
        "CAUSAL_%d_LOCI_%d_h2_%d", causal_code, locus, h2_code
      )
      region_id <- sprintf(
        "c%d_h%d_l%03d", causal_code, h2_code, locus
      )

      main_file <- file.path(
        result_dir, paste0("MESuSiE_", suffix, ".RData")
      )
      metal_file <- file.path(
        result_dir, paste0("Metal_SuSiE_", suffix, "_updated.RData")
      )
      xmap_file <- file.path(
        result_dir, paste0("XMAP_", suffix, ".RData")
      )
      paintor_file <- file.path(
        result_dir, paste0(suffix, ".mcmc.paintor")
      )
      ld1_file <- file.path(summary_dir, paste0(suffix, ".LD1"))
      ld2_file <- file.path(summary_dir, paste0(suffix, ".LD2"))
      susiex_file <- file.path(
        simulation_dir, "susiex_result",
        paste0("SuSiEx_", suffix, "_output_cs95.cs")
      )
      multisusie_cs_file <- file.path(
        simulation_dir, "multi_susie_result",
        paste0("MultiSuSiE_", suffix, "_output_cs.txt")
      )
      multisusie_coverage_file <- file.path(
        simulation_dir, "multi_susie_result",
        paste0("MultiSuSiE_", suffix, "_output_coverage.txt")
      )

      required <- if (paintor_only) {
        c(paintor_file, ld1_file, ld2_file)
      } else {
        c(
          main_file, metal_file, xmap_file, paintor_file,
          ld1_file, ld2_file, susiex_file, multisusie_cs_file
        )
      }
      missing <- required[!file.exists(required)]
      if (length(missing) > 0L) {
        stop(
          "Missing input(s) for ", region_id, ":\n",
          paste(missing, collapse = "\n")
        )
      }

      if (!paintor_only) {
        record_file(region_id, "MESuSiE/SuSiE", main_file)
        record_file(region_id, "SuSiE meta-analysis", metal_file)
        record_file(region_id, "XMAP", xmap_file)
        record_file(region_id, "PAINTOR/truth", paintor_file)
        record_file(region_id, "PAINTOR/XMAP LD1", ld1_file)
        record_file(region_id, "PAINTOR/XMAP LD2", ld2_file)
        record_file(region_id, "SuSiEx", susiex_file)
        record_file(region_id, "MultiSuSiE CS", multisusie_cs_file)
        if (file.exists(multisusie_coverage_file)) {
          record_file(
            region_id, "MultiSuSiE coverage", multisusie_coverage_file
          )
        }
      }

      paintor <- fread(paintor_file)
      if ("RSID" %in% names(paintor) && !"SNP" %in% names(paintor)) {
        setnames(paintor, "RSID", "SNP")
      }
      required_paintor <- c("SNP", "POS", "Signal", "Posterior_Prob")
      if (!all(required_paintor %in% names(paintor))) {
        stop("Unexpected PAINTOR columns in ", paintor_file)
      }
      causal_indices <- which(paintor$Signal != 0)
      truth_counter <- truth_counter + 1L
      truth_rows[[truth_counter]] <- data.table(
        Causal_code = causal_code,
        H2_code = h2_code,
        Locus = locus,
        Region_ID = region_id,
        Causal_SNPs = collapse_sorted(paintor$SNP[causal_indices]),
        N_causal_total = length(causal_indices)
      )

      pip <- as.numeric(paintor$Posterior_Prob)
      positive_pip_indices <- which(is.finite(pip) & pip > 0)
      if (paintor_only) {
        # PAINTOR CS membership can only involve positive-PIP SNPs. Reading
        # their LD submatrix avoids materializing two full dense matrices in
        # the fast refresh path.
        if (length(positive_pip_indices) > 0L) {
          ld1 <- as.matrix(
            fread(ld1_file, select = positive_pip_indices)
          )[positive_pip_indices, , drop = FALSE]
          ld2 <- as.matrix(
            fread(ld2_file, select = positive_pip_indices)
          )[positive_pip_indices, , drop = FALSE]
        } else {
          ld1 <- ld2 <- matrix(numeric(), nrow = 0L, ncol = 0L)
        }
        paintor_index_map <- positive_pip_indices
      } else {
        ld1 <- as.matrix(fread(ld1_file))
        ld2 <- as.matrix(fread(ld2_file))
        paintor_index_map <- seq_along(pip)
      }

      if (!paintor_only) {
        main_env <- new.env(parent = emptyenv())
        load(main_file, envir = main_env)
        metal_env <- new.env(parent = emptyenv())
        load(metal_file, envir = metal_env)
        xmap_env <- new.env(parent = emptyenv())
        load(xmap_file, envir = xmap_env)

      # MESuSiE: retain each reported component CS separately.
      mesusie_cs <- main_env$MESuSiE_res$cs$cs
      mesusie_mass <- main_env$MESuSiE_res$cs$coverage
      if (length(mesusie_cs) > 0L) {
        for (i in seq_along(mesusie_cs)) {
          id <- names(mesusie_cs)[i]
          if (is.null(id) || is.na(id) || !nzchar(id)) id <- paste0("L", i)
          add_cs(
            "MESuSiE", causal_code, h2_code, locus, region_id, "joint",
            id,
            mesusie_cs[[i]], paintor, coverage_value(mesusie_mass, i)
          )
        }
      }

      # SuSiE post-hoc: EUR and AFR CSs are distinct reported sets.
      for (ancestry in c("EUR", "AFR")) {
        fit <- if (ancestry == "EUR") main_env$susie_EU else main_env$susie_BB
        cs_list <- fit$sets$cs
        cs_mass <- fit$sets$coverage
        if (length(cs_list) > 0L) {
          for (i in seq_along(cs_list)) {
            id <- names(cs_list)[i]
            if (is.null(id) || is.na(id) || !nzchar(id)) id <- paste0("L", i)
            add_cs(
              "SuSiE post-hoc", causal_code, h2_code, locus, region_id,
              ancestry, id, cs_list[[i]], paintor,
              coverage_value(cs_mass, i)
            )
          }
        }
      }

      # Meta-analysis SuSiE variants.
      for (spec in list(
        list("SuSiE-weighted", metal_env$susie_weighted),
        list("SuSiE-merged", metal_env$susie_merged)
      )) {
        method <- spec[[1]]
        fit <- spec[[2]]
        cs_list <- fit$sets$cs
        cs_mass <- fit$sets$coverage
        if (length(cs_list) > 0L) {
          for (i in seq_along(cs_list)) {
            id <- names(cs_list)[i]
            if (is.null(id) || is.na(id) || !nzchar(id)) id <- paste0("L", i)
            add_cs(
              method, causal_code, h2_code, locus, region_id,
              "meta-analysis", id, cs_list[[i]], paintor,
              coverage_value(cs_mass, i)
            )
          }
        }
      }

      }

      # PAINTOR has no component/CS label. Apply CARMA's adaptive correlation-
      # threshold construction to its raw marginal PIPs. A pair is correlated
      # if it clears the threshold in either ancestry-specific LD matrix.
      paintor_sets <- paintor_ld_credible_sets(
        pip[paintor_index_map], list(EUR = ld1, AFR = ld2),
        coverage = 0.95
      )
      if (length(paintor_sets) > 0L) {
        for (i in seq_along(paintor_sets)) {
          cs <- paintor_sets[[i]]
          cs$indices <- paintor_index_map[cs$indices]
          cs$lead_index <- paintor_index_map[cs$lead_index]
          cs$cluster_indices <- paintor_index_map[cs$cluster_indices]
          add_cs(
            "PAINTOR", causal_code, h2_code, locus, region_id,
            "CARMA-style LD proxy signal; EUR or AFR",
            paste0("CS", i), cs$indices, paintor,
            fitted_mass = cs$fitted_mass,
            lead_index = cs$lead_index,
            signal_cluster_size = length(cs$cluster_indices),
            signal_posterior_mass = cs$cluster_mass,
            ld_correlation_threshold = cs$correlation_threshold
          )
        }
      }

      if (!paintor_only) {
      # SuSiEx: CS_ID in the same .cs file used by the baseline script.
      sx <- fread(susiex_file)
      if (nrow(sx) > 0L) {
        for (id in unique(sx$CS_ID)) {
          part <- sx[CS_ID == id]
          idx <- which(paintor$POS %in% part$BP)
          add_cs(
            "SuSiEx", causal_code, h2_code, locus, region_id, "joint",
            id, idx, paintor, sum(part$CS_PIP)
          )
        }
      }

      # XMAP: same construction as 5_result_summarzing_meta.R. A CS must pass
      # purity in both ancestry LD matrices; EUR membership is retained.
      # Use every variant when calculating purity. The get_CS() default
      # randomly samples 100 variants from larger CSs, which makes the
      # reported-CS denominator vary across otherwise identical runs.
      xmap_cs1 <- get_CS(
        xmap_env$xmap, Xcorr = ld1, coverage = 0.95, min_abs_corr = 0.5,
        n_purity = Inf
      )
      xmap_cs2 <- get_CS(
        xmap_env$xmap, Xcorr = ld2, coverage = 0.95, min_abs_corr = 0.5,
        n_purity = Inf
      )
      shared_ids <- intersect(names(xmap_cs1$cs), names(xmap_cs2$cs))
      if (length(shared_ids) > 0L) {
        for (id in shared_ids) {
          i <- match(id, names(xmap_cs1$cs))
          add_cs(
            "XMAP", causal_code, h2_code, locus, region_id,
            "purity in EUR+AFR; EUR membership", id,
            xmap_cs1$cs[[id]], paintor,
            coverage_value(xmap_cs1$coverage, i)
          )
        }
      }

      # MultiSuSiE: one observation for every positive CS label.
      ms <- fread(multisusie_cs_file)
      ms_mass <- if (file.exists(multisusie_coverage_file)) {
        fread(multisusie_coverage_file)
      } else {
        data.table()
      }
      if (nrow(ms) > 0L && "CS" %in% names(ms)) {
        for (id in sort(unique(ms[CS > 0]$CS))) {
          part <- ms[CS == id]
          idx <- which(paintor$POS %in% part$POS)
          mass <- if (
            nrow(ms_mass) > 0L &&
            all(c("CS", "coverage") %in% names(ms_mass)) &&
            id %in% ms_mass$CS
          ) {
            ms_mass[CS == id]$coverage[1]
          } else {
            NA_real_
          }
          add_cs(
            "MultiSuSiE", causal_code, h2_code, locus, region_id,
            "joint", id, idx, paintor, mass
          )
        }
      }

      # CARMA-X reports ancestry-specific CSs; count each reported set.
      carmax_results <- main_env$carmax_results
      if (length(carmax_results) >= 2L) {
        for (a in 1:2) {
          source <- c("EUR", "AFR")[a]
          cs_list <- carmax_results[[a]][["Credible set"]][[2]]
          pips <- carmax_results[[a]]$PIPs
          if (length(cs_list) > 0L) {
            for (i in seq_along(cs_list)) {
              idx <- as_indices(cs_list[[i]], nrow(paintor))
              add_cs(
                "CARMA-X", causal_code, h2_code, locus, region_id,
                source, i, idx, paintor, sum(pips[idx], na.rm = TRUE)
              )
            }
          }
        }
      }
      }
    }
    message("Completed causal_code=", causal_code, ", h2_code=", h2_code)
  }
}

coverage_detail <- rbindlist(rows, use.names = TRUE, fill = TRUE)
input_manifest <- rbindlist(manifest, use.names = TRUE)
truth_detail <- rbindlist(truth_rows, use.names = TRUE)

# Older manifests called these matrices XMAP inputs. They are shared by XMAP
# and the new PAINTOR LD-proxy CS construction.
input_manifest[Input_type == "XMAP LD1", Input_type := "PAINTOR/XMAP LD1"]
input_manifest[Input_type == "XMAP LD2", Input_type := "PAINTOR/XMAP LD2"]

if (nrow(truth_detail) != 600L || uniqueN(truth_detail$Region_ID) != 600L) {
  stop("Truth-design validation failed: expected exactly 600 regions.")
}
duplicate_key <- coverage_detail[
  , .N,
  by = .(Method, Causal_code, H2_code, Locus, Source, CS_ID)
][N > 1L]
if (nrow(duplicate_key) > 0L) {
  stop("Duplicate individual-CS keys detected.")
}
if (!setequal(unique(coverage_detail$Method), method_order)) {
  stop("Method validation failed.")
}
if (any(coverage_detail$CS_size <= 0L)) {
  stop("Zero-size rows should not appear in the individual-CS denominator.")
}
paintor_detail <- coverage_detail[Method == "PAINTOR"]
if (
  any(!is.finite(paintor_detail$Fitted_posterior_mass)) ||
  any(paintor_detail$Fitted_posterior_mass <= 0.95) ||
  any(is.na(paintor_detail$Lead_SNP)) ||
  any(paintor_detail$Signal_cluster_size < paintor_detail$CS_size) ||
  any(
    !is.finite(paintor_detail$Signal_posterior_mass) |
    paintor_detail$Signal_posterior_mass <= 0.95
  ) ||
  any(!is.finite(paintor_detail$LD_correlation_threshold)) ||
  any(
    paintor_detail$LD_correlation_threshold < 0.50 |
    paintor_detail$LD_correlation_threshold > 0.99
  )
) {
  stop("PAINTOR LD-proxy credible-set validation failed.")
}
paintor_membership <- paintor_detail[
  , .(SNP = unlist(strsplit(CS_SNPs, ";", fixed = TRUE))),
  by = Region_ID
]
if (anyDuplicated(paintor_membership, by = c("Region_ID", "SNP"))) {
  stop("A PAINTOR SNP was assigned to more than one proxy credible set.")
}

# Confirm that the PAINTOR-derived truth is exactly the truth saved in the
# baseline RData. This is validation only; CS_summary is not
# used to calculate empirical per-CS coverage.
baseline_env <- new.env(parent = emptyenv())
load(baseline_rdata, envir = baseline_env)
if (!exists("data_all", envir = baseline_env, inherits = FALSE)) {
  stop("Baseline RData does not contain data_all.")
}
baseline_truth <- as.data.table(baseline_env$data_all)[
  ,
  .(
    Baseline_Causal_SNPs = collapse_sorted(SNP[Signal != 0]),
    Baseline_N_causal = uniqueN(SNP[Signal != 0])
  ),
  by = .(
    Causal_code = as.integer(causal_num),
    H2_code = as.integer(h2),
    Locus = as.integer(locus)
  )
]
truth_audit <- merge(
  truth_detail, baseline_truth,
  by = c("Causal_code", "H2_code", "Locus"),
  all = TRUE
)
truth_audit[
  ,
  Exact_truth_match :=
    Causal_SNPs == Baseline_Causal_SNPs &
    N_causal_total == Baseline_N_causal
]
if (nrow(truth_audit) != 600L || any(!truth_audit$Exact_truth_match)) {
  stop("Per-CS causal truth does not match baseline data_all.")
}

coverage_detail[, Method := factor(Method, levels = method_order)]
coverage_overall <- coverage_detail[
  ,
  .(
    N_regions = 600L,
    N_regions_with_CS = uniqueN(Region_ID),
    N_CS = .N,
    N_covered = sum(Covered),
    Empirical_coverage = mean(Covered),
    Mean_CS_size = mean(CS_size),
    Median_CS_size = as.numeric(median(CS_size))
  ),
  by = Method
]
coverage_by_setting <- coverage_detail[
  ,
  .(
    N_regions_with_CS = uniqueN(Region_ID),
    N_CS = .N,
    N_covered = sum(Covered),
    Empirical_coverage = mean(Covered),
    Mean_CS_size = mean(CS_size),
    Median_CS_size = as.numeric(median(CS_size))
  ),
  by = .(Method, Causal_code, H2_code)
]
setorder(coverage_detail, Method, Causal_code, H2_code, Locus, Source, CS_ID)
setorder(coverage_overall, Method)
setorder(coverage_by_setting, Method, Causal_code, H2_code)
coverage_detail[, Method := as.character(Method)]
coverage_overall[, Method := as.character(Method)]
coverage_by_setting[, Method := as.character(Method)]

analysis_metadata <- data.table(
  Field = c(
    "Estimand", "Point_estimator", "Resampling", "Baseline_script",
    "Baseline_RData", "Region_settings", "Truth_match",
    "PAINTOR_CS_construction", "PAINTOR_LD_rule", "PAINTOR_refresh_mode"
  ),
  Value = c(
    paste(
      "Among all individually reported credible sets, the proportion",
      "containing at least one simulated causal variant"
    ),
    "sum(Covered) / number of reported individual credible sets",
    "None",
    normalizePath(baseline_script, mustWork = TRUE),
    baseline_rdata,
    "100 loci x 3 causal codes x 2 h2 codes = 600",
    "Exact across all 600 region-setting combinations",
    paste(
      "CARMA credible.set.fun.improved adaptation: scan correlation",
      "thresholds 0.50-0.99; require correlated raw PIP mass >0.95;",
      "retain the smallest raw-PIP prefix >0.95; remove selected CS SNPs;",
      "choose the threshold yielding most CSs (highest threshold on ties)"
    ),
    "At each candidate threshold, SNP joins the lead cluster when abs(r) exceeds the threshold in EUR or AFR LD",
    if (paintor_only) {
      "Only PAINTOR rows recomputed; existing rows for all other methods retained"
    } else {
      "All methods recomputed"
    }
  )
)

fwrite(
  coverage_detail,
  file.path(out_dir, "empirical_cs_coverage_per_set.csv")
)
fwrite(
  coverage_overall,
  file.path(out_dir, "empirical_cs_coverage_overall.csv")
)
fwrite(
  coverage_by_setting,
  file.path(out_dir, "empirical_cs_coverage_by_setting.csv")
)
fwrite(
  input_manifest,
  file.path(out_dir, "empirical_cs_coverage_input_manifest.csv")
)
fwrite(
  truth_audit,
  file.path(out_dir, "empirical_cs_coverage_truth_audit.csv")
)
fwrite(
  analysis_metadata,
  file.path(out_dir, "empirical_cs_coverage_metadata.csv")
)
save(
  coverage_detail, coverage_overall, coverage_by_setting, input_manifest,
  truth_audit, analysis_metadata,
  file = file.path(out_dir, "empirical_cs_coverage_per_set.RData")
)

cat("\nEmpirical coverage of individual reported credible sets\n")
cat("No bootstrap or other resampling was used.\n")
cat("Causal truth matched baseline data_all in 600/600 regions.\n\n")
print(
  coverage_overall[
    ,
    .(
      Method,
      Covered_over_reported = sprintf("%s / %s", N_covered, N_CS),
      Empirical_coverage = sprintf("%.3f", Empirical_coverage)
    )
  ],
  row.names = FALSE
)
