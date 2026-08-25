#!/usr/bin/env Rscript

# Generate the four baseline coverage artifacts from truth-based coverage of
# INDIVIDUAL reported credible sets.
#
# Empirical coverage is:
#
#   number of reported CSs containing >=1 simulated causal variant
#   ----------------------------------------------------------------
#                total number of reported individual CSs
#
# This script does not use fitted posterior mass or the region-level union in
# CS_summary/all_Set_data_dataframe. No bootstrap or other resampling is used.
#
# By default, the script uses the most recent verified per-CS extraction from
# 5_validate_empirical_cs_coverage.R. Pass --refresh to reconstruct that table
# from the raw method outputs before creating the compatibility artifacts.
#
# The output directory can be overridden with COVERAGE_REVISION_OUT_DIR. This
# is useful for staging and validation before replacing published artifacts.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
refresh <- "--refresh" %in% args

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_arg),
    mustWork = TRUE
  )
  simulation_dir <- dirname(script_path)
} else {
  simulation_dir <- normalizePath(getwd(), mustWork = TRUE)
}

detail_file <- file.path(
  simulation_dir,
  "empirical_cs_coverage_results_per_cs",
  "empirical_cs_coverage_per_set.csv"
)
truth_audit_file <- file.path(
  simulation_dir,
  "empirical_cs_coverage_results_per_cs",
  "empirical_cs_coverage_truth_audit.csv"
)
validator_script <- file.path(
  simulation_dir,
  "5_validate_empirical_cs_coverage.R"
)

if (refresh || !file.exists(detail_file) || !file.exists(truth_audit_file)) {
  message("Refreshing empirical per-CS coverage from raw method outputs...")
  sys.source(validator_script, envir = new.env(parent = globalenv()))
}

if (!file.exists(detail_file) || !file.exists(truth_audit_file)) {
  stop("Verified empirical per-CS results were not generated.")
}

coverage_detail <- fread(detail_file)
truth_audit <- fread(truth_audit_file)

required_detail <- c(
  "Method", "Causal_code", "H2_code", "Locus", "Region_ID", "Source",
  "CS_ID", "CS_size", "Covered"
)
missing_detail <- setdiff(required_detail, names(coverage_detail))
if (length(missing_detail) > 0L) {
  stop(
    "Missing required per-CS columns: ",
    paste(missing_detail, collapse = ", ")
  )
}
if (
  nrow(truth_audit) != 600L ||
  !"Exact_truth_match" %in% names(truth_audit) ||
  any(!as.logical(truth_audit$Exact_truth_match))
) {
  stop("Truth audit failed; expected exact agreement for all 600 settings.")
}

coverage_detail[, `:=`(
  Causal_code = as.integer(Causal_code),
  H2_code = as.integer(H2_code),
  Locus = as.integer(Locus),
  CS_size = as.integer(CS_size),
  Covered = as.logical(Covered)
)]

method_order_detail <- c(
  "MESuSiE", "SuSiE post-hoc", "SuSiE-weighted", "SuSiE-merged",
  "PAINTOR", "SuSiEx", "XMAP", "MultiSuSiE", "CARMA-X"
)
if (!setequal(unique(coverage_detail$Method), method_order_detail)) {
  stop("Unexpected methods in empirical per-CS results.")
}
if (anyNA(coverage_detail$Covered)) {
  stop("Covered contains missing values.")
}

duplicate_keys <- coverage_detail[
  ,
  .N,
  by = .(Method, Causal_code, H2_code, Locus, Source, CS_ID)
][N > 1L]
if (nrow(duplicate_keys) > 0L) {
  stop("Duplicate individual credible-set identifiers were detected.")
}

region_grid <- CJ(
  Num_Causal = 1:3,
  H2 = 1:2,
  Locus = 1:100,
  sorted = TRUE
)
setcolorder(region_grid, c("Locus", "H2", "Num_Causal"))
setorder(region_grid, Num_Causal, H2, Locus)

build_method_columns <- function(
    method, n_name, coverage_name, source = NULL) {
  x <- coverage_detail[Method == method]
  if (!is.null(source)) {
    x <- x[Source == source]
  }
  x <- x[
    order(Causal_code, H2_code, Locus, Source, CS_ID),
    .(
      n_cs = .N,
      per_cs_coverage = paste(as.integer(Covered), collapse = ";")
    ),
    by = .(
      Num_Causal = Causal_code,
      H2 = H2_code,
      Locus
    )
  ]
  out <- merge(
    copy(region_grid), x,
    by = c("Num_Causal", "H2", "Locus"),
    all.x = TRUE,
    sort = FALSE
  )
  out[is.na(n_cs), n_cs := 0L]
  setnames(
    out,
    c("n_cs", "per_cs_coverage"),
    c(n_name, coverage_name)
  )
  out
}

coverage_df <- copy(region_grid)
method_specs <- list(
  list(
    "MESuSiE", "MESuSiE_n_cs", "MESuSiE_per_cs_coverage", NULL
  ),
  list(
    "SuSiE post-hoc", "SuSiE_EUR_n_cs",
    "SuSiE_EUR_per_cs_coverage", "EUR"
  ),
  list(
    "SuSiE post-hoc", "SuSiE_AFR_n_cs",
    "SuSiE_AFR_per_cs_coverage", "AFR"
  ),
  list(
    "SuSiE-weighted", "SuSiE_weighted_n_cs",
    "SuSiE_weighted_per_cs_coverage", NULL
  ),
  list(
    "SuSiE-merged", "SuSiE_merged_n_cs",
    "SuSiE_merged_per_cs_coverage", NULL
  ),
  list(
    "PAINTOR", "Paintor_n_cs", "Paintor_per_cs_coverage", NULL
  ),
  list(
    "SuSiEx", "SuSiEx_n_cs", "SuSiEx_per_cs_coverage", NULL
  ),
  list(
    "XMAP", "XMAP_n_cs", "XMAP_per_cs_coverage", NULL
  ),
  list(
    "MultiSuSiE", "MultiSuSiE_n_cs",
    "MultiSuSiE_per_cs_coverage", NULL
  ),
  list(
    "CARMA-X", "CARMAX_n_cs", "CARMAX_per_cs_coverage", NULL
  )
)

for (spec in method_specs) {
  method_table <- build_method_columns(
    method = spec[[1]],
    n_name = spec[[2]],
    coverage_name = spec[[3]],
    source = spec[[4]]
  )
  new_columns <- c(spec[[2]], spec[[3]])
  coverage_df <- merge(
    coverage_df,
    method_table[
      ,
      c("Num_Causal", "H2", "Locus", new_columns),
      with = FALSE
    ],
    by = c("Num_Causal", "H2", "Locus"),
    all.x = TRUE,
    sort = FALSE
  )
}

paintor_sizes <- coverage_detail[
  Method == "PAINTOR"
][
  order(Causal_code, H2_code, Locus, Source, CS_ID),
  .(Paintor_per_cs_size = paste(as.integer(CS_size), collapse = ";")),
  by = .(
    Num_Causal = Causal_code,
    H2 = H2_code,
    Locus
  )
]
coverage_df <- merge(
  coverage_df, paintor_sizes,
  by = c("Num_Causal", "H2", "Locus"),
  all.x = TRUE,
  sort = FALSE
)

output_column_order <- c(
  "Locus", "H2", "Num_Causal",
  "MESuSiE_n_cs", "MESuSiE_per_cs_coverage",
  "SuSiE_EUR_n_cs", "SuSiE_EUR_per_cs_coverage",
  "SuSiE_AFR_n_cs", "SuSiE_AFR_per_cs_coverage",
  "SuSiE_weighted_n_cs", "SuSiE_weighted_per_cs_coverage",
  "SuSiE_merged_n_cs", "SuSiE_merged_per_cs_coverage",
  "Paintor_n_cs", "Paintor_per_cs_coverage", "Paintor_per_cs_size",
  "SuSiEx_n_cs", "SuSiEx_per_cs_coverage",
  "XMAP_n_cs", "XMAP_per_cs_coverage",
  "MultiSuSiE_n_cs", "MultiSuSiE_per_cs_coverage",
  "CARMAX_n_cs", "CARMAX_per_cs_coverage"
)
setcolorder(coverage_df, output_column_order)
setorder(coverage_df, Num_Causal, H2, Locus)

summary_method_labels <- c(
  "MESuSiE" = "MESuSiE",
  "SuSiE post-hoc" = "SuSiE",
  "SuSiE-weighted" = "SuSiE_weighted",
  "SuSiE-merged" = "SuSiE_merged",
  "PAINTOR" = "Paintor",
  "MultiSuSiE" = "MultiSuSiE",
  "SuSiEx" = "SuSiEx",
  "XMAP" = "XMAP",
  "CARMA-X" = "CARMAX"
)
summary_method_order <- unname(summary_method_labels)

coverage_table <- coverage_detail[
  ,
  .(coverage = round(mean(Covered), 3)),
  by = .(Method, H2_code, Causal_code)
]
coverage_table[
  ,
  Method := unname(summary_method_labels[Method])
]
coverage_table[
  ,
  `:=`(
    Method = factor(Method, levels = summary_method_order),
    h2 = fifelse(
      H2_code == 1L,
      "~h^2 == 10^-4",
      "~h^2 == 2%*%10^-4"
    ),
    causal_num = fcase(
      Causal_code == 1L, "Num~Causal == 1",
      Causal_code == 2L, "Num~Causal == 3",
      Causal_code == 3L, "Num~Causal == 5"
    )
  )
]
setorder(coverage_table, Method, H2_code, Causal_code)
coverage_table <- as_tibble(
  coverage_table[
    ,
    .(Method, h2, causal_num, coverage)
  ]
)

out_dir_env <- Sys.getenv("COVERAGE_REVISION_OUT_DIR", unset = "")
out_dir <- if (nzchar(out_dir_env)) {
  normalizePath(out_dir_env, mustWork = FALSE)
} else {
  normalizePath(
    file.path(simulation_dir, "..", "res_summary", "coverage_revision"),
    mustWork = FALSE
  )
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  coverage_df,
  file.path(out_dir, "coverage_simulation_baseline.csv")
)
save(
  coverage_df,
  file = file.path(out_dir, "coverage_simulation_baseline.RData")
)
fwrite(
  coverage_table,
  file.path(
    out_dir,
    "coverage_simulation_baseline_summary_table.csv"
  )
)
save(
  coverage_table,
  file = file.path(
    out_dir,
    "coverage_simulation_baseline_summary_table.RData"
  )
)

overall <- coverage_detail[
  ,
  .(
    N_CS = .N,
    N_covered = sum(Covered),
    Empirical_coverage = mean(Covered)
  ),
  by = Method
]
overall[, Method := factor(Method, levels = method_order_detail)]
setorder(overall, Method)
overall[, Method := as.character(Method)]

cat("Generated empirical per-CS coverage artifacts in:\n", out_dir, "\n")
cat("No fitted posterior mass, region-level union, or bootstrap was used.\n\n")
print(overall, row.names = FALSE)
