# Empirical Credible-Set Coverage

These scripts calculate truth-based empirical coverage for each credible set reported in the baseline simulation. Coverage is the proportion of individual reported credible sets that contain at least one simulated causal variant; credible sets are not first collapsed into a region-level union.

- `validate_empirical_cs_coverage.R` extracts each method's credible sets, matches them to simulated truth, and performs consistency checks.
- `summarize_empirical_cs_coverage.R` creates the coverage summaries used for reporting.

The scripts consume outputs from the baseline `shared_50` workflow. No generated coverage tables or simulation outputs are included.
