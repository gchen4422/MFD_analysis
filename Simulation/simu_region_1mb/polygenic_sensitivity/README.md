# Polygenic Sensitivity

This directory contains the higher-polygenicity simulation used to test MFD and the comparison methods when a locus contains more causal variants. The three simulation codes generate 15, 20, or 25 total causal variants, split between cross-ancestry and ancestry-specific effects. The main fine-mapping runner allows up to 20 effects for MESuSiE, SuSiE, SuSiEx, and XMAP.

The main sequence is:

1. `3_shared_50_simulation_polygenic.R` generates the simulated summary statistics and LD inputs.
2. `4_shared_50_analysis.R` runs the fine-mapping methods.
3. `5_result_summarzing_meta.R` and `6_plot.R` summarize performance and credible-set results.
4. `metal/` contains the corresponding meta-analysis code.

The summarization script evaluates causal-architecture codes 1 and 2; code 3 is retained for the 25-causal-variant extension. Input genotypes, generated summary statistics, and result files are not included. Update local paths before running the scripts.
