# Functional-Annotation Sensitivity

This directory contains the code used to evaluate fine-mapping when causal variants are enriched for functional annotations. The simulation adds five binary annotations and an Evo 2 score, uses annotation-informed sampling probabilities, and compares PAINTOR runs with a constant annotation, the selected annotations, and the full annotation set.

The main sequence is:

1. `3_shared_50_simulation_functional.R` generates the summary statistics and annotation files.
2. `4_shared_50_analysis_functional_w0.R` and `4_shared_50_analysis_functional.R` run the `w0`, `w1`, and `w2` scenarios.
3. `5_result_summarzing_functional.R` and `5_result_summarzing_functional_w1.R` summarize accuracy, calibration, credible sets, and runtime.
4. `6_plot_functional.R` and `7_`–`9_` scripts create the sensitivity and enrichment summaries.
5. `metal/` contains the corresponding meta-analysis code.

Input annotation files, LD matrices, simulated summary statistics, and generated results are not included. Update the project, annotation, software, and output paths before running the scripts.
