# Admixed-Population Simulation

This directory contains the code used to evaluate fine-mapping with an admixed population. The workflow simulates association statistics, constructs the merged LD matrix used in the analysis, runs MFD and the comparison methods, and summarizes performance.

The main sequence is:

1. `3_admix_simulation.R` generates the simulated association data.
2. `5_generate_merged_ld.R` constructs the merged LD input.
3. `4_admix_analysis.R` runs the fine-mapping methods. CARMA-X uses `LD.estimation = TRUE` in this admixed-population scenario.
4. `5_result_summarzing.R` and `6_plot.R` summarize and visualize performance.
5. `metal/` contains the corresponding meta-analysis code.

Reference genotypes, generated summary statistics, LD matrices, and result files are not included. Update local data, software, and output paths before running the workflow.
