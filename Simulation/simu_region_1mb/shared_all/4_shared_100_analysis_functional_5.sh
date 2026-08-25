#!/bin/bash














for gene in {81..100}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/4_shared_100_analysis_functional.R ${gene} ${h2} ${causal};  done;  done;  done
