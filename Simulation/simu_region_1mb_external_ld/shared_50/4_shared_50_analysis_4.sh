#!/bin/bash

  
for gene in {61..80}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/4_shared_50_analysis.R ${gene} ${h2} ${causal};  done;  done;  done
