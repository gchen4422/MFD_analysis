#!/bin/bash











#for gene in {1..70}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_ld/shared_50/4_shared_50_analysis.R ${gene} ${h2} ${causal};  done;  done;  done



#for gene in {72..75}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_ld/shared_50/4_shared_50_analysis.R ${gene} ${h2} ${causal};  done;  done;  done


#for gene in {77..100}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_ld/shared_50/4_shared_50_analysis.R ${gene} ${h2} ${causal};  done;  done;  done



for gene in {91..100}; do     for h2 in {1..2}; do         for causal in {1..3}; do Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/4_shared_100_analysis.R ${gene} ${h2} ${causal};  done;  done;  done
