#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=144:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation
#SBATCH --mem=32000



Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/LD_correct/7_result_summarizing_meta_slalom.R
