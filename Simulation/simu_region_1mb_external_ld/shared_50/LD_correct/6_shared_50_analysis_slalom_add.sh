#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_ld_correct
#SBATCH --mem=50000


Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/LD_correct/6_shared_50_analysis_slalom_add.R 38 2 3

