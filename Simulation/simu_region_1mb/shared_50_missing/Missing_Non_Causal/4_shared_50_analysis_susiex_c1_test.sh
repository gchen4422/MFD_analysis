#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_missing
#SBATCH --mem=8000

bash


Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/4_shared_50_analysis_susiex.R 76 2 1 1 1
Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/4_shared_50_analysis_susiex.R 77 1 1 2 1
