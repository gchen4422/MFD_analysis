#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=4:00:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --job-name=simulation
#SBATCH --mem=32000



Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/5_result_summarzing_meta.R
