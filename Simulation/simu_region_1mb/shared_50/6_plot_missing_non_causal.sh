#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --time=8:00:00
#SBATCH --job-name=MultiSuSiE

Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/6_plot_missing_non_causal.R







