#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=func_simu
#SBATCH --mem=40000
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/logs/slurm-3_simu-%j.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/logs/slurm-3_simu-%j.err


Rscript /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/3_shared_50_simulation_functional.R
