#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation
#SBATCH --mem=50000
#SBATCH --array=38,56
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/out/slurm-%A_%a.err

# Get the SLURM task ID (1-based index)
GENE=${SLURM_ARRAY_TASK_ID}

# Define ranges for h2 and causal
H2_VALUES=(1 2)
CAUSAL_VALUES=(1 2 3)

# Loop over h2 and causal combinations for the current gene
for H2 in "${H2_VALUES[@]}"; do
  for CAUSAL in "${CAUSAL_VALUES[@]}"; do
    # Run the R script for each combination
    Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_BB_50000/shared_50/4_shared_50_analysis.R "${GENE}" "${H2}" "${CAUSAL}"
  done
done







