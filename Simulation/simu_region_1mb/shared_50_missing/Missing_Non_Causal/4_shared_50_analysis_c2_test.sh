#!/bin/bash
#SBATCH -A gao824
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_missing_g68
#SBATCH --mem=16000
#SBATCH --array=1-8%8
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/slurm-%A_%a.err

set -euo pipefail

gene=68
k=0
for h2 in 1 2; do
  for causal_index in 1 2; do
    for External_index in 1 2; do
      k=$((k+1))
      if [ "$k" -eq "${SLURM_ARRAY_TASK_ID}" ]; then
        Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/4_shared_50_analysis.R \
          "$gene" "$h2" 2 "$causal_index" "$External_index"
      fi
    done
  done
done

