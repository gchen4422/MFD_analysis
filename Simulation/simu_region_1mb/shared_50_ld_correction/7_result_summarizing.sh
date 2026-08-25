#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --job-name=Res
#SBATCH --mem=16000
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --array=1-36
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/LD_mismatch_summarize_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/LD_mismatch_summarize_%a.error

set -euo pipefail

k=0
for flip_prop in {1..3}; do
  for rho_num in {1..4}; do
    for ld_correct in {1..3}; do
      k=$((k+1))
      if [ "$k" -eq "$SLURM_ARRAY_TASK_ID" ]; then
        echo "Running flip_prop=${flip_prop}, rho_num=${rho_num}, ld_correct=${ld_correct}"
        Rscript --verbose \
          /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/7_result_summarizing.R \
          "${flip_prop}" "${rho_num}" "${ld_correct}"
        exit 0
      fi
    done
  done
done

echo "No parameter combination matched SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
exit 1
