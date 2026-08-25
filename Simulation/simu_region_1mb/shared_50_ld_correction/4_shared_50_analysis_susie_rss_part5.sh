#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_missing
#SBATCH --mem=60000
#SBATCH --array=0-20
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.err

SCRIPT=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/4_shared_50_analysis_susie_rss.R

# List of 5-arg parameter sets (space-separated)
PARAMS=(
  "68 2 2 1 1"
  "68 2 2 1 3"
  "68 2 2 1 4"
  "68 2 2 2 3"

  "38 2 3 1 1"
  "38 2 3 1 2"
  "38 2 3 1 3"
  "38 2 3 1 4"
  "38 2 3 2 1"
  "38 2 3 2 2"
  "38 2 3 2 3"
  "38 2 3 2 4"
  "38 2 3 3 1"
  "38 2 3 3 2"
  "38 2 3 3 3"
  "38 2 3 3 4"

  "56 2 3 1 1"
  "56 2 3 1 2"
  "56 2 3 1 4"
  "56 2 3 2 4"
  "56 2 3 3 4"
)

# Pick the parameter set for this task
args="${PARAMS[$SLURM_ARRAY_TASK_ID]}"

echo "Task ${SLURM_ARRAY_TASK_ID} running: Rscript --verbose $SCRIPT $args"
Rscript --verbose "$SCRIPT" $args

