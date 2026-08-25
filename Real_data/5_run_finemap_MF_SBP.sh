#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_BMI
#SBATCH --mem=32000
#SBATCH --array=1-46
#SBATCH --output=/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/out/slurm-%A_%a.err


# --- Define Variables ---

# 1. Get the Region ID from the Slurm Array ID (1 to 204)
LD_BLOCK=${SLURM_ARRAY_TASK_ID}

# 2. Set the Trait Index
# Based on the trait-index convention: 1 = BMI, 2 = DBP, 3 = SBP
TRAIT_IDX=3

# --- Run the R Script ---
echo "Running BMI (Index 1) for Block ${LD_BLOCK}"

Rscript --verbose /scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/5_run_finemap_MF.R "${LD_BLOCK}" "${TRAIT_IDX}"