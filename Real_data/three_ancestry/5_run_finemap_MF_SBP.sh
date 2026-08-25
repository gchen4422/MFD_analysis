#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=MFD3_SBP
#SBATCH --mem=32000
#SBATCH --array=1-61
#SBATCH --output=/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/logs/MFD3_SBP-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/logs/MFD3_SBP-%A_%a.err


# --- Define Variables ---

# 1. Get the Region ID from the Slurm Array ID (1 to 61)
LD_BLOCK=${SLURM_ARRAY_TASK_ID}

# 2. Set the Trait Index
# Based on R code: 1 = BMI, 2 = DBP, 3 = SBP
TRAIT_IDX=3

# --- Run the R Script ---
echo "Running SBP (Index 3) for Block ${LD_BLOCK}"

Rscript --verbose "${SLURM_SUBMIT_DIR}/5_run_finemap_MF.R" "${LD_BLOCK}" "${TRAIT_IDX}"
