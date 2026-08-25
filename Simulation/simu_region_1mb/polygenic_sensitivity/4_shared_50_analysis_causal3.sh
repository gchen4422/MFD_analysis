#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=144:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=poli_c3
#SBATCH --mem=32000
#SBATCH --array=1-100
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/out/slurm-%A_%a.err


# Get the SLURM task ID (1-based index)
GENE=${SLURM_ARRAY_TASK_ID}

# Only causal_3
H2_VALUES=(1 2)
CAUSAL_VALUES=(3)

# Loop over h2 and causal combinations for the current gene
for H2 in "${H2_VALUES[@]}"; do
  for CAUSAL in "${CAUSAL_VALUES[@]}"; do
    # Skip if result already exists
    RDATA="/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/causal_num_${CAUSAL}/result/MESuSiE_CAUSAL_${CAUSAL}_LOCI_${GENE}_h2_${H2}.RData"
    if [ ! -f "$RDATA" ]; then
      Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/4_shared_50_analysis.R "${GENE}" "${H2}" "${CAUSAL}"
    fi
  done
done
