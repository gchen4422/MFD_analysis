#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_ld_correct
#SBATCH --mem=16000
#SBATCH --array=1-600
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/out/slurm-%A_%a.err


bash

let k=0
for gene in {1..100}; do
    for h2 in {1..2}; do
			 for num_causal in {1..3}; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/LD_correct/6_shared_50_analysis_slalom.R ${gene} ${h2} ${num_causal}
	
fi

done
done
done
