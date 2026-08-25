#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_missing
#SBATCH --mem=32000
#SBATCH --array=1-600%300
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.err


bash

let k=0
for gene in {51..100}; do
    for h2 in 2; do
			 for flip_prop in {1..3}; do
				for rho_num in {1..4}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/4_shared_50_analysis_susie_rss_runxmap.R ${gene} ${h2} 2 ${flip_prop} ${rho_num}
fi
done
done
done
done
