#!/bin/bash
#SBATCH -A gao824
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation_missing
#SBATCH --mem=8000
#SBATCH --array=1-800%400
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/slurm-%A_%a.err


bash

let k=0
for gene in {1..100}; do
    for h2 in {1..2}; do
        for causal in {1}; do
			 for causal_index in {1..2}; do
				for External_index in {1..2}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/4_shared_50_analysis.R ${gene} ${h2} ${causal} ${causal_index} ${External_index}
fi
done
done
done
done
done
