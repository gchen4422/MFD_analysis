#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --job-name=Res
#SBATCH --mem=16000
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --array=1-4
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/Missing_non_causal_summarize_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/out/Missing_non_causal_summarize_%a.error

bash

let k=0
	for causal_index in {1..2}; do
		for External_index in {1..2}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/5_result_summarizing.R ${causal_index} ${External_index}
fi
done
done

