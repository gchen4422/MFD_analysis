#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --job-name=Res_MFD_sens
#SBATCH --mem=16000
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --array=1-2
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/out/Missing_causal_summarize_MFD_sensitivity_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/out/Missing_causal_summarize_MFD_sensitivity_%a.error

let k=0
	for causal_index in 2; do
		for External_index in {1..2}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript --verbose "${SLURM_SUBMIT_DIR}/summarize_missing_causal_endpoints.R" "${causal_index}" "${External_index}"
fi
done
done
