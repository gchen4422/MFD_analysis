#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --job-name=simulation
#SBATCH --mem=32000
#SBATCH --partition=mulan
#SBATCH --array=1-200
#SBATCH --output=/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/out/Missing_Non_Causal_%a.out
#SBATCH --error=/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/out/Missing_Non_Causal_%a.error

bash

let k=0
for gene in {1..100}; do
	for causal_index in {1..2}; do
let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		Rscript --verbose /net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/code/Missing/Missing_Non_Causal/3_shared_50_simulation.R ${gene} ${causal_index}
fi
done
done
