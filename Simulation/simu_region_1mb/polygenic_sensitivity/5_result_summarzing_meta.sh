#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=4:00:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --job-name=res_poli
#SBATCH --mem=32000
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/out/slurm-res_summarize-%j.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/out/slurm-res_summarize-%j.err

set -e
source ~/.bashrc

Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polygenic_sensitivity/5_result_summarzing_meta.R

echo "Done: result summarizing"
