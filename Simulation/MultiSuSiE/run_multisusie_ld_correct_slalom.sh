#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --time=108:00:00
#SBATCH --job-name=MultiSuSiE

module load anaconda
conda activate MultiSuSiE

python -u run_multisusie_ld_correct_slalom.py
~                                                
