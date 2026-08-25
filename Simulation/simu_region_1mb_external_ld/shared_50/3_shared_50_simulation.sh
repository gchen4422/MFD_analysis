#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation
#SBATCH --mem=40000







Rscript 3_shared_50_simulation.R
