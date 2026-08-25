#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=get_loci
#SBATCH --mem=60000



Rscript 4_generate_loci_sumstats.R
