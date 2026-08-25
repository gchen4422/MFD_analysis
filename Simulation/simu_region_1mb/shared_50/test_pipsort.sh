#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=simulation
#SBATCH --mem=8000



time /home/chen4422/gwas_software/pipsort/PIPSORT -c 3 -l '/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2/summary_data/CAUSAL_2_LOCI_1_h2_2_ldfiles.txt' -z '/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2/summary_data/CAUSAL_2_LOCI_1_h2_2_zfiles.txt' -m '/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2/summary_data/CAUSAL_2_LOCI_1_h2_2_snp_map' -n 300000,300000 -p 0.25 -o '/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2/result/CAUSAL_2_LOCI_1_h2_2_pipsort'
