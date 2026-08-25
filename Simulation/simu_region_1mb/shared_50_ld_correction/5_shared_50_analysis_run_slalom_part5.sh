#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=slalom_run_sel
#SBATCH --mem=32000
#SBATCH --array=1-5%5
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/out/slurm-%A_%a.err

module purge
module load biocontainers
module load hail
module load spark

# (optional but helps apptainer/fs flakiness)
export APPTAINER_CACHEDIR=/scratch/negishi/chen4422/.apptainer_cache
export APPTAINER_TMPDIR=/scratch/negishi/chen4422/.apptainer_tmp
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

export PYSPARK_SUBMIT_ARGS="\
 --driver-memory 8g \
 --executor-memory 8g \
 --jars /home/chen4422/gwas_software/gcs-connector-hadoop3-latest.jar \
 --conf spark.hadoop.fs.gs.impl=com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem \
 --conf spark.hadoop.fs.AbstractFileSystem.gs.impl=com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS \
 --conf spark.hadoop.google.cloud.auth.service.account.enable=true \
 --conf spark.hadoop.google.cloud.auth.service.account.json.keyfile=/home/chen4422/slalom-gwas-a5254bb748e5.json \
 --conf spark.hadoop.fs.gs.project.id=slalom-gwas \
 --conf spark.hadoop.fs.gs.requester.pays.mode=AUTO \
 --conf spark.hadoop.fs.gs.requester.pays.project.id=slalom-gwas \
 pyspark-shell"

SLALOM_PY=~/gwas_software/slalom/slalom.py
BASE_DIR=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction

# fixed parameters
num_causal=3
h2_num=2
rho_num=4

# the selected (flip, loci) pairs, in order of task_id = 1..5
flip_list=(1 1 2 2 2)
loci_list=(6 7 22 23 24)

idx=$((SLURM_ARRAY_TASK_ID - 1))
flip_prop=${flip_list[$idx]}
LD_BLOCK=${loci_list[$idx]}

data_dir_flipped=${BASE_DIR}/causal_num_${num_causal}/summary_data/flipped_${flip_prop}_rho_${rho_num}
output_dir=${data_dir_flipped}/slalom
mkdir -p ${output_dir}

eur_snp=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_eur.snp
eur_out=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_eur.slalom.txt

afr_snp=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_afr.snp
afr_out=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_afr.slalom.txt

echo "[run] task=${SLURM_ARRAY_TASK_ID} causal=${num_causal} h2=${h2_num} flip=${flip_prop} rho=${rho_num} loci=${LD_BLOCK}"

# ---- EUR run (map EUR to gnomAD NFE) ----
if [ -f "${eur_snp}" ]; then
  echo "  [eur] ${eur_snp}"
  python3 ${SLALOM_PY} \
    --snp ${eur_snp} \
    --out ${eur_out} \
    --lead-variant-choice "prob" \
    --weighted-average-r nfe=n_nfe \
    --align-alleles \
    --dentist-s \
    --abf \
    --reference-genome GRCh37 \
    --export-r
else
  echo "  [skip eur] missing ${eur_snp}"
fi

# ---- AFR run ----
if [ -f "${afr_snp}" ]; then
  echo "  [afr] ${afr_snp}"
  python3 ${SLALOM_PY} \
    --snp ${afr_snp} \
    --out ${afr_out} \
    --lead-variant-choice "prob" \
    --weighted-average-r afr=n_afr \
    --align-alleles \
    --dentist-s \
    --abf \
    --reference-genome GRCh37 \
    --export-r
else
  echo "  [skip afr] missing ${afr_snp}"
fi

