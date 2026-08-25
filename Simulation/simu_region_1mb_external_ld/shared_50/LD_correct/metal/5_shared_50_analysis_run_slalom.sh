#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=166:00:00
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --job-name=slalom_run
#SBATCH --mem=32000
#SBATCH --array=1-600
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/out/slurm-%A_%a.err



module purge
module load biocontainers
module load hail
module load spark


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
BASE_DIR=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50

mkdir -p ${BASE_DIR}/out

let k=0
for num_causal in {1..3}; do
  for h2_num in {1..2}; do
        for LD_BLOCK in {1..100}; do

          let k=${k}+1
          if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

            data_dir=${BASE_DIR}/causal_num_${num_causal}/summary_data
            output_dir=${data_dir}/metal_slalom
            mkdir -p ${output_dir}

            metal_snp=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_metal.snp
            metal_out=${output_dir}/CAUSAL_${num_causal}_LOCI_${LD_BLOCK}_h2_${h2_num}_metal.slalom.txt


            echo "[run] task=${SLURM_ARRAY_TASK_ID} causal=${num_causal} h2=${h2_num} loci=${LD_BLOCK}"

            # ---- Metal run ----
            if [ -f "${metal_snp}" ]; then
              echo "  [metal] ${metal_snp}"
              python3 ${SLALOM_PY} \
                --snp ${metal_snp} \
                --out ${metal_out} \
                --lead-variant-choice "prob" \
                --weighted-average-r afr=n_afr nfe=n_nfe \
                --align-alleles \
                --dentist-s \
                --abf \
                --reference-genome GRCh37 \
                --export-r
            else
              echo "  [skip metal] missing ${metal_snp}"
            fi


            exit 0
          fi

    done
  done
done

