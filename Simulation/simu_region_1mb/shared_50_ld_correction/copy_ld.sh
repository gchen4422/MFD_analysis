#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --partition=cpu
#SBATCH --qos=normal
#SBATCH --time=80:00:00
#SBATCH --job-name=MultiSuSiE




root="/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction"

for causal in 1 2 3; do
  base="$root/causal_num_${causal}/summary_data"
  [ -d "$base" ] || continue

  echo "Processing $base"

  # loop over flipped dirs inside this causal_num
  for d in "$base"/flipped_*_rho_*; do
    [ -d "$d" ] || continue

    bn=$(basename "$d")                 # e.g. flipped_1_rho_4
    flip_num=${bn#flipped_}             # 1_rho_4
    flip_num=${flip_num%%_rho_*}        # 1
    rho_num=${bn##*_rho_}               # 4

    # copy all LD1/LD2 files from base into this flipped dir with new names
    for f in "$base"/CAUSAL_${causal}_LOCI_*_h2_*.LD1 \
             "$base"/CAUSAL_${causal}_LOCI_*_h2_*.LD2; do
      [ -e "$f" ] || continue

      fname=$(basename "$f")            # CAUSAL_1_LOCI_1_h2_1.LD1
      stem=${fname%.LD1}
      stem=${stem%.LD2}                # CAUSAL_1_LOCI_1_h2_1
      ext=${fname##*.}                 # LD1 or LD2

      cp "$f" "$d/${stem}_flipped_${flip_num}_rho_${rho_num}.${ext}"
    done
  done
done
