#!/usr/bin/env bash
set -euo pipefail

# Let globs that don't match expand to nothing
shopt -s nullglob

########################################
# Paths
########################################

# CAUSAL_2 base dir
dir2_path="/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/causal_num_2/result"

# CAUSAL_3 base dir
dir3_path="/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/causal_num_3/result"

########################################
# Helper: count how many arguments (files) were passed
########################################
count_files() {
    echo "$#"   # if glob has no matches, $# = 0; when "-" is passed, the placeholder is ignored
}

cd "$dir2_path"

# Header: now includes SuSiEx / MultiSuSiE for both causal_2 and causal_3
printf "%-30s | %-10s | %-10s | %-10s | %-12s | %-10s | %-10s | %-10s | %-12s\n" \
"Directory" "Paintor(2)" "MESuSiE(2)" "SuSiEx(2)" "MultiSuSiE(2)" \
"MESuSiE(3)" "Paintor(3)" "SuSiEx(3)" "MultiSuSiE(3)"

echo "----------------------------------------------------------------------------------------------------------------------------------------"

for sub_dir in flipped_*; do
    [ -d "$sub_dir" ] || continue

    ########################################
    # CAUSAL_2 counts (current tree)
    ########################################
    c_paintor2=$(count_files "${sub_dir}/CAUSAL_2_LOCI_"*_h2_2*.mcmc.paintor)
    c_mesusie2=$(count_files "${sub_dir}/MESuSiE_CAUSAL_2_LOCI_"*_h2_2.RData)
    c_susiex2=$(count_files "${sub_dir}/susiex_result/"*snp)
    c_multisusie2=$(count_files "${sub_dir}/multisusiex_result/"*snp.txt)

    ########################################
    # CAUSAL_3 counts (mirrored subdir under dir3_path)
    ########################################
    if [ -d "${dir3_path}/${sub_dir}" ]; then
        c_mesusie3=$(count_files "${dir3_path}/${sub_dir}/MESuSiE_CAUSAL_3_LOCI_"*_h2_2.RData)
        c_paintor3=$(count_files "${dir3_path}/${sub_dir}/CAUSAL_3_LOCI_"*_h2_2*.mcmc.paintor)
        c_susiex3=$(count_files "${dir3_path}/${sub_dir}/susiex_result/"*snp)
        c_multisusie3=$(count_files "${dir3_path}/${sub_dir}/multisusiex_result/"*snp.txt)
    else
        c_mesusie3="-"
        c_paintor3="-"
        c_susiex3="-"
        c_multisusie3="-"
    fi

    printf "%-30s | %-10s | %-10s | %-10s | %-12s | %-10s | %-10s | %-10s | %-12s\n" \
        "$sub_dir" "$c_paintor2" "$c_mesusie2" "$c_susiex2" "$c_multisusie2" \
        "$c_mesusie3" "$c_paintor3" "$c_susiex3" "$c_multisusie3"
done

