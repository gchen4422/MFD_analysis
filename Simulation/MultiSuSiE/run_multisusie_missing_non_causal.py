#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import time
import MultiSuSiE

# Define ranges
num_causal_range  = range(1, 4)    # 1..3
LD_BLOCK_range    = range(1, 101)  # 1..100
h2_num_range      = range(1, 3)    # 1..2
causal_indexs     = range(1, 3)    # 1..2
External_indexs   = range(1, 3)    # 1..2

for num_causal in num_causal_range:
    for causal_index in causal_indexs:
        for External_index in External_indexs:
            for LD_BLOCK in LD_BLOCK_range:
                for h2_num in h2_num_range:

                    # Names matching R paste0 behavior
                    causal_index_name   = ("Both", "One")[causal_index - 1]
                    External_index_name = ("", "External_")[External_index - 1]

                    # Working dirs
                    wrk_dir = os.path.join(
                        "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb",
                        "shared_50_missing",
                        "Missing_Non_Causal",
                        f"{External_index_name}{causal_index_name}",
                        f"causal_num_{num_causal}",
                    ) + "/"

                    data_dir   = os.path.join(wrk_dir, "summary_data/")
                    result_dir = os.path.join(wrk_dir, "result/")
                    os.makedirs(os.path.join(result_dir, "multisusiex_result"), exist_ok=True)

                    # ----- Read zfile -----
                    zfile_path = f"{data_dir}CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{int(h2_num)}"
                    try:
                        zfile = np.genfromtxt(zfile_path, delimiter=' ', names=True, dtype=None, encoding='utf-8')
                    except Exception as e:
                        print(f"[zfile] CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num} read error: {e}")
                        continue
                    if getattr(zfile, "size", 0) == 0:
                        print(f"[zfile] Empty: CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}")
                        continue

                    # ----- Missingness filter -----
                    if causal_index == 1:
                        missing_mask = (zfile['EU_missing'] != 0) | (zfile['BB_missing'] != 0)
                    elif causal_index == 2:
                        missing_mask = (zfile['EU_missing'] != 0)
                    else:
                        print(f"[missingness] invalid causal_index={causal_index}")
                        continue

                    keep_mask = ~missing_mask
                    keep_idx  = np.nonzero(keep_mask)[0]
                    if keep_idx.size == 0:
                        print(f"[missingness] All SNPs filtered out: CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}")
                        continue

                    zfile_subset = zfile[keep_mask]

                    # Z-scores
                    zscore_1 = zfile_subset['zscore_1']
                    zscore_2 = zfile_subset['zscore_2']
                    N_list   = [300000, 300000]
                    z_list   = [zscore_1, zscore_2]

                    # ----- LD matrices -----
                    try:
                        EU_cov = np.genfromtxt(
                            f"{data_dir}CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{int(h2_num)}.LD1",
                            delimiter=' '
                        )
                        BB_cov = np.genfromtxt(
                            f"{data_dir}CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{int(h2_num)}.LD2",
                            delimiter=' '
                        )
                    except Exception as e:
                        print(f"[LD] read error CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}: {e}")
                        continue

                    m = zfile.shape[0]
                    if EU_cov.shape != (m, m) or BB_cov.shape != (m, m):
                        print(f"[LD] shape mismatch (expected {m}x{m}): EU={EU_cov.shape}, BB={BB_cov.shape}")
                        continue

                    EU_cov_subset = EU_cov[np.ix_(keep_idx, keep_idx)]
                    BB_cov_subset = BB_cov[np.ix_(keep_idx, keep_idx)]
                    R_list = [EU_cov_subset, BB_cov_subset]

                    # ----- MAF -----
                    try:
                        maf_EU = np.loadtxt(
                            f"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/maf_loci_{LD_BLOCK}_maf.txt"
                        )
                        maf_BB = np.loadtxt(
                            f"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/maf_loci_{LD_BLOCK}_maf.txt"
                        )
                    except Exception as e:
                        print(f"[MAF] read error LD_BLOCK_{LD_BLOCK}: {e}")
                        continue

                    if maf_EU.shape[0] != m or maf_BB.shape[0] != m:
                        print(f"[MAF] length mismatch (expected {m}): EU={maf_EU.shape[0]}, BB={maf_BB.shape[0]}")
                        continue
                    maf_EU_subset = maf_EU[keep_idx]
                    maf_BB_subset = maf_BB[keep_idx]
                    maf_list = [maf_EU_subset, maf_BB_subset]

                    # ----- Run MultiSuSiE -----
                    start_time = time.time()
                    try:
                        ss_fit = MultiSuSiE.multisusie_rss(
                            z_list=z_list,
                            R_list=R_list,
                            rho=np.array([[1, 0.8], [0.8, 1]]),
                            population_sizes=N_list,
                            L=10,
                            scaled_prior_variance=0.2,
                            low_memory_mode=False,
                            min_abs_corr=0.5,
                            single_population_mac_thresh=20,
                            maf_list=maf_list,
                            coverage=0.95
                        )
                    except Exception as e:
                        print(f"[MultiSuSiE] error CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}: {e}")
                        continue

                    # ----- Save results -----
                    time_taken = (time.time() - start_time) / 60.0
                    file_path = os.path.join(
                        result_dir,
                        "multisusiex_result",
                        f"MultiSuSiE_CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}_output"
                    )

                    with open(f"{file_path}_runtime.txt", "w") as f:
                        f.write(f"Time taken: {time_taken:.8f} minutes\n")

                    df_zfile = pd.DataFrame(zfile_subset)
                    df_zfile["pip"] = ss_fit.pip
                    df_zfile.to_csv(f"{file_path}_snp.txt", sep="\t", index=False)

                    try:
                        filtered_sets = [
                            ss_fit.sets[0][i]
                            for i in range(len(ss_fit.sets[0]))
                            if ss_fit.sets[3][i]
                        ]
                    except Exception:
                        filtered_sets = []

                    filtered_df = pd.DataFrame()
                    for idx, cs_set in enumerate(filtered_sets):
                        for snp_index in cs_set:
                            if snp_index < len(df_zfile):
                                row = df_zfile.iloc[[snp_index]].copy()
                                row["CS"] = idx + 1
                                filtered_df = pd.concat([filtered_df, row], ignore_index=True)

                    filtered_df.to_csv(f"{file_path}_cs.txt", sep="\t", index=False)
                    print(
                        f"[DONE] num_causal={num_causal} | causal_index={causal_index} ({causal_index_name}) | "
                        f"External_index={External_index} ({External_index_name or 'None'}) | "
                        f"LD_BLOCK={LD_BLOCK} | h2_num={h2_num} | minutes={time_taken:.2f}",
                        flush=True
                    )
