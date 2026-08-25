#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import time
import MultiSuSiE

# Define ranges
num_causal_range  = range(2, 4)    # 2,3
LD_BLOCK_range    = range(1, 101)  # 1..100 (python range end is exclusive)
h2_num_range      = [2]       # 2
flip_props        = range(1, 4)    # 1,2,3
rho_nums          = range(1, 5)    # 1,2,3,4

for num_causal in num_causal_range:
    for rho_num in rho_nums:
        for flip_prop in flip_props:
            for LD_BLOCK in LD_BLOCK_range:
                for h2_num in h2_num_range:

                    # Working dirs
                    wrk_dir = os.path.join(
                        "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb",
                        "shared_50_ld_correction",
                        f"causal_num_{num_causal}",
                    ) + "/"

                    data_dir   = os.path.join(wrk_dir, "summary_data/")
                    data_dir_flipped = os.path.join(data_dir, f"flipped_{flip_prop}_rho_{rho_num}") + "/"
                    result_dir = os.path.join(wrk_dir, f"result/flipped_{flip_prop}_rho_{rho_num}") + "/"
                    os.makedirs(os.path.join(result_dir, "multisusiex_result"), exist_ok=True)

                    # ----- Read zfile (structured array) -----
                    zfile_path = (
                        f"{data_dir_flipped}"
                        f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{int(h2_num)}"
                        f"_flipped_{flip_prop}_rho_{rho_num}"
                    )

                    try:
                        zfile = np.genfromtxt(
                            zfile_path, delimiter=' ', names=True,
                            dtype=None, encoding='utf-8'
                        )
                    except Exception as e:
                        print(f"[zfile] CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num} read error: {e}")
                        continue
                    if zfile.size == 0:
                        print(f"[zfile] Empty: CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}")
                        continue

                    zscore_1 = zfile['zscore_1']
                    zscore_2 = zfile['zscore_2']
                    N_list = [300000, 300000]
                    z_list = [zscore_1, zscore_2]

                    # Reading covariance matrices
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
                        print(f"Error reading covariance matrices for CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}: {e}")
                        continue

                    R_list = [EU_cov, BB_cov]

                    # ----- MAFs -----
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

                    maf_list = [maf_EU, maf_BB]

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

                    # ----- Save outputs -----
                    time_taken = (time.time() - start_time) / 60.0

                    file_path = os.path.join(
                        result_dir,
                        "multisusiex_result",
                        f"MultiSuSiE_CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}_output"
                    )

                    # runtime
                    with open(f"{file_path}_runtime.txt", "w") as f:
                        f.write(f"Time taken: {time_taken:.8f} minutes\n")

                    # SNP-level table with PIP
                    df_zfile = pd.DataFrame(zfile)
                    df_zfile["pip"] = ss_fit.pip
                    df_zfile.to_csv(f"{file_path}_snp.txt", sep="\t", index=False)

                    # Credible sets (keep only significant)
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
                        f"[DONE] num_causal={num_causal} | flipped={flip_prop} | rho={rho_num} | "
                        f"LD_BLOCK={LD_BLOCK} | h2_num={h2_num} | minutes={time_taken:.2f}",
                        flush=True
                    )

