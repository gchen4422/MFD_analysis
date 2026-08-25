#!/usr/bin/env python3
"""
Run MultiSuSiE RSS fine-mapping across simulation grid.

This script loops over:
  - num_causal in {2,3}
  - rho_num in {1,2,3,4}
  - flip_prop in {1,2,3}
  - LD_BLOCK in {1..100}
  - h2_num = 2

For each setting it:
  1) loads flipped z-scores (slalom-filtered + original),
  2) loads LD matrices,
  3) aligns and subselects MAFs to the retained SNP set,
  4) runs MultiSuSiE.multisusie_rss,
  5) writes runtime, SNP PIPs, and credible sets.
"""

import os
import time
import numpy as np
import pandas as pd
import MultiSuSiE


def safe_genfromtxt(path, **kwargs):
    """Read a structured text file safely."""
    try:
        arr = np.genfromtxt(path, **kwargs)
        return arr
    except Exception as e:
        print(f"[read error] {path}: {e}", flush=True)
        return None


def safe_loadtxt(path, **kwargs):
    """Read a numeric text file safely."""
    try:
        arr = np.loadtxt(path, **kwargs)
        return arr
    except Exception as e:
        print(f"[read error] {path}: {e}", flush=True)
        return None


def main():
    # Define ranges
    num_causal_range = range(2, 4)    # 2..3
    LD_BLOCK_range   = range(1, 101)  # 1..100
    h2_num_range     = [2]
    flip_props       = range(1, 4)    # 1..3
    rho_nums         = range(1, 5)    # 1..4

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

                        data_dir = os.path.join(wrk_dir, "summary_data/")
                        data_dir_flipped = os.path.join(
                            data_dir, f"flipped_{flip_prop}_rho_{rho_num}"
                        ) + "/"
                        result_dir = os.path.join(
                            wrk_dir, f"result/flipped_{flip_prop}_rho_{rho_num}_slalom"
                        ) + "/"

                        os.makedirs(os.path.join(result_dir, "multisusiex_result"), exist_ok=True)

                        # ----- Read z files -----
                        zfile_path = (
                            f"{data_dir_flipped}"
                            f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}"
                            f"_flipped_{flip_prop}_rho_{rho_num}_slalom"
                        )
                        zfile_path_ori = (
                            f"{data_dir_flipped}"
                            f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}"
                            f"_flipped_{flip_prop}_rho_{rho_num}"
                        )

                        zfile = safe_genfromtxt(
                            zfile_path,
                            delimiter=" ",
                            names=True,
                            dtype=None,
                            encoding="utf-8",
                        )
                        zfile_ori = safe_genfromtxt(
                            zfile_path_ori,
                            delimiter=" ",
                            names=True,
                            dtype=None,
                            encoding="utf-8",
                        )

                        if zfile is None or zfile_ori is None:
                            print(
                                f"[zfile] CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num} read error",
                                flush=True
                            )
                            continue
                        if zfile.size == 0:
                            print(
                                f"[zfile] Empty: CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}",
                                flush=True
                            )
                            continue

                        zscore_1 = zfile["zscore_1"]
                        zscore_2 = zfile["zscore_2"]
                        N_list = [300000, 300000]
                        z_list = [zscore_1, zscore_2]

                        # ----- Read covariance matrices -----
                        EU_cov = safe_genfromtxt(
                            f"{data_dir_flipped}"
                            f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}"
                            f"_flipped_{flip_prop}_rho_{rho_num}_slalom.LD1",
                            delimiter=" ",
                        )
                        BB_cov = safe_genfromtxt(
                            f"{data_dir_flipped}"
                            f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}"
                            f"_flipped_{flip_prop}_rho_{rho_num}_slalom.LD2",
                            delimiter=" ",
                        )

                        if EU_cov is None or BB_cov is None:
                            print(
                                f"[LD] Error reading covariance matrices for "
                                f"CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}",
                                flush=True
                            )
                            continue

                        R_list = [EU_cov, BB_cov]

                        # ----- MAF alignment -----
                        ori_rsid = np.asarray(zfile_ori["RSID"])
                        new_rsid = np.asarray(zfile["RSID"])

                        keep_mask = np.isin(ori_rsid, new_rsid)

                        ori_subset_rsid = ori_rsid[keep_mask]
                        new_rsid_arr = np.asarray(zfile["RSID"])

                        print("same length:", ori_subset_rsid.size == new_rsid_arr.size, flush=True)
                        print("same order:", np.array_equal(ori_subset_rsid, new_rsid_arr), flush=True)

                        maf_EU = safe_loadtxt(
                            f"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/"
                            f"risk_loci_ld_eur/maf_loci_{LD_BLOCK}_maf.txt"
                        )
                        maf_BB = safe_loadtxt(
                            f"/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/"
                            f"risk_loci_ld_afr/maf_loci_{LD_BLOCK}_maf.txt"
                        )

                        if maf_EU is None or maf_BB is None:
                            print(f"[MAF] read error LD_BLOCK_{LD_BLOCK}", flush=True)
                            continue

                        maf_EU_sub = maf_EU[keep_mask]
                        maf_BB_sub = maf_BB[keep_mask]
                        maf_list = [maf_EU_sub, maf_BB_sub]

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
                                coverage=0.95,
                            )
                        except Exception as e:
                            print(
                                f"[MultiSuSiE] error CAUSAL_{num_causal}_LOCI_{LD_BLOCK}_h2_{h2_num}: {e}",
                                flush=True
                            )
                            continue

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


if __name__ == "__main__":
    main()

