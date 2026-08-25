# MFD Analysis

Analysis code for benchmarking multi-ancestry fine-mapping methods, accompanying the manuscript on [**MFD**](https://github.com/gchen4422/MFD) (Multi-ancestry Fine-mapping Decision).

## Overview

This repository provides code for a systematic evaluation of nine multi-ancestry fine-mapping methods, classified into two broad categories:

- **Joint modeling** — methods that explicitly model multi-ancestry association signals within a unified probabilistic framework: [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0), [SuSiEx](https://github.com/getian107/SuSiEx), [MultiSuSiE](https://github.com/jordanero/MultiSuSiE), [XMAP](https://github.com/YangLabHKUST/XMAP), [MESuSiE](https://github.com/borangao/MESuSiE), and [CARMA-X](https://github.com/Iuliana-Ionita-Laza/CARMAX)
- **Univariate-based modeling** — methods that integrate single-ancestry results, further divided into association-level integration (SuSiE-weighted, SuSiE-merged) and post-hoc integration (SuSiE post-hoc)

Because no single method is uniformly optimal across diverse genetic architectures and LD conditions, we propose **MFD (Multi-ancestry Fine-mapping Decision)**, an adaptive strategy that selects among existing methods based on evidence for ancestry-specific variants informed by local LD structure.

Performance is evaluated across five dimensions: fine-mapping resolution, statistical accuracy (power and calibration), robustness to model misspecification, computational efficiency, and biological validity via functional enrichment.

Two complementary data resources are used:

1. **Simulation** — Imputed genotypes from UK Biobank, benchmarking across diverse genetic architectures and LD patterns across 100 1-Mb loci.
2. **Real_data** — 289 fine-mapping loci for BMI, diastolic blood pressure (DBP), and systolic blood pressure (SBP) from the All of Us Research Program.

## Methods Compared

**Joint modeling:** [MESuSiE](https://github.com/borangao/MESuSiE), [SuSiEx](https://github.com/getian107/SuSiEx), [MultiSuSiE](https://github.com/jordanero/MultiSuSiE), [XMAP](https://github.com/YangLabHKUST/XMAP), [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0), and [CARMA-X](https://github.com/Iuliana-Ionita-Laza/CARMAX)

**Univariate-based:** SuSiE-weighted, SuSiE-merged, SuSiE post-hoc

**Adaptive strategy:** MFD

---

## Repository Structure

```
MFD_analysis/
├── Simulation/
│   ├── simu_region_1mb/              # Main simulation (N=300K EUR + 300K AFR)
│   │   ├── shared_50/                # Baseline: 50% causal variants shared across ancestries
│   │   │   └── empirical_cs_coverage/ # Empirical coverage of individual credible sets
│   │   ├── shared_all/               # All causal variants shared (100%)
│   │   ├── shared_50_missing/        # Robustness: presence of ancestry-specific segregation variants
│   │   │   ├── Missing_Causal/
│   │   │   ├── Missing_Non_Causal/
│   │   │   └── endpoint_sensitivity/ # MFD endpoint sensitivity analyses
│   │   ├── shared_50_ld_correction/  # LD mismatch & correction scenarios
│   │   ├── functional_annotation_sensitivity/ # Functional-prior sensitivity
│   │   ├── polygenic_sensitivity/    # Higher-polygenicity sensitivity
│   │   └── utility.R                 # Shared functions (ROC, credible set metrics)
│   │
│   ├── simu_region_1mb_external_ld/  # Simulation with external LD reference panel
│   ├── simu_region_1mb_BB_50000/     # Smaller sample size scenario (N=50K for AFR)
│   ├── simu_region_admix/             # Admixed-population simulation
│   └── MultiSuSiE/                    # Separate Python code for MultiSuSiE benchmarking
│
└── Real_data/                         # Real GWAS application (BMI, DBP, SBP)
    └── three_ancestry/                # EUR, AFR, and EAS analysis code
```

Only analysis code is included for the added sensitivity and three-ancestry workflows. Input datasets, intermediate files, and generated result tables are not distributed in this repository.

---

## Simulation Workflow

Scripts are numbered sequentially within each scenario folder:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `1_sample_locus.R` | Sample 1-Mb loci from EUR/AFR reference panels |
| 2 | `2_locus_selection.R` | Filter and select loci for simulation |
| 3 | `3_shared_50_simulation.R` | Simulate phenotypes and compute summary statistics |
| 4 | `4_shared_50_analysis.R` | Run fine-mapping methods (MESuSiE, SuSiEx, XMAP, etc.) |
| 5 | `5_result_summarzing.R` | Aggregate results (power, FDR, credible set size) |
| 6 | `6_plot.R` | Generate figures |

SLALOM and METAL runs have dedicated sub-scripts (`metal/` subfolder; `_run_slalom_` scripts).  
SLURM job arrays are handled by corresponding `.sh` scripts.

The revision analyses add four simulation components:

- `functional_annotation_sensitivity/` evaluates the effect of functional priors.
- `polygenic_sensitivity/` evaluates loci with larger numbers of causal variants.
- `simu_region_admix/` evaluates admixed-population data and ancestry-matched or merged LD inputs.
- `shared_50/empirical_cs_coverage/` and `shared_50_missing/endpoint_sensitivity/` evaluate individual credible-set coverage and alternative MFD endpoints.

---

## Real Data Workflow

| Step | Script | Description |
|------|--------|-------------|
| 1 | `1_prepare_summstats.R` | QC and format AFR/EUR GWAS summary statistics |
| 2 | `2_define_region_1mb.R` | Define 1-Mb candidate regions around lead variants |
| 3 | `3_generate_loci.R` | Identify fine-mapping loci |
| 4 | `4_generate_loci_sumstats.R` | Extract per-locus summary statistics |
| 5 | `5_run_finemap_MF.R` | Run Multi-ancestry fine-mapping on each locus |
| 6 | `6_1_format_metal.R` / `6_2_run_metal_finemap.R` | Run METAL meta-analysis + finemap |
| 7 | `7_annotation.R` | Annotate credible set variants (eQTL, functional) |
| 8 | `8_result_summarizing.R` | Summarize results across traits and loci |
| 9–12 | `9_plot_stat.R`, `10–12_plot_*.R` | Generate enrichment and locus zoom plots |

SLURM submission scripts (`*.sh`) parallelize per-locus runs across computing nodes.

The `Real_data/three_ancestry/` directory contains the corresponding EUR/AFR/EAS workflow. Its main runner supplies three matched GWAS and LD inputs to `run_mf_decision`; the same interface extends to any number of ancestries.

---

## Dependencies

**R packages:** `MESuSiE`, `susieR`, `XMAP`, `data.table`, `snpStats`, `mvtnorm`, `ggplot2`, `patchwork`, `dplyr`, etc.

**External tools:** METAL, SuSiEx, SLALOM, MultiSuSiE (Python)

**Input data:** EUR and AFR LD reference panels (1000 Genomes); GWAS summary statistics from All of Us cohorts.
