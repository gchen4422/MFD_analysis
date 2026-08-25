# MFD Endpoint Sensitivity

This directory contains post-processing code for the MFD endpoint sensitivity analysis. It applies MESuSiE, SuSiEx, and XMAP as alternative joint-model endpoints to the missing-causal and missing-noncausal simulation results, then compares their accuracy and LD-distance behavior.

- `summarize_missing_causal_endpoints.R` summarizes missing-causal scenarios.
- `summarize_missing_noncausal_endpoints.R` summarizes missing-noncausal scenarios.
- `plot_endpoint_sensitivity.R` combines the endpoint-specific summaries into figures and tables.

The scripts use the existing outputs under `Missing_Causal/` and `Missing_Non_Causal/`; they do not rerun the fine-mapping methods. Simulation outputs and generated tables are not included.
