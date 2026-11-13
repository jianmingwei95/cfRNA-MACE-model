# cfRNA-MACE-model

Analysis code for the manuscript *Circulating cell-free RNA profiles identify residual cardiovascular risk in patients with advanced atherosclerosis*.

A five-gene cell-free RNA (cfRNA) score (SPINT2, EIF4E, TRIM25, EI24, F11R) was developed for major adverse cardiovascular events (MACE) in preoperative plasma from patients undergoing carotid endarterectomy (discovery, n = 344), locked, and then applied without refitting to an independent serum cohort (replication, n = 277).

## Repository layout

| Path | Content |
|---|---|
| `cfRNA_MACE_nested_LOOCV.py` | Strict nested leave-one-out cross-validation (V1) |
| `Permutation test` | Permutation test for the LOOCV pipeline, HPC version (V1) |
| `Benchmark score2 model` | Benchmark against a clinical score with Cox regression (V1) |
| `model_development/` | Model development and internal validation in the discovery cohort (V2) |
| `external_validation/` | Application of the locked model to the independent replication cohort (V2) |

### `model_development/`

| Script | Purpose |
|---|---|
| `user_loocv_pipeline_fast.py` | Nested LOOCV feature selection and Cox model (Python; vectorised version of the original pipeline) |
| `loocv_topK_reimplementation.R` | Independent R re-implementation of the nested LOOCV with a configurable number of genes |
| `fig4_compute.R` | LOOCV selection frequency, univariable hazard ratios and cross-validated performance of the five-gene model across eight end points |
| `fig4_compute_locked.R` | Final five-gene model refitted on the full discovery cohort and locked (coefficients and risk threshold) |
| `internal_discovery_5gene_endpoints.R` | Apparent and 10-fold cross-validated performance in the discovery cohort |

### `external_validation/`

| Script | Purpose |
|---|---|
| `cfRNA_qc_batchcorrect_3batch277_g70.R` | Quality control, normalisation and batch correction of the replication cohort (three sequencing batches, n = 277) |
| `mace_subtype_validation_3batch.R` | Locked five-gene score applied without refitting; MACE and component end points |
| `fig5b_compute_replication_v2.R` | Discrimination and classification metrics with bootstrap confidence intervals in the replication cohort |
| `generalisation_internal_vs_external.R` | Internal (apparent, cross-validated) versus external performance across eight end points |

## Running

Paths are written as `path/to/project/...`; replace this prefix with your own project root, which is expected to contain `01_raw_data/`, `02_processed_data/` and `03_results/`.

Software: R (survival, glmnet, sva, limma, haven, data.table, survivalROC, pROC) and Python 3 (pandas, numpy, scikit-learn, lifelines, scikit-survival).

## Data

No patient-level data are included. cfRNA sequencing and clinical data from the Athero-Express biobank are available for academic research on reasonable request, subject to approval by the data access committee and a data transfer agreement: https://doi.org/10.34894/4IKE3T. Public data sets used: GSE16561 and GSE58294 (Gene Expression Omnibus).

## Licence

See `LICENSE`.
