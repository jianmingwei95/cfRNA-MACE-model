# cfRNA-MACE-model

Analysis code for the manuscript *Circulating cell-free RNA profiles identify residual cardiovascular risk in patients with advanced atherosclerosis*.

A five-gene cell-free RNA (cfRNA) score (SPINT2, EIF4E, TRIM25, EI24, F11R) was developed for major adverse cardiovascular events (MACE) in preoperative plasma from patients undergoing carotid endarterectomy (discovery, n = 344), locked, and then applied without refitting to an independent serum cohort (replication, n = 277).

## Repository layout

| Path | Content |
|---|---|
| `cfRNA_MACE_nested_LOOCV.py` | Gene selection: nested leave-one-out cross-validation in the discovery cohort (ten genes per fold) |
| `Permutation test` | Permutation test for the nested LOOCV pipeline (HPC version) |
| `model_development/` | Final five-gene model: fitting, locking and internal validation in the discovery cohort |
| `external_validation/` | Application of the locked model to the independent replication cohort |
| `clinical_benchmark/` | Comparison with established clinical risk scores and added value of the cfRNA score |

### `model_development/`

| Script | Purpose |
|---|---|
| `fit_locked_model.R` | Ridge-penalised Cox model with the five genes fitted on the full discovery cohort (penalty by 10-fold cross-validation); saves the locked model, the discovery reference distribution and the score tertiles |
| `fig4_compute.R` | Selection frequency and univariable hazard ratios of the ten most frequently selected genes; leave-one-out cross-validation of the five-gene model across eight end points |
| `fig4_compute_locked.R` | Apparent performance of the locked model in the discovery cohort and the locked risk threshold (upper tertile) |

### `external_validation/`

| Script | Purpose |
|---|---|
| `mace_subtype_validation_3batch.R` | Replication cohort (three sequencing batches, n = 277): quality control, normalisation, batch correction, reference quantile mapping and scoring with the locked model without refitting; C-index, 3-year AUC and hazard ratio per SD for eight end points |
| `fig5b_compute_replication_v2.R` | Discrimination and classification metrics at the locked threshold with bootstrap confidence intervals |

### `clinical_benchmark/`

| Script | Purpose |
|---|---|
| `benchmark_clinical_scores.R` | Construction of the SMART, SMART-REACH, SMART2 and CEA clinical risk scores from published coefficients, and comparison with and without the cfRNA score |
| `benchmark_locked_weights.R` | Clinical score plus cfRNA with combination weights estimated in discovery and applied unchanged to replication |
| `fig6_compute.R` | Risk stratification, calibration and decision-curve analysis of the locked score in the replication cohort |
| `fig6_lrt_nested.R` | Likelihood-ratio tests of nested Cox models (clinical score versus clinical score plus cfRNA) |

## Running

Paths are written as `path/to/project/...`; replace this prefix with your own project root, which is expected to contain `01_raw_data/`, `02_processed_data/` and `03_results/`. Two scripts in `clinical_benchmark/` read `benchmark_clinical_scores.R` or `benchmark_locked_weights.R` from `clinical_benchmark/` under the project root, so keep that folder name.

Software: R (survival, glmnet, sva, limma, haven, mice, data.table, survivalROC, timeROC, pROC, PRROC, dcurves) and Python 3 (pandas, numpy, scikit-learn, lifelines, scikit-survival).

## Data

No patient-level data are included. cfRNA sequencing and clinical data from the Athero-Express biobank are available for academic research on reasonable request, subject to approval by the data access committee and a data transfer agreement: https://doi.org/10.34894/4IKE3T. Public data sets used: GSE16561 and GSE58294 (Gene Expression Omnibus).

## Licence

See `LICENSE`.
