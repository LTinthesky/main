# Bulk RNA-seq Gene Module Scoring and CTA Associations

R code for scoring targeted myeloid gene modules in transformed bulk RNA-seq expression and assessing their associations with CTA remodeling index and an RNA remodeling score.

## Analyses

1. `scripts/01_targeted_myeloid_module_scores.R` maps manually curated neutrophil and plaque-myeloid marker lists to Ensembl IDs, combines them with an existing CTA-derived feature set, calculates GSVA and mean gene-z scores, and summarizes continuous and binary RI associations.
2. `scripts/02_compare_independent_modules.R` compares the curated-module and Hallmark-derived scores with previously calculated CTA-derived scores. It requires the original prior outputs; it is not a standalone reproduction of the elastic-net analysis.

The curated marker lists are defined explicitly in the first script. Their labels describe intended biology; they are not certified reproductions of a particular published signature. The comparison uses the same patients. “Independent” in the original filenames refers to gene-set definitions, not an independent validation cohort. Associations involving the CTA-derived selected genes are descriptive same-cohort analyses.

## Requirements and inputs

R >= 4.2; SummarizedExperiment, AnnotationDbi, org.Hs.eg.db, GSVA with the parameter-object API, BiocParallel, ggplot2, readr, dplyr, tidyr, and tibble.

Set `CTA_PROJECT_DIR` to a local directory with the expected frozen-analysis layout. Default: `data/cta_project`. From the repository root, run:

```sh
Rscript scripts/01_targeted_myeloid_module_scores.R
Rscript scripts/02_compare_independent_modules.R
```

Script 01 expects these relative inputs:

- `results_frozen_analysis_v1/step1_preprocess_1125/vst_1125_allGenes_from_combatseq.rds`: a SummarizedExperiment whose default assay contains finite, transformed expression (genes by samples).
- `results_frozen_analysis_v1/step3_feature_sets/feature_set_fromMergedTable_p001_sameDir_allDirections_ensembl.csv`: `direction`, `symbol`, and `gene` columns. Positive rows have `direction = positive_remodeling`.
- `results_frozen_analysis_v1/step6_apply_1125_results/run_v1/cta_assoc_subset51_from1125score/subset51_scoresFrom1125_with_CTA.csv`: unique `sample`, numeric `RI`, numeric `RI_score`, and `RI_binary` (`Nonpositive_remodeling`/`Positive_remodeling`, missing for the mid-zone).

Expression IDs must match unversioned Ensembl IDs used by the annotation package. Gene overlap and unmapped markers are exported for inspection. Sample IDs are matched case-insensitively; collisions are rejected. Constant expression rows are removed and modules require at least five mapped variable genes. The script uses a Gaussian GSVA kernel for continuous transformed expression. Mean z-scores are standardized across the supplied cohort and are cohort-dependent.

Script 02 additionally expects:

- `results_frozen_analysis_v1/step5_validate_51_results/run_featureSetMerged35_v1/compare_rho_pvalue_custom_set1_set2/compare_5methods_scores_long.csv`: `sample`, `RI`, `model`, `score`.
- `results_frozen_analysis_v1/step7_hallmark_gsva_1125/associations_with_CTA_RI_and_RIscore/subset51_CTA_RI_with_hallmark_GSVA_long.csv`: `sample`, `RI`, `gene_set`, `GSVA_score`. The selected pathways are HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_COMPLEMENT, and HALLMARK_COAGULATION.

## Outputs and interpretation

Outputs go to `results/01_targeted_myeloid_modules` and `results/02_compare_modules`, separate from frozen inputs. Filenames retain the original cohort labels (1125/51) for compatibility; actual sample sizes are calculated from the supplied data.

The code exports effect sizes, available-case sample counts, nominal p-values and BH-adjusted p-values. In script 01, BH adjustment is applied across modules separately within each score method and each correlation/test family. In script 02, it is applied across compared methods. Constant or undersized correlation inputs return NA. Comparisons may use different available-case sample sets; sample counts must be inspected. Negative correlations remain visible in comparison plots.

These exploratory observational associations do not establish causality or out-of-sample predictive performance. The scripts assume one independent observation per patient. No patient-level data or generated outputs are included in the repository, and `.gitignore` excludes local data and results.

Synthetic integration checks passed on 17 November 2026. Original cohort results have not been reproduced by this review. See the accompanying validation report for scope and limitations.

## Method reference

[GSVA package vignette](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html).
