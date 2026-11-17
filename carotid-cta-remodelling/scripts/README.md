# Script Overview

These scripts are simple analysis templates. Each script has editable parameters at the top, followed by the main analysis steps.

## 01_deseq2_RI_analysis.R

Runs the main DESeq2 analyses:

- positive remodeling versus non-positive remodeling
- continuous CTA RI association

Main parameters:

- `positive_cutoff <- 1.05`
- `nonpositive_cutoff <- 0.95`
- `batch_col <- "batch"`

## 02_concordant_gene_sets.R

Merges binary and continuous DESeq2 results and defines concordant genes.

Main parameters:

- `p_cutoff <- 0.01`
- `log2fc_cutoff <- 0.2`

## 03_elastic_net_RNA_RI_score.R

Trains the CTA-derived RNA remodeling score.

Main parameters:

- `alpha_value <- 0.5`
- `lambda_rule <- "lambda.min"`
- `cv_folds <- 5`

## 04_elastic_net_robustness.R

Checks whether the RNA remodeling score depends on one model setting.

It compares:

- mean z-score signature
- PC1 score
- elastic-net scores across alpha and lambda settings

Main parameters:

- `alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)`
- `cv_folds <- 5`

## 05_public_scRNA_module_score.R

Scores remodeling gene sets in a public Seurat single-cell object using `AddModuleScore`.

Main parameter:

- `celltype_col <- "cell_type_level2"`

## 06_public_xenium_spatial_score.py

Scores remodeling gene sets in public Xenium `.h5ad` data using `scanpy.tl.score_genes`.

Main parameters:

- `region_col = "region"`
- `celltype_col = "cell_type"`

## 07_RNA_score_histology_associations.R

Adapted from the original histology-association analysis. Uses the frozen cohort input layout documented in `data/README.md`; outputs go to `results/07_histology_associations`. Failed or potentially separated logistic models are flagged rather than presented with unreliable Wald intervals. Failed models remain in the exported results table and are omitted from the plot.

## Review changes

`utils.R` supplies shared input, sample-alignment, constant-gene and stratified-CV checks. Scripts 03 and 04 use glmnet's internal standardization, identical stratified folds across settings, safe handling of constant scores, and explicit same-cohort interpretation. Script 03 saves its fitted model. Script 04 reports every model setting. The nominal gene-set thresholds in script 02 are exploratory and have not been changed to a different scientific selection rule.
