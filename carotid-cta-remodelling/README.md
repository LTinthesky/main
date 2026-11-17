# CTA-Linked RNA Remodeling Scores and Plaque Histology

This repository contains simplified code templates for the main analyses used to study CTA-defined vascular remodeling and plaque transcriptomic remodeling in human carotid atherosclerosis.

The scripts are intentionally concise. They show the core analysis logic and main parameters, rather than every exploratory analysis or figure-polishing step.

## Main Analysis Steps

1. Differential expression analysis using DESeq2
   - Binary CTA remodeling: positive RI versus non-positive RI.
   - Continuous CTA remodeling: gene association per 1-SD higher CTA RI.

2. Concordant remodeling gene sets
   - Binary and continuous DESeq2 results are merged.
   - Genes are retained when they show the same direction of association in both analyses.

3. CTA-derived RNA remodeling score
   - Concordant positive-remodeling genes are used to train an elastic-net model.
   - The score is evaluated against continuous CTA RI.

4. Elastic-net robustness
   - The elastic-net score is compared with simpler scores, including mean z-score and PC1.
   - Alpha and lambda settings are varied to assess model robustness and sparsity.

5. Public single-cell RNA-seq localization
   - Remodeling gene sets are scored in a public atherosclerosis single-cell reference using Seurat `AddModuleScore`.

6. Public spatial transcriptomics localization
   - Remodeling gene sets are scored in public Xenium spatial transcriptomic data using `scanpy.tl.score_genes`.

7. Histology associations
   - Estimate unadjusted and age/sex-adjusted logistic associations per unit of an already standardized RNA score.
   - Report complete-case sample sizes, event counts, Wald confidence intervals, BH-adjusted p-values, and failed-model flags.

## Interpretation and review status

These are reviewed analysis templates and adapted analysis scripts, not an exact frozen reproduction of a manuscript. Synthetic integration checks passed on 17 November 2026; the original cohort results have not been reproduced by this review.

The elastic-net model uses RNA expression to distinguish CTA remodeling extremes. Its link-scale score measures transcriptomic similarity to the positive-remodeling group; it is not a measured anatomical RI or a calibrated clinical probability. Gene selection and correlation assessment use the same cohort. The correlations and their nominal p-values are descriptive, and the internal cross-validation only tunes lambda conditional on the selected genes. Independent predictive validation would require an external cohort or nested resampling that repeats gene selection and preprocessing within training folds.

Concordant gene sets use nominal p < 0.01 and an absolute binary log2 fold-change threshold of 0.2. These exploratory criteria do not establish genome-wide FDR significance. The scripts retain adjusted p-values for inspection.

The model and CV folds are saved. All alpha/lambda sensitivity settings are reported, without selecting a winner by its same-cohort correlation. Constant genes are removed using training data. The templates assume one independent sample per patient; when `patient_id` is present, duplicate IDs are rejected. Repeated-patient designs require grouped resampling and an appropriate analysis design.

## Data

Patient-level RNA-seq, CTA, clinical, and histology data are not included because they contain sensitive cohort information. The scripts use placeholder paths under `data/`.

Public single-cell and spatial transcriptomic datasets should be downloaded from their original repositories and placed locally before running the corresponding scripts.

## Running the Scripts

Edit the file paths and column names at the top of each script before running.

Example:

```bash
Rscript scripts/01_deseq2_RI_analysis.R
Rscript scripts/02_concordant_gene_sets.R
Rscript scripts/03_elastic_net_RNA_RI_score.R
Rscript scripts/04_elastic_net_robustness.R
Rscript scripts/05_public_scRNA_module_score.R
python scripts/06_public_xenium_spatial_score.py
# Separate cohort-input analysis; see data/README.md.
Rscript scripts/07_RNA_score_histology_associations.R
```

The scripts write results to the `results/` folder.

## Requirements

R >= 4.2, DESeq2, readr, dplyr, tibble, glmnet, Seurat >= 5, ggplot2, and ragg. Spatial scoring requires Python with scanpy, anndata, numpy, pandas, and scipy. See the validation package's session information for the tested environment; other versions have not been certified.

Public single-cell inputs must have a joined, log-normalized RNA `data` layer. Spatial inputs must contain log-normalized expression in `X` across a suitable control-gene background; the script explicitly ignores `adata.raw`. Scores on a targeted Xenium panel summarize the mapped subset of a signature, not the full transcriptome. Cell-type and regional summaries are descriptive and do not treat cells as independent patient replicates.

Histology inputs use the original frozen cohort layout under `CTA_PROJECT_DIR`, defaulting to `data/cta_project`. Outputs go to the local `results/` directory, including a patient-level model-source table. Keep generated results and data out of the public repository. The existing `.gitignore` excludes them.

## Method references

- [DESeq2 input and design guidance](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [glmnet cross-validation](https://glmnet.stanford.edu/reference/cv.glmnet.html)
- [Seurat AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [Scanpy score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html)
