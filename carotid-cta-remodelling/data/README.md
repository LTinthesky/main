# Data Folder

This folder is a placeholder for local input files.

Sensitive patient-level data are not included in this repository.

Expected local input files:

- `raw_counts_matrix.csv`: raw RNA-seq count matrix, with genes as rows and samples as columns.
- `vst_expression_matrix.csv`: transformed expression matrix for scoring, with genes as rows and samples as columns.
- `cta_rna_metadata.csv`: sample metadata with CTA RI, batch, and sample identifiers.
- `gene_annotation.csv`: gene annotation table with `gene_id` and `symbol` columns.
- `public_atherosclerosis_scRNA_seurat_object.rds`: public scRNA-seq Seurat object.
- `public_xenium_panel2.h5ad`: public Xenium spatial transcriptomic AnnData object.

The scripts use these names as examples. They can be changed at the top of each script.

Metadata columns: `sample_id` (unique), `CTA_RI` (numeric), and `batch` (categorical). This template requires at least two batches and non-confounded, non-saturated designs. Missing or nonfinite RI rows and unmatched metadata rows are excluded with a count. Missing batches and duplicate IDs cause a clear error. Raw counts must be finite, nonnegative integers. Each row is assumed to represent a different patient.

For script 07, set `CTA_PROJECT_DIR` to the local cohort root. The input relative path is:
`results_frozen_analysis_v1/step6_apply_1125_results/run_featureSetMerged35_v1/clinical_assoc_set35_whole1093/clinical1093_with_set35_scores.csv`.

Required columns: `RI_score_set35_z` (already standardized), `Age`, `Gender` (`female`/`male`), `hist_type` (`fibrous`/`atheromatous`/`fibroatheromatous`), `Macrophages.bin`, `Calc.bin`, `Collagen.bin`, `SMC.bin` (`no/minor`/`moderate/heavy`), and `IPH.bin` (`no`/`yes`). Missing categories are allowed; unknown labels fail rather than being silently discarded. Fibroatheromatous plaques are excluded from the binary atheromatous-versus-fibrous model. Missingness is handled separately for each model.
