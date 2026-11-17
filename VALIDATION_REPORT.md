# Research code review — 17 November 2026

## Material Passport

- Mode: validate; code review and synthetic integration checks.
- Material: two locally prepared research-code packages, derived from saved Codex project files.
- Verification status: ANALYZED. Synthetic execution verified; original cohort results not reproduced.
- Data used: generated synthetic expression, metadata, histology and public annotation identifiers only.
- Delivery: reviewed copies; original saved project scripts and cohort files unchanged.

## Findings and corrections

| Area | Finding | Correction |
|---|---|---|
| Sample alignment | Sample and gene identifiers were not consistently checked; DESeq2 metadata row names were implicit. | Validate IDs, align by explicit sample IDs, set row names, reject duplicate patients when a patient_id field is supplied. |
| DESeq2 inputs | Counts were silently rounded; numeric batch encoding and missing/nonfinite RI could alter or break models. | Reject noninteger/nonfinite/negative counts, treat batch categorically, report exclusions, reject confounded or saturated designs. |
| Elastic-net scaling | Manually scaling the complete training set before CV reused scale estimates across folds. | Use glmnet internal standardization within each fit. |
| CV comparisons | Alpha settings used different random folds. | Reuse deterministic class-stratified folds, adapting fold count to class size. |
| Constant genes/scores | Zero variance could produce nonfinite values or failed correlations. | Filter constant genes on training data; return NA for undefined correlations. |
| Model reporting | Nonzero gene count subtracted an intercept regardless of whether it was nonzero; the sensitivity summary selected the largest same-cohort correlation. | Exclude the intercept explicitly; report every setting without choosing a winner by the reported correlation. |
| Score reuse | Fitted model and score-scale parameters were not saved. | Save fitted model, gene order, lambda, score mean/SD and training IDs. |
| PC1 comparator | PCA was fitted to all scored samples. | Fit PCA using training samples, apply those loadings to all samples, orient without consulting RI. |
| Histology | Invalid fits, one-level outcomes, separation and unsupported category labels were not explicitly handled. | Add input/category checks, sample and event counts, model-status flags, and omit failed fits from the plot. |
| Histology figure | Fixed axis limits could hide large/small intervals. | Use an unconstrained logarithmic OR axis. |
| Bulk modules | Constant genes could break mean-z scores; duplicate case-normalized IDs and insufficient correlations were unchecked. | Validate IDs and overlap, remove constant expression, use safe correlation/Wilcoxon handling. |
| Multiple testing | Bulk-module tests lacked adjusted p-values. | Add BH adjustment with explicitly documented test families; keep nominal p-values. |
| Comparison plot | X-axis began at zero, hiding negative correlations. | Show negative correlations; avoid a hard-coded sample-size claim. |
| Single-cell scoring | Missing normalized RNA layers, annotation columns and empty gene sets were unchecked. | Require Seurat v5 joined RNA data layer and annotations, explicit seed; empty sets receive numeric NA. |
| Spatial scoring | Scanpy could implicitly use adata.raw although gene matching used X; missing scores were object-valued. | Explicitly score X with a fixed seed, deduplicate genes, require finite input, reject obvious integer-count input, use numeric NaN for empty sets. |
| Output locations | Adapted scripts wrote into frozen original-analysis output paths. | Write generated results into each package's ignored local results directory. |

## Executed checks

- All R files parsed; Python script parsed.
- Synthetic end-to-end execution of DESeq2 binary/continuous analyses, concordant-gene selection, elastic-net score fitting and all 10 alpha/lambda settings.
- Shuffled metadata retained correct sample alignment; reloaded model predictions matched exported scores to numerical tolerance < 1e-10.
- Checks for duplicate metadata, noninteger counts, undersized CV classes, constant expression and reproducible fold assignments.
- Seurat module scoring on synthetic log-normalized expression with finite scores and preserved cell count.
- Histology adjusted logistic beta matched an independently specified direct glm to < 1e-10; zero-event and perfectly separated outcomes were flagged.
- Bulk GSVA and mean-z scoring executed on synthetic expression with public human gene identifiers; module comparisons retained negative correlations and BH columns.
- Spatial scores matched a direct Scanpy reference call using X, despite a different raw matrix; duplicate genes, an empty signature, regional/cell-type summaries and raw-count rejection were checked.

The R and spatial execution logs and runnable synthetic-check scripts are in validation/. R session information was captured. Scanpy 1.11.5 was used. GSVA emitted a package-build-version warning but executed successfully in the tested R environment. Warnings about synthetic Seurat background-gene naming and sparse-matrix coercion did not affect the score assertions.

## Remaining scientific limits

1. The six pre-existing CTA core scripts are simplified templates. They are not certified exact reproductions of the frozen manuscript workflow, sample exclusions, ComBat-seq preprocessing, original feature-set construction or published figures. No biological effect, original sample count or manuscript p-value is certified by these synthetic tests.
2. Concordant genes are selected using nominal p < 0.01. This is an exploratory feature-selection rule, not genome-wide FDR-controlled discovery. Both binary and continuous selection use the same patients and are not independent replications.
3. Elastic-net correlations use the same cohort involved in feature selection and training. Lambda CV conditional on a preselected gene set is not independent validation of the complete pipeline. Full predictive evaluation requires an external cohort or nested resampling that repeats preprocessing, gene selection and tuning within training folds. No such performance is claimed here.
4. The score is an RNA log-odds signature anchored to imaging-defined remodeling. It is not anatomical RI, a validated imaging-only model, or a calibrated clinical probability.
5. The curated neutrophil/myeloid lists require documented biological provenance before being described as a named published signature. Their exact definitions and overlap are available in code and output; genes may overlap each other and the CTA feature set.
6. Histology models assume the input RNA score is already standardized. Complete-case selection, donor independence, label definitions, linearity of continuous predictors on the logit scale, residual confounding and transportability require cohort-level assessment. Separation checks are conservative diagnostics, not a proof that every possible separation issue is excluded.
7. Single-cell/spatial analyses summarize scores descriptively. Normalization provenance and donor annotations must be verified for real inputs. Cells do not supply independent patient replication; targeted-panel overlap limits biological coverage.
8. Source-data outputs can contain patient-level information. Only code, documentation and synthetic tests are prepared for upload; raw data, model files, results, application drafts and local provenance manifests are excluded.

## Statistical fallacy coverage: 11/11 assessed for applicability

| Check | Assessment in this code review |
|---|---|
| Simpson's paradox | Cannot resolve without cohort/subgroup results; adjustment does not establish direction consistency across groups. |
| Ecological fallacy | Patient-level models and descriptive cell summaries have distinct units; no patient-level inference from cell counts is claimed. |
| Berkson's paradox | Endarterectomy-cohort selection limits generalizability; not quantified by synthetic tests. |
| Collider bias | Causal adjustment assumptions cannot be established from age/sex covariates alone. |
| Base-rate neglect | No clinical probability or diagnostic accuracy claim is made. |
| Regression to the mean | Extreme RI groups define score training; no longitudinal improvement claim is made. |
| Survivorship bias | Complete-case exclusions and follow-up selection need cohort-level review; no survival model is included in this package. |
| Look-elsewhere effect | Nominal gene selection is identified as exploratory; all sensitivity settings and adjusted bulk/histology p-values are retained. |
| Garden of forking paths | The reviewed copies retain explicit thresholds/settings and report all tested elastic-net configurations. No preregistration is claimed. |
| Correlation versus causation | Outputs are labeled as observational associations, not causal effects. |
| Reverse causality | Cross-sectional associations cannot establish temporal or causal direction. |

## Primary technical references

- [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [glmnet CV reference](https://glmnet.stanford.edu/reference/cv.glmnet.html)
- [GSVA vignette](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html)
- [Seurat AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [Scanpy score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html)
