# Tian Lan | Vascular Research Code

Selected analysis code for carotid plaque transcriptomics, CTA-linked RNA remodeling scores, and associations with plaque histology.

| Project | Methods and scope |
| --- | --- |
| [Bulk RNA-seq Gene Module Scoring](bulk-rna-module-scoring/) | Targeted myeloid GSVA and mean-z module scores; associations with CTA remodeling index; comparisons with other RNA scores. |
| [CTA-Linked RNA Remodeling Scores and Plaque Histology](carotid-cta-remodelling/) | Differential expression, concordant gene sets, elastic-net RNA scores, histology associations, and single-cell/spatial localization. Also available in the [dedicated CTA repository](https://github.com/LTinthesky/carotid_CTA-remodelling). |
| [Carotid Plaque Transcriptomics Across Clinical Presentations](https://github.com/LTinthesky/carotid-plaque-omics-AFX-TIA) | Code for bulk RNA-seq analysis across amaurosis fugax, cerebral TIA, stroke and asymptomatic presentations. |

## Code review and validation

The module-scoring and CTA packages were reviewed and tested with synthetic data on 17 November 2026. See the [validation report](VALIDATION_REPORT.md) and [runnable synthetic checks](validation/).

The CTA package includes simplified analysis templates and adapted scripts. The tests verify code behavior; they do not reproduce original cohort results or establish independent predictive accuracy. The elastic-net score is an RNA signature anchored to imaging-defined remodeling. Same-cohort associations are labeled as descriptive.

No patient-level data, Rubicon application drafts, or generated cohort outputs are included. Each project documents its required local inputs and interpretation limits.

## Earlier smoking analyses

The existing [plaque-type smoking analysis](plaque%20type) and [smoking script](smoker_for_Tian.R) are retained as historical files. They were not covered by the validation of the two packages above.
