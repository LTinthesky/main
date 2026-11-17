#!/usr/bin/env Rscript

# Train a CTA-derived RNA remodeling index score using concordant RI genes.
# The model is trained on binary RI extremes and evaluated against continuous CTA RI.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glmnet)
})

# -------------------------
# User parameters
# -------------------------
expression_file <- "data/vst_expression_matrix.csv"  # rows = genes, columns = samples; VST/log-like expression
metadata_file <- "data/cta_rna_metadata.csv"
gene_set_file <- "results/02_concordant_gene_sets/concordant_positive_RI_genes.csv"
output_dir <- "results/03_elastic_net_RNA_RI_score"

sample_id_col <- "sample_id"
ri_col <- "CTA_RI"

positive_cutoff <- 1.05
nonpositive_cutoff <- 0.95

alpha_value <- 0.5
lambda_rule <- "lambda.min"  # use "lambda.min" or "lambda.1se"
cv_folds <- 5
set.seed(2026)

source("scripts/utils.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load data
# -------------------------
expr <- read.csv(expression_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
gene_set <- read_csv(gene_set_file, show_col_types = FALSE)

require_columns(gene_set, "gene_id")
prepared <- prepare_score_data(expr, metadata, gene_set$gene_id, sample_id_col,
                               ri_col, positive_cutoff, nonpositive_cutoff)
metadata <- prepared$metadata
x_train <- prepared$x_train
x_all <- prepared$x_all
y_train <- prepared$y
genes <- colnames(x_train)
foldid <- stratified_folds(y_train, cv_folds)
write_csv(tibble(sample_id = rownames(x_train), fold = foldid),
          file.path(output_dir, "CV_fold_assignments.csv"))
# glmnet standardizes within each fit, including CV training folds.
# Gene selection was performed upstream on this cohort: CV is conditional on
# that selected set and is not independent validation of the full pipeline.
cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = alpha_value,
                    foldid = foldid, standardize = TRUE, keep = TRUE)
lambda_value <- cv_fit[[lambda_rule]]
rna_ri_score <- as.numeric(predict(cv_fit, newx = x_all, s = lambda_value, type = "link"))

score_table <- metadata |>
  mutate(
    RNA_RI_score = rna_ri_score,
    RNA_RI_score_z = if (sd(rna_ri_score) > 0) as.numeric(scale(rna_ri_score)) else NA_real_,
    used_for_training = prepared$train
  )

score_cor <- safe_correlation(score_table$RNA_RI_score, score_table[[ri_col]])

coef_table <- as.matrix(coef(cv_fit, s = lambda_value)) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_id") |>
  rename(coefficient = 2) |>
  filter(coefficient != 0)

write_csv(score_table, file.path(output_dir, "RNA_RI_score_CTA_RNA_overlap.csv"))
write_csv(coef_table, file.path(output_dir, "elastic_net_nonzero_coefficients.csv"))
write_csv(
  tibble(
    alpha = alpha_value,
    lambda_rule = lambda_rule,
    lambda = lambda_value,
    n_input_genes = length(genes),
    n_nonzero_coefficients = sum(coef_table$gene_id != "(Intercept)"),
    evaluation = "Descriptive same-cohort association; not independent validation",
    spearman_rho_CTA_RI = unname(score_cor$estimate),
    pvalue_CTA_RI = score_cor$p.value
  ),
  file.path(output_dir, "elastic_net_score_summary.csv")
)

message("Saved RNA RI score results to: ", output_dir)

# The binomial link is a log-odds score, not continuous anatomical RI.
# Preserve the trained model and training score scale for application to new samples.
saveRDS(list(fit = cv_fit, genes = genes, lambda = lambda_value,
             score_center = mean(rna_ri_score), score_scale = sd(rna_ri_score),
             training_sample_ids = rownames(x_train), seed = 2026L),
        file.path(output_dir, "RNA_RI_model.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
