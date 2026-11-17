#!/usr/bin/env Rscript

# Elastic-net robustness analysis.
# Compares elastic-net scores with simple z-score and PC1 scores,
# and checks alpha/lambda sensitivity and model sparsity.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glmnet)
})

# -------------------------
# User parameters
# -------------------------
expression_file <- "data/vst_expression_matrix.csv"
metadata_file <- "data/cta_rna_metadata.csv"
gene_set_file <- "results/02_concordant_gene_sets/concordant_positive_RI_genes.csv"
output_dir <- "results/04_elastic_net_robustness"

sample_id_col <- "sample_id"
ri_col <- "CTA_RI"

positive_cutoff <- 1.05
nonpositive_cutoff <- 0.95

alpha_grid <- c(0, 0.25, 0.5, 0.75, 1)
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
# -------------------------
# Simple score definitions
# -------------------------
training_center <- colMeans(x_train)
training_scale <- apply(x_train, 2, sd)
x_train_scaled <- scale(x_train, center = training_center, scale = training_scale)
x_all_scaled <- scale(x_all, center = training_center, scale = training_scale)
z_score_signature <- rowMeans(x_all_scaled)
pc <- prcomp(x_train_scaled, center = FALSE, scale. = FALSE)
pc1_score <- as.numeric(predict(pc, newdata = x_all_scaled)[, 1])
# Choose a deterministic loading orientation without using the RI outcome.
if (sum(pc$rotation[, 1]) < 0) pc1_score <- -pc1_score

cor_z <- safe_correlation(z_score_signature, metadata[[ri_col]])
cor_pc1 <- safe_correlation(pc1_score, metadata[[ri_col]])

score_method_results <- tibble(
  method = c("Mean z-score signature", "PC1 score"),
  spearman_rho = c(unname(cor_z$estimate), unname(cor_pc1$estimate)),
  pvalue = c(cor_z$p.value, cor_pc1$p.value),
  n_genes = length(genes)
)

# -------------------------
# Elastic-net alpha/lambda sensitivity
# -------------------------
robustness_results <- data.frame()
score_table <- metadata[, c(sample_id_col, ri_col, "RI_group")]

for (alpha_value in alpha_grid) {
  cv_fit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_value,
    foldid = foldid,
    standardize = TRUE
  )

  for (lambda_rule in c("lambda.min", "lambda.1se")) {
    lambda_value <- if (lambda_rule == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min
    score <- as.numeric(predict(cv_fit, newx = x_all, s = lambda_value, type = "link"))
    ct <- safe_correlation(score, metadata[[ri_col]])

    coef_table <- as.matrix(coef(cv_fit, s = lambda_value)) |>
      as.data.frame() |>
      tibble::rownames_to_column("gene_id") |>
      rename(coefficient = 2)

    n_nonzero <- sum(coef_table$coefficient != 0 & coef_table$gene_id != "(Intercept)")

    robustness_results <- bind_rows(
      robustness_results,
      tibble(
        alpha = alpha_value,
        lambda_rule = lambda_rule,
        lambda = lambda_value,
        spearman_rho = unname(ct$estimate),
        pvalue = ct$p.value,
        n_nonzero_genes = n_nonzero
      )
    )

    score_col <- paste0("elastic_net_alpha", alpha_value, "_", lambda_rule)
    score_table[[score_col]] <- score
  }
}

# Report every setting; do not choose a winner by its training-cohort correlation.
elastic_comparison <- robustness_results |>
  mutate(method = paste0("Elastic net, alpha=", alpha, ", ", lambda_rule)) |>
  select(method, spearman_rho, pvalue, n_genes = n_nonzero_genes)
score_method_results <- bind_rows(score_method_results, elastic_comparison) |>
  mutate(evaluation = "Descriptive same-cohort association; not independent validation")

write_csv(score_method_results, file.path(output_dir, "score_method_comparison.csv"))
write_csv(robustness_results, file.path(output_dir, "elastic_net_alpha_lambda_robustness.csv"))
write_csv(score_table, file.path(output_dir, "elastic_net_robustness_scores_by_sample.csv"))

message("Saved elastic-net robustness results to: ", output_dir)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
