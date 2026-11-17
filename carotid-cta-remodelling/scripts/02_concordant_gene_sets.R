#!/usr/bin/env Rscript

# Define concordant remodeling-associated genes from binary and continuous RI DESeq2.
# Edit paths and thresholds below before running.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# -------------------------
# User parameters
# -------------------------
binary_file <- "results/01_deseq2_RI_analysis/DESeq2_binary_positive_vs_nonpositive_RI.csv"
continuous_file <- "results/01_deseq2_RI_analysis/DESeq2_continuous_RI_scaled.csv"
annotation_file <- "data/gene_annotation.csv"  # columns: gene_id, symbol
output_dir <- "results/02_concordant_gene_sets"

p_cutoff <- 0.01
log2fc_cutoff <- 0.2

source("scripts/utils.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load DE tables
# -------------------------
binary <- read_csv(binary_file, show_col_types = FALSE) |>
  rename(
    log2FC_binary = log2FoldChange,
    pvalue_binary = pvalue,
    padj_binary = padj
  )

continuous <- read_csv(continuous_file, show_col_types = FALSE) |>
  rename(
    beta_RI_scaled = log2FoldChange,
    pvalue_continuous = pvalue,
    padj_continuous = padj
  )

annotation <- read_csv(annotation_file, show_col_types = FALSE)

require_columns(binary, c("gene_id", "log2FC_binary", "pvalue_binary", "padj_binary"))
require_columns(continuous, c("gene_id", "beta_RI_scaled", "pvalue_continuous", "padj_continuous"))
require_columns(annotation, c("gene_id", "symbol"))
check_ids(binary$gene_id, "Binary-result gene IDs")
check_ids(continuous$gene_id, "Continuous-result gene IDs")
check_ids(annotation$gene_id, "Annotation gene IDs")
# Nominal p-values define an exploratory feature set, not FDR-significant discoveries.
merged <- binary |>
  select(gene_id, log2FC_binary, pvalue_binary, padj_binary) |>
  inner_join(
    continuous |> select(gene_id, beta_RI_scaled, pvalue_continuous, padj_continuous),
    by = "gene_id"
  ) |>
  left_join(annotation, by = "gene_id") |>
  mutate(
    same_direction = sign(log2FC_binary) == sign(beta_RI_scaled),
    concordant_positive = pvalue_binary < p_cutoff &
      pvalue_continuous < p_cutoff &
      log2FC_binary > log2fc_cutoff &
      beta_RI_scaled > 0,
    concordant_negative = pvalue_binary < p_cutoff &
      pvalue_continuous < p_cutoff &
      log2FC_binary < -log2fc_cutoff &
      beta_RI_scaled < 0
  )

write_csv(merged, file.path(output_dir, "binary_continuous_RI_DESeq2_merged.csv"))

positive_genes <- merged |>
  filter(concordant_positive) |>
  arrange(pvalue_binary)

negative_genes <- merged |>
  filter(concordant_negative) |>
  arrange(pvalue_binary)

write_csv(positive_genes, file.path(output_dir, "concordant_positive_RI_genes.csv"))
write_csv(negative_genes, file.path(output_dir, "concordant_negative_RI_genes.csv"))

writeLines(na.omit(unique(positive_genes$symbol)), file.path(output_dir, "concordant_positive_RI_symbols.txt"))
writeLines(na.omit(unique(negative_genes$symbol)), file.path(output_dir, "concordant_negative_RI_symbols.txt"))

summary_table <- tibble(
  threshold = paste0("p_binary < ", p_cutoff, ", p_continuous < ", p_cutoff, ", |log2FC_binary| > ", log2fc_cutoff),
  n_positive_rows = nrow(positive_genes),
  n_positive_symbols = length(unique(na.omit(positive_genes$symbol))),
  n_negative_rows = nrow(negative_genes),
  n_negative_symbols = length(unique(na.omit(negative_genes$symbol)))
)

write_csv(summary_table, file.path(output_dir, "concordant_gene_set_summary.csv"))

message("Saved concordant gene sets to: ", output_dir)
