#!/usr/bin/env Rscript

# Core DESeq2 analysis for CTA-defined carotid remodeling.
# Edit the paths and column names below before running.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
})

# -------------------------
# User parameters
# -------------------------
count_file <- "data/raw_counts_matrix.csv"      # rows = genes, columns = samples
metadata_file <- "data/cta_rna_metadata.csv"    # one row per sample
output_dir <- "results/01_deseq2_RI_analysis"

sample_id_col <- "sample_id"
ri_col <- "CTA_RI"
batch_col <- "batch"

positive_cutoff <- 1.05
nonpositive_cutoff <- 0.95

source("scripts/utils.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load data
# -------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)

counts <- as.matrix(counts)
check_matrix(counts, counts = TRUE)
metadata <- align_metadata(metadata, counts, sample_id_col, ri_col)
require_columns(metadata, batch_col)
if (anyNA(metadata[[batch_col]])) stop("Batch must not contain missing values")
metadata[[batch_col]] <- factor(metadata[[batch_col]])
if (nlevels(metadata[[batch_col]]) < 2L) stop("This batch-adjusted template requires at least two batches")
metadata <- metadata |>
  mutate(
    RI_group = factor(case_when(
      .data[[ri_col]] >= positive_cutoff ~ "positive",
      .data[[ri_col]] <= nonpositive_cutoff ~ "non_positive",
      TRUE ~ NA_character_
    ), levels = c("non_positive", "positive")),
    RI_scaled = as.numeric(scale(.data[[ri_col]]))
  )
rownames(metadata) <- metadata[[sample_id_col]]
counts <- counts[, metadata[[sample_id_col]], drop = FALSE]
if (any(colSums(counts) == 0)) stop("Raw-count samples must have nonzero library size")

# -------------------------
# Binary RI DESeq2 model
# Positive RI >= 1.05 versus non-positive RI <= 0.95.
# Mid-zone samples are excluded.
# -------------------------
metadata_binary <- metadata |>
  filter(RI_group %in% c("non_positive", "positive")) |>
  droplevels()

if (nlevels(metadata_binary$RI_group) != 2L || any(table(metadata_binary$RI_group) < 2L))
  stop("Need at least two samples in both RI extremes")
rownames(metadata_binary) <- metadata_binary[[sample_id_col]]
for (term in c("RI_group", "RI_scaled")) {
  dat <- if (term == "RI_group") metadata_binary else metadata
  mm <- model.matrix(reformulate(c(batch_col, term)), data = dat)
  if (qr(mm)$rank < ncol(mm) || nrow(mm) <= ncol(mm))
    stop("Confounded or saturated DESeq2 design: ", term)
}
counts_binary <- counts[, metadata_binary[[sample_id_col]], drop = FALSE]

dds_binary <- DESeqDataSetFromMatrix(
  countData = counts_binary,
  colData = metadata_binary,
  design = as.formula(paste("~", batch_col, "+ RI_group"))
)

dds_binary <- dds_binary[rowSums(counts(dds_binary)) >= 10, ]
dds_binary <- DESeq(dds_binary)

res_binary <- results(dds_binary, contrast = c("RI_group", "positive", "non_positive")) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_id") |>
  arrange(pvalue)

write_csv(res_binary, file.path(output_dir, "DESeq2_binary_positive_vs_nonpositive_RI.csv"))

# -------------------------
# Continuous RI DESeq2 model
# All CTA-RNA overlap samples are retained.
# Effect is per 1-SD higher CTA RI.
# -------------------------
dds_continuous <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = as.formula(paste("~", batch_col, "+ RI_scaled"))
)

dds_continuous <- dds_continuous[rowSums(counts(dds_continuous)) >= 10, ]
dds_continuous <- DESeq(dds_continuous)

res_continuous <- results(dds_continuous, name = "RI_scaled") |>
  as.data.frame() |>
  tibble::rownames_to_column("gene_id") |>
  arrange(pvalue)

write_csv(res_continuous, file.path(output_dir, "DESeq2_continuous_RI_scaled.csv"))

message("Saved DESeq2 results to: ", output_dir)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
