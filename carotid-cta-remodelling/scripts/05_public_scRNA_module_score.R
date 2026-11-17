#!/usr/bin/env Rscript

# Public atherosclerosis scRNA-seq module score analysis.
# Requires a Seurat object with published cell-type annotations.

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
})

# -------------------------
# User parameters
# -------------------------
seurat_rds_file <- "data/public_atherosclerosis_scRNA_seurat_object.rds"
positive_gene_file <- "results/02_concordant_gene_sets/concordant_positive_RI_symbols.txt"
negative_gene_file <- "results/02_concordant_gene_sets/concordant_negative_RI_symbols.txt"
output_dir <- "results/05_public_scRNA_module_score"

celltype_col <- "cell_type_level2"
assay_name <- "RNA"
set.seed(2026)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load object and gene sets
# -------------------------
obj <- readRDS(seurat_rds_file)
if (!assay_name %in% Assays(obj)) stop("Missing assay: ", assay_name)
if (!celltype_col %in% names(obj[[]])) stop("Missing cell-type annotation")
if (length(Layers(obj[[assay_name]], search = "^data$")) != 1L)
  stop("Provide one joined, log-normalized RNA data layer before scoring")
DefaultAssay(obj) <- assay_name

positive_genes <- readLines(positive_gene_file)
negative_genes <- readLines(negative_gene_file)

positive_genes <- intersect(unique(trimws(positive_genes)), rownames(obj))
negative_genes <- intersect(unique(trimws(negative_genes)), rownames(obj))

writeLines(positive_genes, file.path(output_dir, "mapped_positive_RI_genes.txt"))
writeLines(negative_genes, file.path(output_dir, "mapped_negative_RI_genes.txt"))

# -------------------------
# Seurat AddModuleScore
# -------------------------
# Empty gene sets cannot be scored; retain numeric NA and report mapped genes.
for (entry in list(list(genes = positive_genes, name = "Positive_RI_module"),
                   list(genes = negative_genes, name = "Negative_RI_module"))) {
  if (length(entry$genes)) {
    obj <- AddModuleScore(object = obj, features = list(entry$genes), name = entry$name,
                          assay = assay_name, slot = "data", seed = 2026)
  } else {
    warning("No mapped genes for ", entry$name)
    obj[[paste0(entry$name, "1")]] <- rep(NA_real_, ncol(obj))
  }
}

positive_score_col <- "Positive_RI_module1"
negative_score_col <- "Negative_RI_module1"

score_df <- obj@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  select(
    cell_id,
    cell_type = all_of(celltype_col),
    positive_RI_module_score = all_of(positive_score_col),
    negative_RI_module_score = all_of(negative_score_col)
  )

summary_df <- score_df |>
  group_by(cell_type) |>
  summarize(
    n_cells = n(),
    positive_median = median(positive_RI_module_score, na.rm = TRUE),
    positive_mean = mean(positive_RI_module_score, na.rm = TRUE),
    negative_median = median(negative_RI_module_score, na.rm = TRUE),
    negative_mean = mean(negative_RI_module_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(positive_median))

write_csv(score_df, file.path(output_dir, "public_scRNA_RI_module_scores_per_cell.csv"))
write_csv(summary_df, file.path(output_dir, "public_scRNA_RI_module_scores_by_celltype.csv"))

message("Saved public scRNA module score results to: ", output_dir)
