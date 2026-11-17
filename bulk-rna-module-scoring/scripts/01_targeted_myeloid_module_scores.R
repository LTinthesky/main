suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GSVA)
  library(BiocParallel)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

source("scripts/utils.R")
set.seed(2026)
project_dir <- Sys.getenv("CTA_PROJECT_DIR", unset = "data/cta_project")
vst_path <- file.path(
  project_dir,
  "results_frozen_analysis_v1/step1_preprocess_1125/vst_1125_allGenes_from_combatseq.rds"
)
feature_set_path <- file.path(
  project_dir,
  "results_frozen_analysis_v1/step3_feature_sets/feature_set_fromMergedTable_p001_sameDir_allDirections_ensembl.csv"
)
cta_file <- file.path(
  project_dir,
  "results_frozen_analysis_v1/step6_apply_1125_results/run_v1/cta_assoc_subset51_from1125score/subset51_scoresFrom1125_with_CTA.csv"
)
out_dir <- "results/01_targeted_myeloid_modules"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading VST object: ", vst_path)
vst <- readRDS(vst_path)
expr <- assay(vst)
storage.mode(expr) <- "double"
if (any(!is.finite(expr))) stop("VST expression must be finite")
check_ids(rownames(expr), "Expression gene IDs")
check_ids(tolower(colnames(expr)), "Case-normalized expression sample IDs")
if (ncol(expr) < 3L) stop("Need at least three samples")
variable <- apply(expr, 1, sd) > 0
message("Removing constant genes: ", sum(!variable))
expr <- expr[variable, , drop = FALSE]
expr_genes <- rownames(expr)

module_symbols <- list(
  Literature_canonical_neutrophil = c(
    "CSF3R", "CXCR1", "CXCR2", "FCGR3B", "CEACAM8", "S100A8", "S100A9",
    "S100A12", "MPO", "ELANE", "PRTN3", "CTSG", "LCN2", "MMP8", "MMP9",
    "OLFM4", "CAMP", "DEFA1", "DEFA3", "AZU1", "LTF", "BPI", "SELL",
    "ITGAM", "ITGB2", "FPR1", "FPR2", "VNN2", "MCEMP1"
  ),
  Single_cell_plaque_myeloid = c(
    "LST1", "AIF1", "TYROBP", "FCER1G", "LYZ", "CD14", "FCGR3A", "MS4A7",
    "CD68", "C1QA", "C1QB", "C1QC", "FCN1", "VCAN", "S100A8", "S100A9",
    "S100A12", "IL1B", "NLRP3", "TREM1", "CCL2", "CCL3", "CCL4", "CCL8",
    "CXCL2", "PLIN2", "TREM2", "APOE", "LPL", "LGALS3", "SPP1", "GPNMB",
    "CTSB", "CTSL", "MMP9"
  )
)

feature_set <- read_csv(feature_set_path, show_col_types = FALSE)
require_columns(feature_set, c("direction", "symbol", "gene"))
positive_set35 <- feature_set |>
  filter(direction == "positive_remodeling") |>
  transmute(
    module = "CTA_RI_positive_set35",
    symbol = symbol,
    ensembl_gene = gene,
    source = "Frozen CTA RI positive-remodeling feature-set rows"
  )

map_symbols_to_ensembl <- function(module_name, symbols) {
  symbols <- unique(symbols)
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL")
  )
  mapped <- mapped |>
    filter(!is.na(ENSEMBL)) |>
    distinct(SYMBOL, ENSEMBL) |>
    transmute(
      module = module_name,
      symbol = SYMBOL,
      ensembl_gene = ENSEMBL,
      source = dplyr::case_when(
        module_name == "Literature_canonical_neutrophil" ~
          "Curated canonical neutrophil abundance/activation markers from standard immunology literature",
        module_name == "Single_cell_plaque_myeloid" ~
          "Curated plaque scRNA myeloid markers covering monocyte, inflammatory macrophage, foamy/TREM2 macrophage and phagocyte programs",
        TRUE ~ "Custom module"
      )
    )
  missing <- setdiff(symbols, mapped$symbol)
  if (length(missing) > 0) {
    mapped <- bind_rows(
      mapped,
      tibble(
        module = module_name,
        symbol = missing,
        ensembl_gene = NA_character_,
        source = paste("Curated module:", module_name)
      )
    )
  }
  mapped
}

custom_mapped <- bind_rows(lapply(names(module_symbols), function(nm) {
  map_symbols_to_ensembl(nm, module_symbols[[nm]])
}))

module_gene_map <- bind_rows(custom_mapped, positive_set35) |>
  mutate(
    in_vst = !is.na(ensembl_gene) & ensembl_gene %in% expr_genes,
    symbol = as.character(symbol)
  ) |>
  arrange(module, symbol, ensembl_gene)

write_csv(module_gene_map, file.path(out_dir, "targeted_module_gene_map_full.csv"))
write_csv(
  filter(module_gene_map, in_vst),
  file.path(out_dir, "targeted_module_gene_map_overlap_vst.csv")
)

overlap_report <- module_gene_map |>
  group_by(module, source) |>
  summarize(
    n_input_symbols = n_distinct(symbol),
    n_input_ensembl_rows = sum(!is.na(ensembl_gene)),
    n_overlap_ensembl_rows = n_distinct(ensembl_gene[in_vst]),
    n_overlap_symbols = n_distinct(symbol[in_vst]),
    missing_symbols = paste(sort(unique(symbol[!in_vst])), collapse = ";"),
    .groups = "drop"
  )
write_csv(overlap_report, file.path(out_dir, "targeted_module_overlap_report.csv"))

gene_sets <- split(module_gene_map$ensembl_gene[module_gene_map$in_vst], module_gene_map$module[module_gene_map$in_vst])
gene_sets <- lapply(gene_sets, function(x) sort(unique(na.omit(x))))
gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 5]

if (!length(gene_sets)) stop("No modules have at least five variable mapped genes")
message("Running GSVA for ", length(gene_sets), " targeted modules")
param <- gsvaParam(
  exprData = expr,
  geneSets = gene_sets,
  kcdf = "Gaussian",
  minSize = 5,
  maxSize = Inf
)
gsva_scores <- gsva(param, BPPARAM = SerialParam())

zexpr <- t(scale(t(expr)))
z_scores <- do.call(rbind, lapply(names(gene_sets), function(module_name) {
  genes <- gene_sets[[module_name]]
  colMeans(zexpr[genes, , drop = FALSE], na.rm = TRUE)
}))
rownames(z_scores) <- names(gene_sets)

write_csv(
  data.frame(module = rownames(gsva_scores), gsva_scores, check.names = FALSE),
  file.path(out_dir, "targeted_module_GSVA_scores_1125.csv")
)
write_csv(
  data.frame(module = rownames(z_scores), z_scores, check.names = FALSE),
  file.path(out_dir, "targeted_module_meanZ_scores_1125.csv")
)
saveRDS(gsva_scores, file.path(out_dir, "targeted_module_GSVA_scores_1125.rds"))
saveRDS(z_scores, file.path(out_dir, "targeted_module_meanZ_scores_1125.rds"))

scores_long <- bind_rows(
  as.data.frame(gsva_scores) |>
    tibble::rownames_to_column("module") |>
    pivot_longer(-module, names_to = "sample", values_to = "module_score") |>
    mutate(score_method = "GSVA"),
  as.data.frame(z_scores) |>
    tibble::rownames_to_column("module") |>
    pivot_longer(-module, names_to = "sample", values_to = "module_score") |>
    mutate(score_method = "mean_gene_zscore")
) |>
  mutate(sample = tolower(sample))
write_csv(scores_long, file.path(out_dir, "targeted_module_scores_1125_long.csv"))

cta <- read_csv(cta_file, show_col_types = FALSE) |>
  mutate(
    sample = tolower(sample),
    RI_binary = factor(RI_binary, levels = c("Nonpositive_remodeling", "Positive_remodeling"))
  )

require_columns(cta, c("sample", "RI", "RI_score", "RI_binary"))
check_ids(cta$sample, "CTA sample IDs after case normalization")
if (!is.numeric(cta$RI) || !is.numeric(cta$RI_score)) stop("RI and RI_score must be numeric")
cta_merged <- inner_join(cta, scores_long, by = "sample")
if (!nrow(cta_merged)) stop("No matching CTA and expression samples")
write_csv(cta_merged, file.path(out_dir, "subset51_CTA_RI_with_targeted_module_scores_long.csv"))

cor_stats <- cta_merged |>
  group_by(module, score_method) |>
  summarize(
    n = sum(is.finite(RI) & is.finite(module_score)),
    spearman_rho_CTA_RI = safe_correlation(RI, module_score)$estimate,
    spearman_p_CTA_RI = safe_correlation(RI, module_score)$p.value,
    pearson_r_CTA_RI = safe_correlation(RI, module_score, "pearson")$estimate,
    pearson_p_CTA_RI = safe_correlation(RI, module_score, "pearson")$p.value,
    spearman_rho_transcriptomic_RI_score = safe_correlation(RI_score, module_score)$estimate,
    spearman_p_transcriptomic_RI_score = safe_correlation(RI_score, module_score)$p.value,
    .groups = "drop"
  ) |>
  arrange(score_method, desc(abs(spearman_rho_CTA_RI)))
cor_stats <- cor_stats |>
  group_by(score_method) |>
  mutate(across(starts_with("spearman_p_"), ~ p.adjust(.x, method = "BH"), .names = "{.col}_BH"),
         pearson_p_CTA_RI_BH = p.adjust(pearson_p_CTA_RI, method = "BH")) |>
  ungroup()
write_csv(cor_stats, file.path(out_dir, "subset51_targeted_module_vs_CTA_RI_stats.csv"))

box_stats <- cta_merged |>
  filter(!is.na(RI_binary)) |>
  group_by(module, score_method) |>
  summarize(
    n_nonpositive = sum(RI_binary == "Nonpositive_remodeling" & !is.na(module_score)),
    n_positive = sum(RI_binary == "Positive_remodeling" & !is.na(module_score)),
    wilcox_p_RI_binary = safe_wilcox(module_score, RI_binary),
    delta_median_positive_minus_nonpositive =
      median(module_score[RI_binary == "Positive_remodeling"], na.rm = TRUE) -
      median(module_score[RI_binary == "Nonpositive_remodeling"], na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(score_method, wilcox_p_RI_binary)
box_stats <- box_stats |>
  group_by(score_method) |>
  mutate(wilcox_p_RI_binary_BH = p.adjust(wilcox_p_RI_binary, method = "BH")) |>
  ungroup()
write_csv(box_stats, file.path(out_dir, "subset51_targeted_module_by_RI_binary_stats.csv"))

module_labels <- c(
  Literature_canonical_neutrophil = "Canonical neutrophil",
  Single_cell_plaque_myeloid = "Plaque scRNA myeloid",
  CTA_RI_positive_set35 = "CTA RI positive set35"
)

plot_data <- cta_merged |>
  mutate(
    module_label = recode(module, !!!module_labels, .default = module),
    module_label = factor(module_label, levels = module_labels),
    score_method = factor(score_method, levels = c("GSVA", "mean_gene_zscore"))
  )

p1 <- ggplot(plot_data, aes(x = RI, y = module_score)) +
  geom_point(aes(color = RI_binary), size = 2.2, alpha = 0.85, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.45, na.rm = TRUE) +
  facet_grid(score_method ~ module_label, scales = "free_y") +
  scale_color_manual(values = c(Nonpositive_remodeling = "#2166ac", Positive_remodeling = "#b2182b"), na.value = "grey55") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Targeted myeloid module scores versus CTA remodeling index",
    subtitle = "CTA overlap subset; GSVA primary, mean z-score sensitivity",
    x = "CTA RI",
    y = "Module score",
    color = "CTA RI binary"
  )
ggsave(file.path(out_dir, "subset51_targeted_module_scores_vs_CTA_RI_scatter.pdf"), p1, width = 11, height = 6.5)
ggsave(file.path(out_dir, "subset51_targeted_module_scores_vs_CTA_RI_scatter.png"), p1, width = 11, height = 6.5, dpi = 300)

p2 <- ggplot(filter(plot_data, !is.na(RI_binary)), aes(x = RI_binary, y = module_score, fill = RI_binary)) +
  geom_boxplot(width = 0.58, alpha = 0.86, outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.75, na.rm = TRUE) +
  facet_grid(score_method ~ module_label, scales = "free_y") +
  scale_fill_manual(values = c(Nonpositive_remodeling = "#2166ac", Positive_remodeling = "#b2182b")) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    title = "Targeted myeloid module scores by CTA RI binary group",
    subtitle = "Binary CTA labels are shown as sensitivity analysis",
    x = "",
    y = "Module score"
  )
ggsave(file.path(out_dir, "subset51_targeted_module_scores_by_RI_binary_boxplot.pdf"), p2, width = 11, height = 6.5)
ggsave(file.path(out_dir, "subset51_targeted_module_scores_by_RI_binary_boxplot.png"), p2, width = 11, height = 6.5, dpi = 300)

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "targeted_module_sessionInfo.txt"))

message("Saved targeted module outputs to: ", out_dir)
