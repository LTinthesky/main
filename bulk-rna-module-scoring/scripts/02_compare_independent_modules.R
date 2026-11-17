suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("scripts/utils.R")
project_dir <- Sys.getenv("CTA_PROJECT_DIR", unset = "data/cta_project")
base_dir <- file.path(
  project_dir,
  "results_frozen_analysis_v1/step5_validate_51_results/run_featureSetMerged35_v1/compare_rho_pvalue_custom_set1_set2"
)
targeted_dir <- "results/01_targeted_myeloid_modules"
hallmark_assoc_dir <- file.path(
  project_dir,
  "results_frozen_analysis_v1/step7_hallmark_gsva_1125/associations_with_CTA_RI_and_RIscore"
)

existing_stats_file <- file.path(base_dir, "compare_5methods_spearman_rho_pvalue_improved.csv")
existing_scores_file <- file.path(base_dir, "compare_5methods_scores_long.csv")
targeted_scores_file <- file.path(targeted_dir, "subset51_CTA_RI_with_targeted_module_scores_long.csv")
hallmark_scores_file <- file.path(hallmark_assoc_dir, "subset51_CTA_RI_with_hallmark_GSVA_long.csv")

out_dir <- "results/02_compare_modules"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- file.path(out_dir, "compare_8methods_with_independent_modules")

calc_stats <- function(df) {
  df |>
    group_by(model) |>
    summarize(
      rho = safe_correlation(RI, score)$estimate,
      pvalue = safe_correlation(RI, score)$p.value,
      n = sum(is.finite(RI) & is.finite(score)),
      .groups = "drop"
    ) |>
    mutate(
      pvalue_BH = p.adjust(pvalue, method = "BH"),
      neglog10p = -log10(pmax(pvalue, .Machine$double.xmin)),
      text_label = paste0("rho = ", sprintf("%.3f", rho), ", p = ", format.pval(pvalue, digits = 2, eps = 1e-99)),
      p_label = paste0("p=", format.pval(pvalue, digits = 2, eps = 1e-99))
    )
}

existing_scores <- read_csv(existing_scores_file, show_col_types = FALSE) |>
  transmute(
    sample = tolower(sample),
    RI = RI,
    model = model,
    score = score,
    comparison_group = "CTA-derived score"
  )

targeted_scores <- read_csv(targeted_scores_file, show_col_types = FALSE) |>
  filter(
    score_method == "GSVA",
    module %in% c("Literature_canonical_neutrophil", "Single_cell_plaque_myeloid")
  ) |>
  transmute(
    sample = tolower(sample),
    RI = RI,
    model = recode(
      module,
      Literature_canonical_neutrophil = "Canonical neutrophil module GSVA",
      Single_cell_plaque_myeloid = "Plaque myeloid module GSVA"
    ),
    score = module_score,
    comparison_group = "Independent biology module"
  )

hallmark_scores <- read_csv(hallmark_scores_file, show_col_types = FALSE) |>
  filter(gene_set %in% c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION"
  )) |>
  group_by(sample, RI, gene_set) |>
  summarize(GSVA_score = mean(GSVA_score, na.rm = TRUE), .groups = "drop") |>
  group_by(gene_set) |>
  mutate(GSVA_score_z = as.numeric(scale(GSVA_score))) |>
  ungroup() |>
  group_by(sample, RI) |>
  summarize(
    score = mean(GSVA_score_z, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    sample = tolower(sample),
    model = "Hallmark inflammation/complement/coagulation GSVA",
    comparison_group = "Independent pathway module"
  ) |>
  select(sample, RI, model, score, comparison_group)

expanded_scores <- bind_rows(existing_scores, targeted_scores, hallmark_scores)
# Repeated rows would overweight patients and invalidate the correlation p-value.
if (anyDuplicated(expanded_scores[c("sample", "model")]))
  stop("Duplicate sample/model rows; resolve repeated measurements before comparison")
expanded_stats <- calc_stats(expanded_scores)

method_labels <- c(
  "Simple Z-score Signature" = "Simple Z-score Signature",
  "PC1 Score" = "PC1 Score",
  "Elastic Net Score (main)" = "Elastic Net Score\n(concordant binary+continuous DE, p<0.01)",
  "Elastic Net Score (Set1: continuous p<0.01 up)" = "Elastic Net Score\n(continuous RI DE, p<0.01 up)",
  "Elastic Net Score (Set2: binary p<0.01 up)" = "Elastic Net Score\n(binary RI DE, p<0.01 up)",
  "Canonical neutrophil module GSVA" = "Canonical neutrophil\nGSVA",
  "Plaque myeloid module GSVA" = "Plaque scRNA myeloid\nGSVA",
  "Hallmark inflammation/complement/coagulation GSVA" = "Hallmark inflammation+\ncomplement+coagulation GSVA"
)

expanded_stats <- expanded_stats |>
  mutate(
    method_label = recode(model, !!!method_labels, .default = model),
    comparison_group = case_when(
      model %in% c(
        "Canonical neutrophil module GSVA",
        "Plaque myeloid module GSVA"
      ) ~ "Independent biology module",
      model == "Hallmark inflammation/complement/coagulation GSVA" ~ "Independent pathway module",
      TRUE ~ "CTA-derived score"
    ),
    model_order = match(model, names(method_labels))
  ) |>
  arrange(model_order) |>
  select(model, rho, pvalue, pvalue_BH, n, neglog10p, text_label, method_label, p_label, comparison_group)

write_csv(expanded_scores, paste0(out_prefix, "_scores_long.csv"))
write_csv(expanded_stats, paste0(out_prefix, "_spearman_rho_pvalue.csv"))

plot_df <- expanded_stats |>
  mutate(method_label = factor(method_label, levels = rev(method_label)))

p <- ggplot(plot_df, aes(x = rho, y = method_label, color = comparison_group, size = neglog10p)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey65") +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = p_label), hjust = -0.08, size = 3.1, color = "black") +
  scale_color_manual(values = c(
    "CTA-derived score" = "#2166ac",
    "Independent biology module" = "#b2182b",
    "Independent pathway module" = "#4d9221"
  )) +
  scale_size_continuous(range = c(2.5, 7.5), name = "-log10 p") +
  coord_cartesian(xlim = c(-1, 1.15), clip = "off") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 45, 8, 8)
  ) +
  labs(
    title = "CTA RI association across score definitions",
    subtitle = "Descriptive same-cohort associations; sample size reported separately for each method",
    x = "Spearman rho with CTA RI",
    y = "",
    color = "Score type"
  )

ggsave(paste0(out_prefix, "_spearman_rho_pvalue.pdf"), p, width = 10.5, height = 5.8)
ggsave(paste0(out_prefix, "_spearman_rho_pvalue.png"), p, width = 10.5, height = 5.8, dpi = 300)

print(expanded_stats)
message("Saved expanded comparison to: ", base_dir)
