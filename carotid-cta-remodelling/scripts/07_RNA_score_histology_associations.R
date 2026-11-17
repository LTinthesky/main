source("scripts/utils.R")
project_dir <- Sys.getenv("CTA_PROJECT_DIR", unset = "data/cta_project")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

input_file <- file.path(project_dir, "results_frozen_analysis_v1/step6_apply_1125_results/run_featureSetMerged35_v1/clinical_assoc_set35_whole1093/clinical1093_with_set35_scores.csv")
output_dir <- "results/07_histology_associations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

require_columns(d, c("RI_score_set35_z", "Age", "Gender", "hist_type", "Macrophages.bin",
                     "Calc.bin", "Collagen.bin", "SMC.bin", "IPH.bin"))
for (field in c("RI_score_set35_z", "Age")) {
  if (!is.numeric(d[[field]]) || any(!is.finite(d[[field]]) & !is.na(d[[field]])))
    stop(field, " must be numeric and finite when observed")
}
expected <- list(Gender = c("female", "male"), hist_type = c("atheromatous", "fibrous", "fibroatheromatous"),
                 Macrophages.bin = c("moderate/heavy", "no/minor"), Calc.bin = c("moderate/heavy", "no/minor"),
                 Collagen.bin = c("moderate/heavy", "no/minor"), SMC.bin = c("moderate/heavy", "no/minor"),
                 IPH.bin = c("yes", "no"))
for (field in names(expected)) {
  if (any(!is.na(d[[field]]) & !d[[field]] %in% expected[[field]]))
    stop("Unexpected category in ", field, "; explicitly map labels before analysis")
}
if ("patient_id" %in% names(d)) check_ids(as.character(d$patient_id), "Patient IDs")
analysis <- data.frame(
  score_z = as.numeric(d$RI_score_set35_z),
  age = as.numeric(d$Age),
  sex = factor(d$Gender, levels = c("female", "male")),
  plaque_atheromatous = ifelse(d$hist_type == "atheromatous", 1,
                               ifelse(d$hist_type == "fibrous", 0, NA)),
  macrophage_high = ifelse(d$Macrophages.bin == "moderate/heavy", 1,
                           ifelse(d$Macrophages.bin == "no/minor", 0, NA)),
  calcification_high = ifelse(d$Calc.bin == "moderate/heavy", 1,
                              ifelse(d$Calc.bin == "no/minor", 0, NA)),
  collagen_high = ifelse(d$Collagen.bin == "moderate/heavy", 1,
                         ifelse(d$Collagen.bin == "no/minor", 0, NA)),
  smc_high = ifelse(d$SMC.bin == "moderate/heavy", 1,
                    ifelse(d$SMC.bin == "no/minor", 0, NA)),
  iph = ifelse(d$IPH.bin == "yes", 1, ifelse(d$IPH.bin == "no", 0, NA))
)

outcomes <- c(
  plaque_atheromatous = "Atheromatous plaque",
  macrophage_high = "Macrophage-rich plaque",
  calcification_high = "Moderate/heavy calcification",
  collagen_high = "Collagen-rich plaque",
  smc_high = "SMC-rich plaque",
  iph = "Intraplaque hemorrhage"
)

fit_one <- function(outcome, label, adjusted = FALSE) {
  vars <- c(outcome, "score_z", if (adjusted) c("age", "sex"))
  dat <- analysis[, vars, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  result <- data.frame(
    outcome = outcome, label = label,
    model = if (adjusted) "Age- and sex-adjusted" else "Unadjusted",
    n = nrow(dat), events = sum(dat[[outcome]] == 1), nonevents = sum(dat[[outcome]] == 0),
    beta = NA_real_, se = NA_real_, odds_ratio = NA_real_, conf_low = NA_real_,
    conf_high = NA_real_, p_value = NA_real_, status = "not_fitted", stringsAsFactors = FALSE)
  form <- reformulate(c("score_z", if (adjusted) c("age", "sex")), response = outcome)
  if (nrow(dat) < 5L || length(unique(dat[[outcome]])) < 2L || length(unique(dat$score_z)) < 2L) {
    result$status <- "insufficient_data_or_variation"
    return(result)
  }
  if (adjusted && nlevels(droplevels(dat$sex)) < 2L) {
    result$status <- "sex_has_one_level"
    return(result)
  }
  mm <- model.matrix(form, data = dat)
  if (qr(mm)$rank < ncol(mm) || nrow(mm) <= ncol(mm)) {
    result$status <- "singular_or_saturated_design"
    return(result)
  }
  issues <- character()
  fit <- withCallingHandlers(glm(form, data = dat, family = binomial()),
    warning = function(w) { issues <<- c(issues, conditionMessage(w)); invokeRestart("muffleWarning") })
  sm <- summary(fit)$coefficients
  if (!fit$converged || length(issues) || any(fitted(fit) < 1e-8 | fitted(fit) > 1 - 1e-8)) {
    result$status <- "nonconvergence_or_possible_separation"
    return(result)
  }
  beta <- sm["score_z", "Estimate"]; se <- sm["score_z", "Std. Error"]
  if (!all(is.finite(c(beta, se))) || se <= 0) {
    result$status <- "invalid_standard_error"
    return(result)
  }
  result$beta <- beta; result$se <- se
  result$odds_ratio <- exp(beta); result$conf_low <- exp(beta - 1.96 * se)
  result$conf_high <- exp(beta + 1.96 * se); result$p_value <- sm["score_z", "Pr(>|z|)"]
  result$status <- "ok"
  result
}

results <- bind_rows(lapply(names(outcomes), function(v) {
  bind_rows(
    fit_one(v, outcomes[[v]], FALSE),
    fit_one(v, outcomes[[v]], TRUE)
  )
})) %>%
  group_by(model) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    estimate_ci = sprintf("%.2f (%.2f-%.2f)", odds_ratio, conf_low, conf_high),
    p_display = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    fdr_display = ifelse(p_fdr < 0.001, "<0.001", sprintf("%.3f", p_fdr))
  )

write.csv(results,
          file.path(output_dir, "RNA_remodeling_score_histology_logistic_models.csv"),
          row.names = FALSE)
write.csv(analysis,
          file.path(output_dir, "RNA_remodeling_score_histology_model_source_data.csv"),
          row.names = FALSE)

label_order <- rev(unname(outcomes))
plot_df <- results %>%
  filter(status == "ok") %>%
  mutate(
    label = factor(label, levels = label_order),
    model = factor(model, levels = c("Unadjusted", "Age- and sex-adjusted"))
  )

pd <- position_dodge(width = 0.48)
p <- ggplot(plot_df, aes(x = odds_ratio, y = label,
                         linetype = model, shape = model)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high),
                width = 0, linewidth = 0.6, position = pd,
                orientation = "y", colour = "black") +
  geom_point(size = 2.2, stroke = 0.65, position = pd, colour = "black") +
  scale_x_log10() +
  scale_linetype_manual(values = c("Unadjusted" = "dashed",
                                   "Age- and sex-adjusted" = "solid")) +
  scale_shape_manual(values = c("Unadjusted" = 1,
                                "Age- and sex-adjusted" = 16)) +
  labs(
    x = "Odds ratio per 1-SD increase in RNA remodeling score",
    y = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 8.5, base_family = "Helvetica") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.4, colour = "black"),
    axis.text.x = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8.5, colour = "black"),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 7)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(12, "pt"),
    plot.margin = margin(5, 9, 5, 5),
    panel.grid = element_blank()
  )

pdf_file <- file.path(output_dir, "RNA_remodeling_score_histology_logistic_forestplot.pdf")
grDevices::pdf(pdf_file, width = 145 / 25.4, height = 92 / 25.4,
               family = "Helvetica", useDingbats = FALSE)
print(p)
dev.off()

# R-generated raster used only for visual quality assurance.
qa_file <- file.path(output_dir, "RNA_remodeling_score_histology_logistic_forestplot_QA.png")
ragg::agg_png(qa_file, width = 145 / 25.4, height = 92 / 25.4,
              units = "in", res = 180)
print(p)
dev.off()

writeLines(c(
  "Purpose: estimate associations between RNA remodeling score and histological plaque features.",
  "The input score is assumed already standardized; no additional scaling is applied.",
  "Models with insufficient data, singular design or possible separation are flagged and omitted from the plot.",
  "Models: binary logistic regression with RNA remodeling score standardized per 1 SD.",
  "Unadjusted: histology ~ RNA score.",
  "Adjusted: histology ~ RNA score + age + sex.",
  "Reference outcomes: fibrous plaque; no/minor macrophages, calcification, collagen or SMCs; and no intraplaque hemorrhage.",
  "Fibroatheromatous plaques were excluded from the plaque-phenotype model.",
  "Intervals: Wald 95% confidence intervals.",
  "Multiple testing: Benjamini-Hochberg FDR calculated separately within each model family.",
  paste0("Input: ", input_file)
), file.path(output_dir, "RNA_remodeling_score_histology_logistic_forestplot_QA_notes.txt"))

cat("PDF:", pdf_file, "\n")
cat("QA:", qa_file, "\n")
print(results[, c("label", "model", "n", "events", "odds_ratio", "conf_low", "conf_high", "p_value", "p_fdr")])

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
