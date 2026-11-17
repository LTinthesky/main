# Shared input checks. Run entry points from the repository root.
require_columns <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
check_ids <- function(ids, label) {
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop(label, " must be present, nonempty, and unique")
}
check_matrix <- function(x, counts = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x))) stop("Expression must contain finite numeric values")
  check_ids(rownames(x), "Gene IDs")
  check_ids(colnames(x), "Expression sample IDs")
  if (counts && (any(x < 0) || any(x != round(x))))
    stop("DESeq2 requires nonnegative integer raw counts; do not supply normalized expression")
}
align_metadata <- function(metadata, expr, sample_col, ri_col) {
  require_columns(metadata, c(sample_col, ri_col))
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  check_ids(metadata[[sample_col]], "Metadata sample IDs")
  if (!is.numeric(metadata[[ri_col]])) stop("RI must be numeric")
  keep <- metadata[[sample_col]] %in% colnames(expr) & is.finite(metadata[[ri_col]])
  message("Excluded metadata rows without matching expression or finite RI: ", sum(!keep))
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) < 4L || sd(metadata[[ri_col]]) == 0) stop("Insufficient samples or RI variation")
  # These templates assume one independent observation per patient.
  if ("patient_id" %in% names(metadata)) check_ids(as.character(metadata$patient_id), "Patient IDs")
  rownames(metadata) <- metadata[[sample_col]]
  metadata
}
safe_correlation <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0)
    return(list(estimate = NA_real_, p.value = NA_real_, n = length(x)))
  result <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  list(estimate = unname(result$estimate), p.value = result$p.value, n = length(x))
}
stratified_folds <- function(y, requested = 5L, seed = 2026L) {
  sizes <- table(y)
  if (length(sizes) != 2L || min(sizes) < 3L)
    stop("Need at least three independent patients in each RI extreme for stratified CV")
  k <- min(as.integer(requested), min(sizes))
  if (k < 3L) stop("Cross-validation requires at least three folds")
  set.seed(seed)
  foldid <- integer(length(y))
  for (level in names(sizes)) {
    idx <- which(as.character(y) == level)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}
prepare_score_data <- function(expr, metadata, genes, sample_col, ri_col, positive, nonpositive) {
  if (positive <= nonpositive) stop("Positive cutoff must exceed nonpositive cutoff")
  expr <- as.matrix(expr)
  check_matrix(expr)
  metadata <- align_metadata(metadata, expr, sample_col, ri_col)
  metadata$RI_group <- ifelse(metadata[[ri_col]] >= positive, "positive",
                            ifelse(metadata[[ri_col]] <= nonpositive, "non_positive", "mid"))
  genes <- intersect(unique(na.omit(genes)), rownames(expr))
  if (length(genes) < 2L) stop("Need at least two mapped genes for glmnet")
  x_all <- t(expr[genes, metadata[[sample_col]], drop = FALSE])
  train <- metadata$RI_group != "mid"
  if (!any(train)) stop("No RI-extreme training samples")
  variability <- apply(x_all[train, , drop = FALSE], 2, sd)
  keep <- is.finite(variability) & variability > 0
  message("Removed constant training genes: ", sum(!keep))
  x_all <- x_all[, keep, drop = FALSE]
  if (ncol(x_all) < 2L) stop("Need at least two variable training genes")
  list(metadata = metadata, x_all = x_all, x_train = x_all[train, , drop = FALSE],
       train = train, y = as.integer(metadata$RI_group[train] == "positive"))
}
