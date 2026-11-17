# Shared checks for bulk module analysis.
require_columns <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
check_ids <- function(ids, label) {
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop(label, " must be present, nonempty, and unique")
}
safe_correlation <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0)
    return(list(estimate = NA_real_, p.value = NA_real_, n = length(x)))
  test <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  list(estimate = unname(test$estimate), p.value = test$p.value, n = length(x))
}
safe_wilcox <- function(score, group) {
  keep <- is.finite(score) & !is.na(group)
  if (length(unique(group[keep])) != 2L) return(NA_real_)
  suppressWarnings(wilcox.test(score[keep] ~ droplevels(group[keep]), exact = FALSE)$p.value)
}
