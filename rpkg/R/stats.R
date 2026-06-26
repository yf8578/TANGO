safe_pvalue <- function(p, eps = 1e-15) {
  if (!is.finite(p)) stop("p-value must be finite")
  min(max(p, eps), 1 - eps)
}

#' Combine p-values using ACAT
#' @export
acat <- function(pvalues, weights = NULL) {
  p <- as.numeric(pvalues)
  if (length(p) == 0) stop("pvalues cannot be empty")
  p <- vapply(p, safe_pvalue, numeric(1))
  if (is.null(weights)) {
    w <- rep(1 / length(p), length(p))
  } else {
    w <- as.numeric(weights)
    if (length(w) != length(p)) stop("weights must have same length as pvalues")
    if (any(w < 0) || any(!is.finite(w))) stop("weights must be finite and non-negative")
    if (sum(w) <= 0) stop("at least one weight must be positive")
    w <- w / sum(w)
  }
  cauchy_stat <- sum(w * tan((0.5 - p) * pi))
  safe_pvalue(0.5 - atan(cauchy_stat) / pi)
}

two_sided_z_pvalue <- function(z) {
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}

li_ji_effective_tests <- function(corr) {
  eig <- eigen(corr, symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  sum(pmin(eig, 1))
}

minp_sidak <- function(pvalues, effective_tests = NULL) {
  p <- as.numeric(pvalues)
  if (length(p) == 0) stop("pvalues cannot be empty")
  pmin <- safe_pvalue(min(p))
  m <- if (is.null(effective_tests)) length(p) else effective_tests
  m <- max(as.numeric(m), 1)
  safe_pvalue(1 - (1 - pmin)^m)
}
