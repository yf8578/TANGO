validate_inputs <- function(z, corr = NULL) {
  z <- as.numeric(z)
  if (length(z) < 2) stop("at least two features are required")
  if (any(!is.finite(z))) stop("z contains non-finite values")
  if (is.null(corr)) {
    corr <- diag(length(z))
  } else {
    corr <- as.matrix(corr)
    if (!all(dim(corr) == c(length(z), length(z)))) stop("corr must match z length")
  }
  list(z = z, corr = shrink_correlation(corr))
}

dense_burden_test <- function(z, corr) {
  w <- rep(1, length(z))
  denom <- as.numeric(t(w) %*% corr %*% w)
  denom <- max(denom, 1e-12)
  stat <- as.numeric((sum(z))^2 / denom)
  p <- stats::pchisq(stat, df = 1, lower.tail = FALSE)
  list(stat = stat, p = p)
}

variance_component_test <- function(z, corr) {
  eig <- eigen(corr, symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  mean_q <- sum(eig)
  var_q <- 2 * sum(eig^2)
  df <- 2 * mean_q^2 / var_q
  scale <- var_q / (2 * mean_q)
  stat <- sum(z^2)
  p <- stats::pchisq(stat / scale, df = df, lower.tail = FALSE)
  list(stat = stat, p = p)
}

sparse_minp_test <- function(z, corr) {
  p_feature <- two_sided_z_pvalue(z)
  meff <- li_ji_effective_tests(corr)
  p <- minp_sidak(p_feature, effective_tests = meff)
  list(stat = max(abs(z)), p = p)
}

network_smoothed_test <- function(z, corr, network) {
  s <- normalize_network(network)
  if (!all(dim(s) == dim(corr))) stop("network and corr must have same dimension")
  z2 <- as.numeric(s %*% z)
  c2 <- s %*% corr %*% t(s)
  c2 <- shrink_correlation(c2)
  variance_component_test(z2, c2)
}

#' Run the TANGO adaptive omnibus test
#' @export
tango_test <- function(z, corr = NULL, network = NULL, component_weights = NULL) {
  input <- validate_inputs(z, corr)
  z <- input$z
  corr <- input$corr
  dense <- dense_burden_test(z, corr)
  vc <- variance_component_test(z, corr)
  sparse <- sparse_minp_test(z, corr)
  pvals <- c(dense = dense$p, vc = vc$p, sparse = sparse$p)
  stats <- c(dense = dense$stat, vc = vc$stat, sparse = sparse$stat)
  network_p <- NA_real_
  if (!is.null(network)) {
    net <- network_smoothed_test(z, corr, network)
    pvals <- c(pvals, network = net$p)
    stats <- c(stats, network = net$stat)
    network_p <- net$p
  }
  weights <- NULL
  if (!is.null(component_weights)) {
    weights <- component_weights[names(pvals)]
    weights[is.na(weights)] <- 0
  }
  list(
    pvalue = acat(pvals, weights = weights),
    dense_pvalue = dense$p,
    vc_pvalue = vc$p,
    sparse_pvalue = sparse$p,
    network_pvalue = network_p,
    component_pvalues = pvals,
    statistics = stats,
    n_features = length(z)
  )
}
