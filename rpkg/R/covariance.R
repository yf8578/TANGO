symmetrize <- function(mat) {
  (mat + t(mat)) / 2
}

#' Regularize a correlation matrix
#' @export
shrink_correlation <- function(corr, shrinkage = 0.05, min_eigen = 1e-6) {
  c <- as.matrix(corr)
  if (nrow(c) != ncol(c)) stop("corr must be square")
  k <- nrow(c)
  c <- symmetrize(c)
  c <- (1 - shrinkage) * c + shrinkage * diag(k)
  ed <- eigen(c, symmetric = TRUE)
  vals <- pmax(ed$values, min_eigen)
  c2 <- ed$vectors %*% diag(vals, length(vals)) %*% t(ed$vectors)
  c2 <- symmetrize(c2)
  d <- sqrt(pmax(diag(c2), min_eigen))
  out <- c2 / outer(d, d)
  diag(out) <- 1
  symmetrize(out)
}

#' Estimate correlation from null SNP Z scores
#' @export
estimate_null_correlation <- function(null_z, shrinkage = 0.05) {
  z <- as.matrix(null_z)
  if (nrow(z) < 3) stop("at least three rows are required")
  c <- stats::cor(z, use = "pairwise.complete.obs")
  c[!is.finite(c)] <- 0
  diag(c) <- 1
  shrink_correlation(c, shrinkage = shrinkage)
}

normalize_network <- function(network, add_identity = TRUE) {
  w <- as.matrix(network)
  if (nrow(w) != ncol(w)) stop("network must be square")
  w <- pmax(symmetrize(w), 0)
  if (add_identity) w <- w + diag(nrow(w))
  deg <- rowSums(w)
  deg[deg <= 0] <- 1
  d <- 1 / sqrt(deg)
  (d * w) * rep(d, each = nrow(w))
}
