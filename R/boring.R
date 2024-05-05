.boring <- function(x, w) {
  covar <- .wcov(x, w)
  d <- sqrt(.Mahalanobis(x, covar$center, covar$cov))
  ord <- order(d)
  cor(d[ord], cumsum(w[ord]) / (d[ord]^ncol(x)), method = "spearman")^2
}

#' @title Multivariate Unimodality Index
#'
#' @description \code{boring} estimates how boring (i.e., unimodal) an empirical
#'  distribution is.
#'
#' @param x A matrix with \eqn{m} columns and \eqn{n} rows, where each column
#'  represents a different variable and each row a different observation.
#'
#' @param w A optional non-negative and non-zero vector of weights for each
#'  observation. Its length must equal the number of rows of x.
#'
#' @param na_rm A logical indicating whether NA values should be stripped before
#'  the computation proceeds (default: FALSE).
#'
#' @param ci A logical indidicating whether a confidence interval should be
#'  estimated by bootstrapping (default: FALSE).
#'
#' @param boot_rep If \code{ci = TRUE}, an integer number of bootstrap
#'  replicates (default: 1000).
#'
#' @param boot_conf A scalar indicating the confidence level of the required
#'  confidence interval (default: 0.95).
#'
#' @param ... Additional parameters to be passed to \code{\link[boot]{boot}}.
#'
#' @return A numeric value indicating how boring (i.e., unimodal) the empirical
#'  distribution is. Values close to 1 indicate that the distribution is likely
#'  unimodal. Values close to 0 indicate that it is likely not.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @examples
#' m1 <- matrix(c(rnorm(500, 6), rnorm(500, 11, 3)), ncol = 2)
#' m2 <- matrix(c(rnorm(500, 2), rnorm(500, 1, 1)), ncol = 2)
#' m3 <- matrix(c(rnorm(500, -13), rnorm(500, -3, 2)), ncol = 2)
#' X <- rbind(m1, m3)
#' boring(X)
#'
#' @export
boring <- function(
    x, w = rep(1, nrow(x)), na_rm = FALSE,
    ci = FALSE, boot_rep = 1000, boot_conf = 0.95, ...) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != length(w)) {
    stop("The length of 'w' does not equal the number of rows of 'x'.")
  }

  if (na_rm) {
    na_ix <- apply(x, 1, function(r) any(is.na(r))) | is.na(w)
    x <- x[!na_ix, , drop = FALSE]
    w <- w[!na_ix]
  }

  if (ci) {
    bt <- boot::boot(x, function(xi, i) .boring(xi[i, , drop = FALSE], w[i]), R = boot_rep, ...)
    bt_ci <- boot::boot.ci(bt, boot_conf, "perc")
    c(estimate = bt_ci$t0, lower = bt_ci$percent[4], upper = bt_ci$percent[5])
  } else {
    c(estimate = .boring(x, w))
  }
}
