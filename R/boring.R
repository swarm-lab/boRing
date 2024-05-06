.cor <- function(x, y, conf_level) {
  rho <- cor(x, y, method = "spearman")

  if (is.na(conf_level)) {
    c(estimate = rho)
  } else {
    n <- length(x)
    e_fx <- exp(2 * ((0.5 * log((1 + rho) / (1 - rho))) - c(1, -1) *
      (abs(qnorm((1 - conf_level) / 2))) * (1 / sqrt(sum(n) - 3))))
    ci <- (e_fx - 1) / (e_fx + 1)
    c(
      estimate = rho,
      lower_ci = max(ci[1], -1),
      upper_ci = min(ci[2], 1)
    )
  }
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
#' @param conf_level A scalar indicating the confidence level of the required
#'  confidence interval (default: NA, no confidence interval is returned). The
#'  confidence intervals are calculated via z-Transformation.
#'
#' @return A numeric value indicating how boring (i.e., unimodal) the empirical
#'  distribution is. Values close to 1 indicate that the distribution is likely
#'  unimodal.
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
boring <- function(x, w = rep(1, nrow(x)), na_rm = FALSE, conf_level = NA) {
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

  covar <- .wcov(x, w)
  d <- sqrt(.Mahalanobis(x, covar$center, covar$cov))
  ord <- order(d)
  .cor(
    x = -d[ord],
    y = cumsum(w[ord]) / (d[ord]^ncol(x)),
    conf_level = conf_level
  )
}
