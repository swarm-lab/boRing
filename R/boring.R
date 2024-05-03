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
boring <- function(x, w = rep(1, nrow(x))) {
  covar <- wcov(x, w)
  d <- Mahalanobis(x, covar$center, covar$cov)
  ord <- order(d)
  cor(d[ord], cumsum(w[ord]) / (d[ord]^2), method = "spearman")^2
}