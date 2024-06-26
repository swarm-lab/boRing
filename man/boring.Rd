% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boring.R
\name{boring}
\alias{boring}
\title{Multivariate Unimodality Index}
\usage{
boring(x, w = rep(1, nrow(x)), na_rm = FALSE, conf_level = NA)
}
\arguments{
\item{x}{A matrix with \eqn{m} columns and \eqn{n} rows, where each column
represents a different variable and each row a different observation.}

\item{w}{A optional non-negative and non-zero vector of weights for each
observation. Its length must equal the number of rows of x.}

\item{na_rm}{A logical indicating whether NA values should be stripped before
the computation proceeds (default: FALSE).}

\item{conf_level}{A scalar indicating the confidence level of the required
confidence interval (default: NA, no confidence interval is returned). The
confidence intervals are calculated via z-Transformation.}
}
\value{
A numeric value indicating how boring (i.e., unimodal) the empirical
 distribution is. Values close to 1 indicate that the distribution is likely
 unimodal.
}
\description{
\code{boring} estimates how boring (i.e., unimodal) an empirical
 distribution is.
}
\examples{
m1 <- matrix(c(rnorm(500, 6), rnorm(500, 11, 3)), ncol = 2)
m2 <- matrix(c(rnorm(500, 2), rnorm(500, 1, 1)), ncol = 2)
m3 <- matrix(c(rnorm(500, -13), rnorm(500, -3, 2)), ncol = 2)
X <- rbind(m1, m3)
boring(X)

}
\author{
Simon Garnier, \email{garnier@njit.edu}
}
