#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

//' @name wcov
//' 
//' @title Weighted Covariance Matrix
//' 
//' @description \code{wcov} computes the estimates of the weighted covariance 
//'  matrix and the weighted mean of the data. 
//'
//' @param x A matrix with \eqn{m} columns and \eqn{n} rows, where each column 
//'  represents a different variable and each row a different observation. 
//'
//' @param w A non-negative and non-zero vector of weights for each observation. 
//'  Its length must equal the number of rows of x.
//'
//' @return A list with two components:
//'  \itemize{
//'   \item \code{center}: an estimate for the center (mean) of the data.
//'   \item \code{cov}: the estimated (weighted) covariance matrix.
//'  }
//'
//' @author Simon Garnier, \email{garnier@@njit.edu}
//'
//' @examples
//' m <- matrix(c(rnorm(500, 6), rnorm(500, 11, 3)), ncol = 2)
//' w <- runif(500)
//' wcov(m, w)
//'
// [[Rcpp::export]]
Rcpp::List wcov(Eigen::MatrixXd &x, Eigen::VectorXd &w)
{
    Eigen::VectorXd center(x.cols());
    double ws = w.sum();
    for (int i = 0; i < x.cols(); i++)
    {
        center(i) = (x.col(i).array() * w.array()).sum() / ws;
    }

    int p = x.cols();
    Eigen::VectorXd sqw = (w.array() / ws).cwiseSqrt();
    Eigen::MatrixXd X = x.array().rowwise() - center.transpose().array();
    Eigen::MatrixXd cov = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.transpose() * sqw.asDiagonal());

    return Rcpp::List::create(
        Rcpp::_["center"] = center,
        Rcpp::_["cov"] = cov);
}


//' @name Mahalanobis
//' 
//' @title Mahalanobis Distance
//' 
//' @description \code{Mahalanobis} computes the squared Mahalanobis distance of
//'  all rows in \code{x} and the vector \eqn{\mu}{mu} = \code{center} with 
//'  respect to \eqn{\Sigma}{Sigma} = \code{cov}. This is (for vector \code{x}) 
//'  defined as \deqn{D^2 = (x - \mu)' \Sigma^{-1} (x - \mu)}{D^2 = (x - \mu)' 
//'  \Sigma^-1 (x - \mu)}
//'
//' @param x A matrix with \eqn{m} columns and \eqn{n} rows, where each column 
//'  represents a different variable and each row a different observation.
//'
//' @param center The mean vector of the distribution.
//'
//' @param cov The covariance matrix (\eqn{p \times p}{p x p}) of the 
//'  distribution.
//'
//' @return A matrix with 1 column and \eqn{n} rows reprensenting the estimated
//'  distance for each observation.
//'
//' @author Simon Garnier, \email{garnier@@njit.edu}
//'
//' @examples
//' m <- matrix(c(rnorm(500, 6), rnorm(500, 11, 3)), ncol = 2)
//' w <- runif(500)
//' covar <- wcov(m, w)
//' Mahalanobis(m, covar$center, covar$cov)
//'
// [[Rcpp::export]]
Eigen::MatrixXd Mahalanobis(Eigen::MatrixXd &x, Eigen::VectorXd &center,
                            Eigen::MatrixXd &cov)
{
    Eigen::MatrixXd cx = x.rowwise() - center.transpose();
    Eigen::MatrixXd out = (cx * cov.inverse()).array() * cx.array();
    return out.rowwise().sum();
}
