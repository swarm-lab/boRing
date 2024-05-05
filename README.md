
# boRing <a href="https://swarm-lab.github.io/boRing/"><img src="man/figures/logo.png" align="right" height="138" alt="boRing website" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/boRing)](https://CRAN.R-project.org/package=boRing)
[![R-CMD-check](https://github.com/swarm-lab/boRing/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/swarm-lab/boRing/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/swarm-lab/boRing/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/swarm-lab/boRing/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Description

The goal of `boRing` is to provide an easy metric to determine whether a 
multivariate empirical distribution is likely unimodal (that is, a boring 
distribution) or not. The "boringness" index is based on the idea that the 
density of observations should decrease monotonically with the the distance to 
the center of mass of an unimodal empirical distribution.

The computation of the "boringness" index is done as follows: 

1. We compute the center of mass and the covariance matrix of the empirical 
distribution (weighted observations are allowed in `boRing`).

2. From the center of mass and the covariance matrix, we compute the Mahalanobis
distance of each observation to the center of mass of the distribution.

3. We order the observations based on their Mahalanobis distance (from closest
to furthest away).

4. We compute the squared Spearman correlation between the Mahalanobis distances 
of the ordered observations and the weighted density of the observations up to 
each computed distance. If the correlation is close to 1, it indicates that the
density of observations decreases monotonically with the distance to the center
of mass of the distribution and, hence, that the distribution is likely unimodal
(or, boring). 

---

## Installation

At this time, `boRing` is not yet available on CRAN. 

You can install the development version of boRing from 
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("swarm-lab/boRing")
```

---

## Performance tip

`boRing` depends on [`RcppEigen`](https://github.com/RcppCore/RcppEigen) for 
computing efficiently the weighted covariance matrix of the empirical 
distribution as well as the Mahalanobis distances of the observations. Provided
that your compiler allows it, a significant performance increase can be obtained
by adding the `-O3` optimization flag to the `CXXFLAGS` line in your local
`.R/Makevars` file. If you don't know where your local `.R/Makevars` file is 
located, install the [`usethis`](https://usethis.r-lib.org/) package and then 
run `usethis::edit_r_makevars()`.

---

## Example

``` r
library(boRing)

m1 <- matrix(c(rnorm(500, 6), rnorm(500, 11, 3)), ncol = 2)
m2 <- matrix(c(rnorm(500, 2), rnorm(500, 1, 1)), ncol = 2)
m3 <- matrix(c(rnorm(500, -13), rnorm(500, -3, 2)), ncol = 2)

X <- rbind(m1, m2, m3)
plot(X, asp = 1)
boring(X)
```
