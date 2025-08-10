
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LaCDS

<!-- badges: start -->

<!-- badges: end -->

In survival analysis, accurate identification of latent classes is
essential to effectively account for potential hidden population
heterogeneity. In response to this challenge, we introduce the Latent
Class Discrete Survival (LaCDS) model. LaCDS employs a finite-mixture
model structure within the context of the discrete failure time model
and implements the expectation-maximization algorithm for efficient
optimization.

## Installation

You can install the development version of LaCDS from
[GitHub](https://github.com/Qinmengge/LaCDS) with:

``` r
install.packages("LaCDS")
```

## Example

This is a basic example which shows you how to fit the LaCDS model:

``` r
library(LaCDS)
LaCDS(t,ind,X,G=2,tol=0.0001,max_iter = 1000, min_iter=1,alpha=0.5)
```

### Input parameters

t: a integer vector of the event time values starting from 1. ind: a
integer vector of the event indicator taking values 0 and 1 with 0
corresponding to censoring and 1 corresponding to event. X: a numeric
matrix of the covariates of interest. G: a integer scalar of the maximum
number of subgroups to be fitted. tol: a numeric scalar of the model
convergence criteria defined as the l2 norm of covariate estimates
between iterations. max_iter: a integer scalar of the maximum rounds of
iterations. min_iter: a integer scalar of the minimum rounds of
iterations. alpha: a numeric scalar of the step size control parameter.

### Output models

A list of fitted models corresponding to each g smaller or equal to G.
