
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiCure

[![Build
Status](https://github.com/zhenyu1029/multiCure/workflows/R-CMD-check/badge.svg)](https://github.com/zhenyu1029/multiCure/actions)
[![Codecov test
coverage](https://codecov.io/gh/zhenyu1029/multiCure/branch/main/graph/badge.svg)](https://codecov.io/gh/zhenyu1029/multiCure?branch=main)
[![R-CMD-check](https://github.com/zhenyu1029/multiCure/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhenyu1029/multiCure/actions/workflows/R-CMD-check.yaml)

The goal of **multiCure** is to provide estimation and
hypothesis-testing tools for cure rates in multivariate survival data.
It supports:

- Two‐sample (paired and independent) Wald‐type tests for comparing cure
  probabilities  
- $K$‐sample Wald tests across multiple margins  
- Covariance‐weighted, simple‐average, and pooled Kaplan-Meier
  estimators of the overall cure rate  
- Ying-Wei variance adjustments for the pooled Kaplan-Meier curve

## Installation

You can install the development version of **multiCure** from GitHub
with:

``` r
# install.packages("remotes")  # if you don't already have remotes
# remotes::install_github("zhenyu1029/multiCure")
```

Or using **devtools**:

``` r
# install.packages("devtools")
# devtools::install_github("zhenyu1029/multiCure")
```

## Basic Example

``` r
# devtools::load_all(".")
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(multiCure)
#> Registered S3 method overwritten by 'future':
#>   method               from      
#>   all.equal.connection parallelly
library(copula)

set.seed(123)
# Define data generator
genData2_long<- function(n, tau, u_max1, u_max2, p1, p2, 
                         cure_time1, cure_time2, k1, k2, 
                         missing_prob = 0.1) {
  ## Compute Clayton theta from Kendall's tau
  theta <- copula::iTau(copula::claytonCopula(dim = 2), tau)
  
  ## Draw two independent n×2 copula uniforms:
  ##  - u_event: for potential event-times
  ##  - u_cure: for cure-status indicators
  u_event <- copula::rCopula(n, copula::claytonCopula(theta, dim = 2))
  u_cure  <- copula::rCopula(n, copula::claytonCopula(theta, dim = 2))
  
  ## Pack margin-specific parameters into length-2 vectors
  u_max  <- c(u_max1,      u_max2)
  p      <- c(p1,          p2)
  cure_t <- c(cure_time1,  cure_time2)
  k      <- c(k1,          k2)
  
  ## Loop over the two margins (m = 1, 2) to build each data.frame
  dfs <- lapply(1:2, function(m) {
    ## Inverse transform to get potential event-times for margin m
    t_u <- (
      -log(
        1 - (1 - exp(-cure_t[m]^k[m])) * u_event[, m]
      )
    )^(1 / k[m])
    
    ## Determine cure status for margin m based on p[m]
    cured <- as.integer(u_cure[, m] <= p[m])
    
    ## Assign event-time: Inf if cured, otherwise t_u
    t_all <- ifelse(cured == 1, Inf, t_u)
    
    ## Generate censoring times for margin m in [0, u_max[m]]
    c_time <- runif(n, 0, u_max[m])
    
    ## Calculate observed time and status for margin m
    o_time <- pmin(t_all, c_time)
    d_flag <- as.integer(t_all <= c_time)
    
    ## Introduce missingness at random for margin m
    miss_idx <- runif(n) < missing_prob
    o_time[miss_idx] <- NA
    d_flag[miss_idx] <- NA
    
    ## Create and return a data.frame for margin m
    data.frame(
      id     = seq_len(n),
      margin = factor(as.character(m), levels = c("1", "2")),
      time   = o_time,
      status = d_flag
    )
  })
  
  ## Combine the two margin data.frames into one long-format data.frame
  long_df <- do.call(rbind, dfs)
  rownames(long_df) <- NULL
  return(long_df)
}

## Simulate two‐sample cure‐rate data
set.seed(123)
long_df2 <- genData2_long(
  n = 200,
  tau = 0.5,
  u_max1 = 5,
  u_max2 = 6,
  p1 = 0.4,
  p2 = 0.4,
  cure_time1 = 10,
  cure_time2 = 12,
  k1 = 1.5,
  k2 = 2.0,
  missing_prob = 0.05
)
head(long_df2)
#>   id margin      time status
#> 1  1      1 0.4862635      1
#> 2  2      1 1.8292274      0
#> 3  3      1 0.6063603      0
#> 4  4      1 0.2349684      0
#> 5  5      1 1.3139815      0
#> 6  6      1 0.1295460      1

## Paired two‐sample test
two_paired <- test.2.cure(
  long_df2,
  time,
  status,
  margin,
  type = "paired"
)
print(two_paired)
#> Two-Sample Cure Rate Test (Paired) 
#> ----------------------------------
#> data:  1 vs 2
#> Test Statistic:  -0.1232
#> P-Value:         0.9019
#> Alternative Hypothesis: cure probabilities differ
summary(two_paired)
#> Two-Sample Cure Rate Test (Paired, 95% CI) 
#> ------------------------------------------
#> Estimated cure rate:
#>   1: 0.3948
#>   2: 0.4001
#> 
#> Difference ( 1 - 2 ) = -0.0052
#> 95% CI for difference: [ -0.0886 , 0.0781 ]
#> 
#> Test statistic = -0.1232   p-value = 0.9019
#> Alternative Hypothesis: cure probabilities differ

## Independent two‐sample test
two_indep <- test.2.cure(
  long_df2,
  time,
  status,
  margin,
  type = "indep"
)
print(two_indep)
#> Two-Sample Cure Rate Test (Independent) 
#> ---------------------------------------
#> data:  1 vs 2
#> Test Statistic:  -0.0907
#> P-Value:         0.9277
#> Alternative Hypothesis: cure probabilities differ
summary(two_indep)
#> Two-Sample Cure Rate Test (Independent, 95% CI) 
#> -----------------------------------------------
#> Estimated cure rate:
#>   1: 0.3948
#>   2: 0.4001
#> 
#> Difference ( 1 - 2 ) = -0.0052
#> 95% CI for difference: [ -0.1185 , 0.1080 ]
#> 
#> Test statistic = -0.0907   p-value = 0.9277
#> Alternative Hypothesis: cure probabilities differ

## K‐sample Wald test (here with 2 margins yields same as independent)
wald_k <- test.k.cure(long_df2, time, status, margin)
summary(wald_k)
#> K-Sample Cure Rate Wald Test (Summary) 
#> --------------------------------------
#> Margins:  1, 2 
#> 
#> Test Statistic (Chi^2):  0.0152
#> Degrees of Freedom:      1
#> P-Value:                 0.9019
#> 
#> Alternative Hypothesis: At least one margin's cure probability differs.

## Estimation of cure rates
estimates <- est.cure(long_df2, time, status, margin)
print(estimates)
#> Cure-Rate Estimates 
#> -------------------
#> Weighted estimator:     0.3977 (SE = 0.0348)
#> Simple average:        0.3974 (SE = 0.0349)
#> Pooled KM estimator:   0.3988 (SE = 0.0291)
#> Ying-Wei estimator:    0.3988 (SE = 0.0289)
#> 
#> Per-margin cure rates:
#>   1: 0.3948  (SE = 0.0419)
#>   2: 0.4001  (SE = 0.0397)
```
