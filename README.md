
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
- Covariance‐weighted, simple‐average, and pooled Kaplan–Meier
  estimators of the overall cure rate  
- Ying–Wei variance adjustments for the pooled Kaplan–Meier curve

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

# Define data generator
genData2_long <- function(n, tau, u_max1, u_max2, p1, p2, cure_time1, cure_time2, k1, k2, missing_prob = 0.1) {
  # Calculate the copula parameter theta from Kendall's tau
  theta <- copula::iTau(copula::claytonCopula(dim = 2), tau)
  
  # Generate U(0,1) variates from a Clayton copula for survival times
  u <- copula::rCopula(n, copula::claytonCopula(theta, dim = 2))
  # Generate U(0,1) variates from a Clayton copula for cure status
  u2 <- copula::rCopula(n, copula::claytonCopula(theta, dim = 2))
  
  # --- Margin 1 Calculations ---
  # Inverse transform to get potential event times for non-cured in margin 1
  t1_u <- (-log(1 - (1 - exp(-cure_time1^k1)) * u[, 1]))^(1 / k1)
  # Determine cure status for margin 1 based on p1
  cured1 <- as.integer(u2[, 1] <= p1)
  # Assign event time: Inf if cured, otherwise t1_u
  t1 <- ifelse(cured1 == 1, Inf, t1_u)
  # Generate censoring times for margin 1
  c1 <- runif(n, 0, u_max1)
  # Calculate observed time and status for margin 1
  o1 <- pmin(t1, c1)
  d1 <- as.integer(t1 <= c1)
  
  # Introduce missingness for Margin 1
  miss1 <- runif(n) < missing_prob
  o1[miss1] <- NA
  d1[miss1] <- NA
  
  # --- Margin 2 Calculations ---
  # Inverse transform to get potential event times for non-cured in margin 2
  t2_u <- (-log(1 - (1 - exp(-cure_time2^k2)) * u[, 2]))^(1 / k2)
  # Determine cure status for margin 2 based on p2
  cured2 <- as.integer(u2[, 2] <= p2)
  # Assign event time: Inf if cured, otherwise t2_u
  t2 <- ifelse(cured2 == 1, Inf, t2_u)
  # Generate censoring times for margin 2
  c2 <- runif(n, 0, u_max2)
  # Calculate observed time and status for margin 2
  o2 <- pmin(t2, c2)
  d2 <- as.integer(t2 <= c2)
  
  # Introduce missingness for Margin 2
  miss2 <- runif(n) < missing_prob
  o2[miss2] <- NA
  d2[miss2] <- NA
  
  # Create data frames for each margin
  df_g1 <- data.frame(
    id = seq_len(n),
    margin = factor("1", levels = c("1", "2")),
    time = o1,
    status = d1
  )
  
  df_g2 <- data.frame(
    id = seq_len(n),
    margin = factor("2", levels = c("1", "2")),
    time = o2,
    status = d2
  )
  
  # Combine the two margin data frames into a single long format data frame
  long_df <- rbind(df_g1, df_g2)
  
  return(long_df)
}

# Simulate two‐sample cure‐rate data
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

# Paired two‐sample test
two_paired <- test.2.cure(
  long_df2,
  time,
  status,
  margin,
  type = "paired"
)
print(two_paired)
#> $test_stat
#> [1] -0.1232064
#> 
#> $p_value
#> [1] 0.9019437

# Independent two‐sample test
two_indep <- test.2.cure(
  long_df2,
  time,
  status,
  margin,
  type = "indep"
)
print(two_indep)
#> $test_stat
#> [1] -0.09068138
#> 
#> $p_value
#> [1] 0.9277458

# K‐sample Wald test (here with 2 margins yields same as independent)
wald_k <- test.k.cure(long_df2, time, status, margin)
print(wald_k)
#> $test_stat
#> [1] 0.01517981
#> 
#> $p_value
#> [1] 0.9019437
#> 
#> $df
#> [1] 1
#> 
#> $invertible
#> [1] TRUE

# Estimation of cure rates
estimates <- est.cure(long_df2, time, status, margin)
print(estimates)
#> $weight_vec
#> [1] 0.4506688 0.5493312
#> 
#> $weighted_est
#> [1] 0.3976998
#> 
#> $var_weighted_est
#> [1] 0.001212122
#> 
#> $ave_est
#> [1] 0.3974414
#> 
#> $var_ave_est
#> [1] 0.001216521
#> 
#> $pool_est
#> [1] 0.3987973
#> 
#> $var_pool_est
#> [1] 0.0008448245
#> 
#> $yw_est
#> [1] 0.3987973
#> 
#> $var_yw_est
#> [1] 0.0008333953
#> 
#> $cure_prob_list
#> [1] 0.3948223 0.4000606
#> 
#> $var_cure_prob
#> [1] 0.001757597 0.001579253
#> 
#> $invertible
#> [1] TRUE
```
