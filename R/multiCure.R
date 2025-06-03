#' @importFrom stats    na.omit pnorm pchisq qnorm
#' @importFrom dplyr    %>% group_by group_split pull
#' @importFrom rlang    enquo !!
#' @importFrom purrr    map map_int
#' @importFrom survival survfit Surv cluster
#' @importFrom prodlim  prodlim Hist
NULL



#' @title Solve a linear system with invertible check
#'
#' @description
#' Attempts to solve \(A x = vec\).  If \(A\) is singular, returns a vector of NAs and sets `invertible = FALSE`.
#'
#' @param A A square numeric matrix.
#' @param vec A numeric vector of length `ncol(A)`.
#'
#' @return A list with components:
#' \describe{
#'   \item{solve}{Numeric vector giving the solution or NAs.}
#'   \item{invertible}{Logical, `TRUE` if `solve(A, vec)` succeeded.}
#' }
#'
#' @keywords internal
solve.system <- function(A, vec) {
  tryCatch(
    list(solve = solve(A, vec), invertible = TRUE),
    error = function(e) list(solve = rep(NA, ncol(A)), invertible = FALSE)
  )
}



#' @title Construct adjacent-difference contrast matrix
#'
#' @description
#' Builds a \\((k-1)\\times k\\) matrix with rows \\((1,-1,0,\\dots)\\), \\((0,1,-1,0,\\dots)\\), etc.
#'
#' @param k Integer number of margins.
#'
#' @return A \\((k-1)\\times k\\) numeric matrix.
#'
#' @keywords internal
construct.C <- function(k) {
  C <- matrix(0, nrow = k - 1, ncol = k)
  for (i in seq_len(k - 1)) {
    C[i, i]     <-  1
    C[i, i + 1] <- -1
  }
  C
}



#' @title Influence function for Kaplan-Meier cure probability
#'
#' @description
#' Computes the influence-function values \\(\\xi_j\\) for a single margin at a specified time point.
#'
#' @param times1 Numeric vector of observed times (may include NAs).
#' @param status1 Integer vector (0/1) indicating censoring/event, same length as `times1`.
#' @param surv_prob_at_last_event_time1 Numeric: Kaplan-Meier survival estimate at `time_point1`.
#' @param time_point1 Numeric time at which the survival is evaluated.
#'
#' @return Numeric vector of the same length as `times1`, giving influence contributions (zeros for NAs).
#'
#' @keywords internal
est.inf.fun <- function(times1, status1, surv_prob_at_last_event_time1, time_point1) {
  times_nonna <- na.omit(times1)
  risk_nonna  <- sapply(times_nonna, function(t) sum(times_nonna >= t))

  risk_set <- rep(NA_integer_, length(times1))
  risk_set[!is.na(times1)] <- risk_nonna

  n  <- length(times1)
  xi <- numeric(n)

  for (j in seq_len(n)) {
    if (!is.na(times1[j])) {
      first_term <- 0
      if (status1[j] == 1 && times1[j] <= time_point1) {
        first_term <- 1 / (risk_set[j] / n)
      }
      second_term <- 0
      idx <- which(times1 <= time_point1 & !is.na(times1))
      for (i in idx) {
        if (times1[i] <= times1[j]) {
          second_term <- second_term + status1[i] / (risk_set[i] / n)^2
        }
      }
      second_term <- second_term / n

      xi[j] <- -surv_prob_at_last_event_time1 * (first_term - second_term)
    }
  }

  xi
}



#' @title Covariance matrix of cure-probability influence functions
#'
#' @description
#' Given lists of per-margin times and statuses, computes the covariance matrix of the influence functions.
#'
#' @param times_list List of numeric vectors (one per margin), all the same length (padded with NAs).
#' @param status_list List of integer (0/1) vectors, same structure as `times_list`.
#' @param surv_prob_at_last_event_time_list Numeric vector of KM estimates (length = number of margins).
#' @param time_points_list Numeric vector of evaluation times (length = number of margins).
#'
#' @return A \\(k \\times k\\) covariance matrix.
#'
#' @keywords internal
est.cov <- function(
    times_list,
    status_list,
    surv_prob_at_last_event_time_list,
    time_points_list
) {
  k     <- length(times_list)
  n_obs <- length(times_list[[1]])

  if (any(map_int(times_list, length) != n_obs) ||
      any(map_int(status_list, length) != n_obs)) {
    stop("All margins must have the same number of observations (use .long_to_lists to pad).")
  }

  influence_mat <- matrix(NA_real_, nrow = n_obs, ncol = k)
  for (i in seq_len(k)) {
    influence_mat[, i] <- est.inf.fun(
      times1  = times_list[[i]],
      status1 = status_list[[i]],
      surv_prob_at_last_event_time1 = surv_prob_at_last_event_time_list[[i]],
      time_point1 = time_points_list[[i]]
    )
  }

  (t(influence_mat) %*% influence_mat) / (n_obs^2)
}



#' @title Split and pad a long data frame by margin
#'
#' @description
#' Splits a long-format data frame by `margin`, extracts `time` and `status` into lists,
#' and pads each list with NAs so they are all equal length.
#'
#' @param df A data frame.
#' @param time Unquoted column name of follow-up time.
#' @param status Unquoted column name of event indicator (0/1).
#' @param margin Unquoted column name of grouping factor.
#'
#' @return A list with elements:
#' \describe{
#'   \item{times_list}{List of numeric vectors (padded).}
#'   \item{status_list}{List of integer vectors (padded).}
#' }
#'
#' @keywords internal
.long_to_lists <- function(df, time, status, margin) {
  time   <- enquo(time)
  status <- enquo(status)
  margin  <- enquo(margin)

  raw_split <- df %>%
    group_by(!!margin) %>%
    group_split()

  times_list  <- map(raw_split, ~ pull(.x, !!time))
  status_list <- map(raw_split, ~ pull(.x, !!status))

  max_n <- max(lengths(times_list))
  times_list  <- map(times_list,  ~ c(.x, rep(NA, max_n - length(.x))))
  status_list <- map(status_list, ~ c(.x, rep(NA, max_n - length(.x))))

  list(times_list = times_list, status_list = status_list)
}



#' @title Two-sample cure-rate test
#'
#' @description
#' Performs a Wald-type test comparing cure rates between two margins,
#' in either paired or independent fashion, and stores enough information
#' so that a 95 % CI for (p̂₁ − p̂₂) can be computed.
#'
#' @param df      A data frame in long format (one row per subject*margin).
#' @param time    Unquoted name of the follow-up time column.
#' @param status  Unquoted name of the event indicator column (0/1).
#' @param margin  Unquoted name of the grouping variable (two levels).
#' @param type    Character: `"paired"` (default) or `"indep"`.
#'
#' @return An object of class `"cure_test"`, which is a list containing:
#' \describe{
#'   \item{test_stat}{Numeric Wald statistic.}
#'   \item{p_value}{Two-sided p-value.}
#'   \item{type}{Either `"paired"` or `"indep"`.}
#'   \item{margin_names}{Character vector of the two margin levels.}
#'   \item{estimates}{Length-2 numeric vector of cure estimates (p̂₁, p̂₂).}
#'   \item{var_diff}{Variance of (p̂₁ − p̂₂).}
#'   \item{V}{2×2 covariance matrix of (p̂₁, p̂₂) if `type == "paired"`, otherwise NULL.}
#'   \item{v1}{Var(p̂₁) if `type == "indep"`, otherwise NULL.}
#'   \item{v2}{Var(p̂₂) if `type == "indep"`, otherwise NULL.}
#' }
#'
#' @export
test.2.cure <- function(df, time, status, margin, type = c("paired", "indep")) {
  type <- match.arg(type)

  # 1. Split & pad by margin (same as before)
  ls           <- .long_to_lists(df, time, status, margin)
  times_list   <- ls$times_list
  status_list  <- ls$status_list
  k            <- length(times_list)
  if (k != 2) stop("test.2.cure requires exactly two margins.")

  # 2. Compute per-margin KM cure estimate at the maximum follow-up (tau_hat):
  tau_hat <- numeric(2)
  p_hat   <- numeric(2)
  for (i in 1:2) {
    fit         <- survfit(Surv(times_list[[i]], status_list[[i]]) ~ 1)
    tau_hat[i]  <- max(times_list[[i]], na.rm = TRUE)
    p_hat[i]    <- summary(fit, times = tau_hat[i], extend = TRUE)$surv
  }

  # 3. Depending on paired vs. independent, compute variance of p̂₁ − p̂₂:
  if (type == "paired") {
    # 3a. For paired: get full 2×2 covariance matrix V
    V <- est.cov(
      times_list                        = times_list,
      status_list                       = status_list,
      surv_prob_at_last_event_time_list = p_hat,
      time_points_list                  = tau_hat
    )
    var_diff <- V[1,1] + V[2,2] - 2 * V[1,2]
    v1 <- v2 <- NULL
  } else {
    # 3b. For independent: compute each margin’s IF variance, v1 & v2
    xi1 <- est.inf.fun(times_list[[1]], status_list[[1]],
                       surv_prob_at_last_event_time1 = p_hat[1],
                       time_point1 = tau_hat[1])
    xi2 <- est.inf.fun(times_list[[2]], status_list[[2]],
                       surv_prob_at_last_event_time1 = p_hat[2],
                       time_point1 = tau_hat[2])
    v1 <- sum(xi1^2) / (length(xi1)^2)
    v2 <- sum(xi2^2) / (length(xi2)^2)
    var_diff <- v1 + v2
    V <- NULL
  }

  # 4. Wald statistic and p-value:
  stat <- (p_hat[1] - p_hat[2]) / sqrt(var_diff)
  pval <- 2 * pnorm(-abs(stat))

  # 5. Capture margin (group) names for printing:
  margin_var   <- deparse(substitute(margin))
  margin_names <- levels(factor(df[[margin_var]]))

  # 6. Return an S3 object of class "cure_test"
  structure(
    list(
      test_stat    = stat,
      p_value      = pval,
      type         = type,
      margin_names = margin_names,
      estimates    = p_hat,
      var_diff     = var_diff,
      V            = V,
      v1           = v1,
      v2           = v2
    ),
    class = "cure_test"
  )
}




#' @title Wald test across K margins for cure probabilities
#'
#' @description
#' Tests equality of cure probabilities across K margins using a multivariate Wald statistic.
#'
#' @param df      A data frame in long format.
#' @param time    Unquoted name of the time column.
#' @param status  Unquoted name of the event indicator column (0/1).
#' @param margin  Unquoted name of the grouping factor (>= 2 levels).
#'
#' @return An object of class `"cure_k_test"`, a list containing:
#' \describe{
#'   \item{test_stat}{Wald statistic (χ²).}
#'   \item{p_value}{p-value of χ²ₖ₋₁ test.}
#'   \item{df}{Degrees of freedom (k−1).}
#'   \item{invertible}{Logical, was the contrast covariance invertible?}
#'   \item{margin_names}{Character vector of all margin levels.}
#' }
#'
#' @export
test.k.cure <- function(df, time, status, margin) {
  # 1. Split & pad
  ls           <- .long_to_lists(df, time, status, margin)
  times_list   <- ls$times_list
  status_list  <- ls$status_list
  k            <- length(times_list)
  if (k < 2) stop("test.k.cure requires at least two margins.")

  # 2. Compute per-margin KM cure estimates:
  tau_hat <- p_hat <- numeric(k)
  for (i in seq_len(k)) {
    fit         <- survfit(Surv(times_list[[i]], status_list[[i]]) ~ 1)
    tau_hat[i]  <- max(times_list[[i]], na.rm = TRUE)
    p_hat[i]    <- summary(fit, times = tau_hat[i], extend = TRUE)$surv
  }

  # 3. Covariance matrix of all p̂’s:
  V <- est.cov(
    times_list                        = times_list,
    status_list                       = status_list,
    surv_prob_at_last_event_time_list = p_hat,
    time_points_list                  = tau_hat
  )

  # 4. Contrast matrix C and solve C V Cᵀ:
  C   <- construct.C(k)          # (k−1)×k adjacent-differences
  sol <- solve.system(C %*% V %*% t(C), C %*% p_hat)
  W   <- as.numeric(t(C %*% p_hat) %*% sol$solve)  # Wald statistic
  pval <- 1 - pchisq(W, df = k - 1)

  # 5. Capture margin (group) names:
  margin_var   <- deparse(substitute(margin))
  margin_names <- levels(factor(df[[margin_var]]))

  # 6. Return an S3 object of class "cure_k_test"
  structure(
    list(
      test_stat    = W,
      p_value      = pval,
      df           = k - 1,
      invertible   = sol$invertible,
      margin_names = margin_names
    ),
    class = "cure_k_test"
  )
}



#' @title Estimate cure probabilities (weighted, average, pooled, Ying-Wei)
#'
#' @description
#' Computes three overall estimators of the cure rate—covariance-weighted,
#' simple average, and pooled Kaplan–Meier—as well as the Ying-Wei variance-adjusted KM.
#'
#' @param df      A data frame in long format.
#' @param time    Unquoted name of the time column.
#' @param status  Unquoted name of the event indicator column (0/1).
#' @param margin  Unquoted name of the grouping factor (≥ 2 levels).
#'
#' @return An object of class `"cure_est"`, which is a list containing:
#' \describe{
#'   \item{weight_vec}{Weights for the covariance-weighted estimator.}
#'   \item{weighted_est}{Weighted cure-rate estimate.}
#'   \item{var_weighted_est}{Variance of the weighted estimate.}
#'   \item{ave_est}{Simple-average cure-rate estimate.}
#'   \item{var_ave_est}{Variance of the simple average.}
#'   \item{pool_est}{Pooled Kaplan–Meier cure probability.}
#'   \item{var_pool_est}{Variance of the pooled estimate.}
#'   \item{yw_est}{Ying-Wei variance-corrected KM estimate.}
#'   \item{var_yw_est}{Variance of the Ying-Wei estimate.}
#'   \item{cure_prob_list}{Per-margin KM cure probabilities (length = K).}
#'   \item{var_cure_prob}{Per-margin variances (diagonal of V).}
#'   \item{invertible}{Logical—was V invertible?}
#'   \item{margin_names}{Character vector of margin levels.}
#' }
#'
#' @export
est.cure <- function(df, time, status, margin) {
  # 1. Split & pad by margin
  ls           <- .long_to_lists(df, time, status, margin)
  times_list   <- ls$times_list
  status_list  <- ls$status_list
  k            <- length(times_list)
  if (k < 1) stop("est.cure: need at least one margin.")

  # 2. Compute per-margin KM cure estimates at tau_hat[i]:
  tau_hat <- p_hat <- numeric(k)
  for (i in seq_len(k)) {
    fit         <- survfit(Surv(times_list[[i]], status_list[[i]]) ~ 1)
    tau_hat[i]  <- max(times_list[[i]], na.rm = TRUE)
    p_hat[i]    <- summary(fit, times = tau_hat[i], extend = TRUE)$surv
  }

  # 3. Covariance matrix V of all p̂’s:
  V <- est.cov(
    times_list                        = times_list,
    status_list                       = status_list,
    surv_prob_at_last_event_time_list = p_hat,
    time_points_list                  = tau_hat
  )

  # 4. Covariance-weighted estimator:
  one_vec          <- rep(1, k)
  sol1             <- solve.system(V, one_vec)
  inv_one          <- sol1$solve
  weight_vec       <- inv_one / sum(inv_one)     # weights
  weighted_est     <- sum(weight_vec * p_hat)     # weighted cure estimate
  var_weighted_est <- 1 / sum(inv_one)           # Var(weighted)

  # 5. Simple average estimator:
  ave_est     <- mean(p_hat)
  var_ave_est <- as.numeric((one_vec / k) %*% V %*% (one_vec / k))

  # 6. Pooled Kaplan–Meier:
  pooled_times  <- unlist(times_list)
  pooled_status <- unlist(status_list)
  fit_pool      <- survfit(Surv(pooled_times, pooled_status) ~ 1)
  pool_est      <- summary(fit_pool, times = max(tau_hat), extend = TRUE)$surv
  var_pool_est  <- (summary(fit_pool, times = max(tau_hat), extend = TRUE)$std.err)^2

  # 7. Ying-Wei (clustered) estimator:
  n_all    <- length(pooled_times)
  pooled_df <- data.frame(
    time    = pooled_times,
    status  = pooled_status,
    subject = seq_len(n_all)
  )
  km_yw     <- prodlim::prodlim(Hist(time, status) ~ cluster(subject),
                                data = pooled_df)
  yw_sum    <- summary(km_yw, times = max(tau_hat), extend = TRUE)
  yw_est    <- yw_sum$surv
  var_yw_est <- (yw_sum$se.surv)^2

  # 8. Capture margin (group) names:
  margin_var   <- deparse(substitute(margin))
  margin_names <- levels(factor(df[[margin_var]]))

  # 9. Return an S3 object of class "cure_est"
  structure(
    list(
      weight_vec        = weight_vec,
      weighted_est      = weighted_est,
      var_weighted_est  = var_weighted_est,
      ave_est           = ave_est,
      var_ave_est       = var_ave_est,
      pool_est          = pool_est,
      var_pool_est      = var_pool_est,
      yw_est            = yw_est,
      var_yw_est        = var_yw_est,
      cure_prob_list    = p_hat,
      var_cure_prob     = diag(V),
      invertible        = sol1$invertible,
      margin_names      = margin_names
    ),
    class = "cure_est"
  )
}

