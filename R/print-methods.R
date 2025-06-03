################################################################################
# R/print-methods.R
################################################################################

#---------------------------------------------
# 1) print.cure_test: simple print for test.2.cure()
#---------------------------------------------
#' @export
print.cure_test <- function(x, ...) {
  # x is an object of class "cure_test"
  header <- switch(
    x$type,
    paired = "Two-Sample Cure Rate Test (Paired)",
    indep  = "Two-Sample Cure Rate Test (Independent)",
    "Two-Sample Cure Rate Test"
  )
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n", sep = "")

  # Show which margins were compared
  if (!is.null(x$margin_names) && length(x$margin_names) == 2) {
    cat(sprintf("data:  %s vs %s\n",
                x$margin_names[1],
                x$margin_names[2]))
  }

  # Test statistic & p-value
  cat(sprintf("Test Statistic:  %.4f\n", x$test_stat))
  cat(sprintf("P-Value:         %.4g\n", x$p_value))

  # Basic alternative statement
  cat("Alternative Hypothesis: cure probabilities differ\n")

  invisible(x)
}

#---------------------------------------------
# 2a) summary.cure_test: compute 95% CI for p1-p2
#---------------------------------------------
#' @export
summary.cure_test <- function(object, ...) {
  if (!inherits(object, "cure_test")) {
    stop("`summary.cure_test()` only applies to objects of class 'cure_test'.")
  }
  p_hat    <- object$estimates      # length-2 vector
  diff_hat <- p_hat[1] - p_hat[2]
  var_diff <- object$var_diff
  se_diff  <- sqrt(var_diff)

  # 95% CI:
  alpha   <- 0.05
  z       <- qnorm(1 - alpha/2)
  ci_low  <- diff_hat - z * se_diff
  ci_high <- diff_hat + z * se_diff

  out <- list(
    margin_names = object$margin_names,
    p_hat        = p_hat,
    diff_hat     = diff_hat,
    se_diff      = se_diff,
    conf_int     = c(lower = ci_low, upper = ci_high, level = 0.95),
    test_stat    = object$test_stat,
    p_value      = object$p_value,
    type         = object$type
  )
  class(out) <- c("summary.cure_test", "list")
  out
}

#---------------------------------------------
# 2b) print.summary.cure_test: nice display of estimates + CI
#---------------------------------------------
#' @export
print.summary.cure_test <- function(x, ...) {
  # x is a "summary.cure_test"
  margin1 <- x$margin_names[1]
  margin2 <- x$margin_names[2]

  # Header
  header <- switch(
    x$type,
    paired = "Two-Sample Cure Rate Test (Paired, 95% CI)",
    indep  = "Two-Sample Cure Rate Test (Independent, 95% CI)",
    "Two-Sample Cure Rate Test (95% CI)"
  )
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n", sep = "")

  # Show each margin's estimated cure rate
  cat(sprintf("Estimated cure rate:\n  %s: %.4f\n  %s: %.4f\n\n",
              margin1, x$p_hat[1],
              margin2, x$p_hat[2]))

  # Show difference + CI
  cat(sprintf("Difference ( %s - %s ) = %.4f\n",
              margin1, margin2, x$diff_hat))
  cat(sprintf("95%% CI for difference: [ %.4f , %.4f ]\n\n",
              x$conf_int["lower"], x$conf_int["upper"]))

  # Show Wald test result
  cat(sprintf("Test statistic = %.4f   p-value = %.4g\n",
              x$test_stat, x$p_value))

  # Basic alternative statement
  cat("Alternative Hypothesis: cure probabilities differ\n")

  invisible(x)

}

#---------------------------------------------
# 3) print.cure_k_test: print for test.k.cure()
#---------------------------------------------
#' @export
print.cure_k_test <- function(x, ...) {
  # x is an object of class "cure_k_test"
  cat("K-Sample Cure Rate Wald Test\n")
  cat("----------------------------\n")

  # List the margins
  if (!is.null(x$margin_names)) {
    cat("Margins: ", paste(x$margin_names, collapse = ", "), "\n")
  }

  # Test statistic, df, p-value
  cat(sprintf("Test Statistic (Chi^2):  %.4f\n", x$test_stat))
  cat(sprintf("Degrees of Freedom:      %d\n", x$df))
  cat(sprintf("P-Value:                 %.4g\n", x$p_value))

  # State the alternative hypothesis
  cat("Alternative Hypothesis: All margins have different cure probabilities\n")

  # Warn if covariance was singular
  if (!x$invertible) {
    cat("WARNING: Covariance matrix was not invertible\n")
  }
  invisible(x)
}

#---------------------------------------------
# 4) summary.cure_k_test: compute and print details for test.k.cure()
#---------------------------------------------
#' @export
summary.cure_k_test <- function(object, ...) {
  if (!inherits(object, "cure_k_test")) {
    stop("`summary.cure_k_test()` only applies to objects of class 'cure_k_test'.")
  }
  out <- list(
    margin_names = object$margin_names,
    test_stat    = object$test_stat,
    df           = object$df,
    p_value      = object$p_value,
    invertible   = object$invertible
  )
  class(out) <- c("summary.cure_k_test", "list")
  out
}

#---------------------------------------------
# 5) print.summary.cure_k_test: nice display of K-sample test
#---------------------------------------------
#' @export
print.summary.cure_k_test <- function(x, ...) {
  # x is a "summary.cure_k_test"
  header <- "K-Sample Cure Rate Wald Test (Summary)"
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n", sep = "")

  # Show margin names
  cat("Margins: ", paste(x$margin_names, collapse = ", "), "\n\n")

  # Print test statistic, df, p-value
  cat(sprintf("Test Statistic (Chi^2):  %.4f\n", x$test_stat))
  cat(sprintf("Degrees of Freedom:      %d\n", x$df))
  cat(sprintf("P-Value:                 %.4g\n\n", x$p_value))

  # State the alternative hypothesis
  cat("Alternative Hypothesis: At least one margin's cure probability differs.\n\n")

  # Warn if covariance was singular
  if (!x$invertible) {
    cat("WARNING: Covariance matrix was not invertible\n")
  }
  invisible(x)
}


#---------------------------------------------
# 6) print.cure_est: print for est.cure()
#---------------------------------------------
#' @export
print.cure_est <- function(x, ...) {
  # x is an object of class "cure_est"
  header <- "Cure-Rate Estimates"
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n", sep = "")

  # 4a) Weighted estimator
  cat(sprintf("Weighted estimator:     %.4f (Var = %.4f)\n",
              x$weighted_est, x$var_weighted_est))

  # 4b) Simple average estimator
  cat(sprintf("Simple average:        %.4f (Var = %.4f)\n",
              x$ave_est, x$var_ave_est))

  # 4c) Pooled KM estimator
  cat(sprintf("Pooled KM estimator:   %.4f (Var = %.4f)\n",
              x$pool_est, x$var_pool_est))

  # 4d) Ying-Wei estimator
  cat(sprintf("Ying-Wei estimator:    %.4f (Var = %.4f)\n",
              x$yw_est, x$var_yw_est))

  # 4e) Per-margin cure rates
  if (!is.null(x$margin_names) && length(x$margin_names) == length(x$cure_prob_list)) {
    cat("\nPer-margin cure rates:\n")
    for (i in seq_along(x$margin_names)) {
      cat(sprintf("  %s: %.4f  (Var = %.4f)\n",
                  x$margin_names[i],
                  x$cure_prob_list[i],
                  x$var_cure_prob[i]))
    }
  }

  # 4f) Warn if weighted covariance was singular
  if (!x$invertible) {
    cat("\nWARNING: Covariance matrix was not invertible-weighted estimator may be unreliable.\n")
  }
  invisible(x)
}

################################################################################
# End of R/print-methods.R
################################################################################
