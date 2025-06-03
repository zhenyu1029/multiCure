# tests/testthat/test-estCure.R

test_that("est.cure returns an object of class 'cure_est' with all expected components", {
  df <- data.frame(
    id      = rep(1:4, each = 2),
    margin  = factor(rep(c("G1", "G2"), each = 4)),
    time    = c(1, 2, 3, 4, 1, 2, 3, 4),
    status  = c(1, 1, 0, 0, 1, 0, 1, 0)
  )

  res <- est.cure(df, time, status, margin)

  # 1) It should be an S3 object of class "cure_est"
  expect_s3_class(res, "cure_est")

  # 2) It should contain exactly these named components (in this order):
  expected_names <- c(
    "weight_vec", "weighted_est", "var_weighted_est",
    "ave_est",     "var_ave_est",    "pool_est",     "var_pool_est",
    "yw_est",      "var_yw_est",     "cure_prob_list",
    "var_cure_prob", "invertible",  "margin_names"
  )
  expect_named(res, expected_names)

  # 3) All global estimates should lie in [0, 1]
  expect_true(res$weighted_est >= 0 && res$weighted_est <= 1)
  expect_true(res$ave_est      >= 0 && res$ave_est      <= 1)
  expect_true(res$pool_est     >= 0 && res$pool_est     <= 1)
  expect_true(res$yw_est       >= 0 && res$yw_est       <= 1)

  # 4) Per-margin cure probabilities ("cure_prob_list") also lie in [0,1]
  expect_true(all(res$cure_prob_list >= 0 & res$cure_prob_list <= 1))

  # 5) margin_names should match the factor levels in df
  expect_equal(res$margin_names, levels(df$margin))
})
