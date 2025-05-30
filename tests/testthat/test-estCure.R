# tests/testthat/test-estCure.R

test_that("est.cure returns all expected components and valid estimates", {
  df <- data.frame(
    id     = rep(1:4, each = 2),
    group  = factor(rep(c("G1","G2"), each = 4)),
    time   = c(1,2,3,4, 1,2,3,4),
    status = c(1,1,0,0, 1,0,1,0)
  )
  
  res <- est.cure(df, time, status, group)
  
  expected <- c(
    "weight_vec","weighted_est","var_weighted_est",
    "ave_est","var_ave_est","pool_est","var_pool_est",
    "yw_est","var_yw_est","cure_prob_list",
    "var_cure_prob","invertible"
  )
  expect_named(res, expected)
  
  # All estimates should be within [0,1]
  expect_true(res$weighted_est >= 0 && res$weighted_est <= 1)
  expect_true(res$ave_est      >= 0 && res$ave_est      <= 1)
  expect_true(res$pool_est     >= 0 && res$pool_est     <= 1)
})
