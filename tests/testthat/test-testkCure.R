# tests/testthat/test-testkCure.R

test_that("test.k.cure works for 3 margins and handles non-invertibility", {
  df <- data.frame(
    id     = rep(1:3, each = 3),
    margin  = factor(rep(c("G1","G2","G3"), each = 3)),
    time   = rep(1:3, times = 3),
    status = c(
      1,0,1,  # G1
      1,1,0,  # G2
      0,1,1   # G3
    )
  )

  res <- test.k.cure(df, time, status, margin)

  expect_named(res, c("test_stat","p_value","df","invertible"))
  expect_equal(res$df, 2)
  expect_type(res$invertible, "logical")

  if (res$invertible) {
    expect_true(is.numeric(res$p_value) && res$p_value >= 0 && res$p_value <= 1)
  } else {
    expect_true(is.na(res$p_value))
  }
})
