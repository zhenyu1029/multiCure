# tests/testthat/test-test2cure.R

test_that("paired two-sample cure-rate test returns correct names and p-value bounds", {
  df <- data.frame(
    id     = rep(1:2, each = 2),
    margin  = factor(rep(c("A","B"), times = 2)),
    time   = c(5,10, 5,10),
    status = c(1,0, 1,0)
  )

  res <- test.2.cure(df, time, status, margin, type = "paired")

  expect_named(res, c("test_stat","p_value"))
  expect_true(is.numeric(res$test_stat))
  expect_true(is.numeric(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("paired test errors when there are not exactly two margins", {
  df3 <- data.frame(
    id     = 1:3,
    margin  = factor(c("A","B","C")),
    time   = c(1,2,3),
    status = c(1,0,1)
  )

  expect_error(
    test.2.cure(df3, time, status, margin, type = "paired"),
    "requires exactly two margins"
  )
})
