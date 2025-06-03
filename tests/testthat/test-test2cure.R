# tests/testthat/test-test2cure.R

test_that("paired two‐sample cure‐rate test returns an object of class 'cure_test' with correct components", {
  df <- data.frame(
    id     = rep(1:2, each = 2),
    margin = factor(rep(c("A", "B"), times = 2)),
    time   = c(5, 10, 5, 10),
    status = c(1, 0, 1, 0)
  )

  res <- test.2.cure(df, time, status, margin, type = "paired")

  # 1) Should be an S3 object of class "cure_test"
  expect_s3_class(res, "cure_test")

  # 2) It should contain at least these named components:
  expected_fields <- c(
    "test_stat", "p_value", "type", "margin_names",
    "estimates", "var_diff", "V", "v1", "v2"
  )
  expect_true(all(expected_fields %in% names(res)))

  # 3) test_stat and p_value should be numeric, p_value in [0,1]
  expect_true(is.numeric(res$test_stat))
  expect_true(is.numeric(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)

  # 4) type should equal "paired", and margin_names should be exactly c("A","B")
  expect_equal(res$type, "paired")
  expect_equal(res$margin_names, c("A", "B"))

  # 5) estimates should be a length-2 numeric vector in [0,1]
  expect_true(is.numeric(res$estimates) && length(res$estimates) == 2)
  expect_true(all(res$estimates >= 0 & res$estimates <= 1))
})

test_that("paired test errors when there are not exactly two margins", {
  df3 <- data.frame(
    id     = 1:3,
    margin = factor(c("A", "B", "C")),
    time   = c(1, 2, 3),
    status = c(1, 0, 1)
  )

  expect_error(
    test.2.cure(df3, time, status, margin, type = "paired"),
    "requires exactly two margins"
  )
})
