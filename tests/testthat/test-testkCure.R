# tests/testthat/test-testkCure.R

test_that("test.k.cure returns a 'cure_k_test' object and handles invertibility", {
  df <- data.frame(
    id     = rep(1:3, each = 3),
    margin = factor(rep(c("G1", "G2", "G3"), each = 3)),
    time   = rep(1:3, times = 3),
    status = c(
      1, 0, 1,  # G1
      1, 1, 0,  # G2
      0, 1, 1   # G3
    )
  )

  res <- test.k.cure(df, time, status, margin)

  # 1) Should be an S3 object of class "cure_k_test"
  expect_s3_class(res, "cure_k_test")

  # 2) It should contain at least these components:
  expected_fields <- c("test_stat", "p_value", "df", "invertible", "margin_names")
  expect_true(all(expected_fields %in% names(res)))

  # 3) Degrees of freedom should be k-1 = 2
  expect_equal(res$df, 2)

  # 4) invertible should be logical
  expect_type(res$invertible, "logical")

  # 5) margin_names should be exactly the factor levels
  expect_equal(res$margin_names, levels(df$margin))

  # 6) If invertible, p_value must be numeric in [0,1]; otherwise p_value is NA
  if (res$invertible) {
    expect_true(is.numeric(res$p_value))
    expect_true(res$p_value >= 0 && res$p_value <= 1)
  } else {
    expect_true(is.na(res$p_value))
  }
})
