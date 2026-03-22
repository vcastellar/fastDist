test_that("fdist correlation matches base cor-derived distance", {
  A <- matrix(c(
    1, 2, 3,
    2, 5, 6,
    3, 1, 0
  ), nrow = 3, byrow = TRUE)
  B <- matrix(c(
    2, 1, 4,
    4, 4, 1
  ), nrow = 2, byrow = TRUE)

  expected <- 1 - cor(t(A), t(B))

  expect_equal(
    fdist(A, B, method = "correlation"),
    expected,
    tolerance = 1e-8
  )
})

test_that("fdist correlation returns a symmetric zero-diagonal matrix for identical input", {
  A <- matrix(c(
    1, 2, 3,
    2, 5, 6,
    3, 1, 0
  ), nrow = 3, byrow = TRUE)

  result <- fdist(A, method = "correlation")

  expect_equal(result, t(result), tolerance = 1e-8)
  expect_equal(diag(result), rep(0, nrow(A)), tolerance = 1e-8)
})
