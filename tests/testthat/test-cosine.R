test_that("fdist cosine matches manual cosine distance", {
  A <- matrix(c(
    1, 0, 0,
    1, 1, 0,
    0, 1, 1
  ), nrow = 3, byrow = TRUE)
  B <- matrix(c(
    1, 0, 1,
    0, 1, 0
  ), nrow = 2, byrow = TRUE)

  expected <- matrix(0, nrow(A), nrow(B))
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(nrow(B))) {
      denom <- sqrt(sum(A[i, ]^2)) * sqrt(sum(B[j, ]^2))
      cos_sim <- if (denom > 0) sum(A[i, ] * B[j, ]) / denom else 0
      expected[i, j] <- 1 - cos_sim
    }
  }

  expect_equal(
    fdist(A, B, method = "cosine"),
    expected,
    tolerance = 1e-8
  )
})

test_that("fdist cosine returns a symmetric zero-diagonal matrix for identical input", {
  A <- matrix(c(
    1, 2, 3,
    2, 5, 6,
    3, 1, 0
  ), nrow = 3, byrow = TRUE)

  result <- fdist(A, method = "cosine")

  expect_equal(result, t(result), tolerance = 1e-8)
  expect_equal(diag(result), rep(0, nrow(A)), tolerance = 1e-8)
})
