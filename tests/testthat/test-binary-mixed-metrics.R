# Validation of the Hamming, Jaccard and Gower backends against base-R
# references.

pairwise <- function(A, B, f) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  m <- nrow(A)
  n <- nrow(B)
  R <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      R[i, j] <- f(A[i, ], B[j, ])
    }
  }
  R
}

set.seed(7)
# binary data
Xb <- matrix(rbinom(6 * 10, 1, 0.4), 6, 10)
Yb <- matrix(rbinom(4 * 10, 1, 0.4), 4, 10)
# general numeric data
A <- matrix(rnorm(5 * 4), 5, 4)
B <- matrix(rnorm(3 * 4), 3, 4)

tol <- 1e-8

test_that("hamming counts differing positions", {
  ref <- pairwise(Xb, Yb, function(x, y) sum(x != y))
  expect_equal(fdist(Xb, Yb, method = "hamming"), ref, tolerance = tol)
  # also works on categorical (integer-coded) data
  Xc <- matrix(sample(1:3, 5 * 6, replace = TRUE), 5, 6)
  Yc <- matrix(sample(1:3, 4 * 6, replace = TRUE), 4, 6)
  refc <- pairwise(Xc, Yc, function(x, y) sum(x != y))
  expect_equal(fdist(Xc, Yc, method = "hamming"), refc, tolerance = tol)
})

test_that("jaccard matches 1 - intersection/union on non-zero patterns", {
  jac <- function(x, y) {
    u <- sum(x != 0 | y != 0)
    if (u == 0) 0 else 1 - sum(x != 0 & y != 0) / u
  }
  ref <- pairwise(Xb, Yb, jac)
  expect_equal(fdist(Xb, Yb, method = "jaccard"), ref, tolerance = tol)
  # two all-zero rows are at distance zero
  Z <- rbind(rep(0, 4), rep(0, 4))
  expect_equal(fdist(Z, method = "jaccard"), matrix(0, 2, 2), tolerance = tol)
})

test_that("gower matches range-scaled mean absolute difference", {
  rng <- apply(A, 2, function(col) {
    r <- max(col) - min(col)
    if (r == 0) 1 else r
  })
  gw <- function(x, y) mean(abs(x - y) / rng)
  ref <- pairwise(A, B, gw)
  expect_equal(fdist(A, B, method = "gower"), ref, tolerance = tol)
})

test_that("gower is invariant to per-feature linear rescaling of both inputs", {
  scale_v <- c(1, 100, 0.01, 10)
  A2 <- sweep(A, 2, scale_v, "*")
  B2 <- sweep(B, 2, scale_v, "*")
  expect_equal(fdist(A2, B2, method = "gower"),
               fdist(A, B, method = "gower"), tolerance = tol)
})

test_that("new metrics return symmetric matrices with zero diagonal", {
  for (m in c("hamming", "jaccard", "gower")) {
    M <- if (m == "gower") A else Xb
    d <- fdist(M, method = m)
    expect_equal(d, t(d), tolerance = tol, info = m)
    expect_equal(diag(d), rep(0, nrow(M)), tolerance = tol, info = m)
  }
})

test_that("minkowski defaults to p = 2 (euclidean)", {
  expect_equal(fdist(A, B, method = "minkowski"),
               fdist(A, B, method = "euclidean"), tolerance = tol)
})

test_that("unknown methods raise an error", {
  expect_error(fdist(A, B, method = "no_such_method"))
})
