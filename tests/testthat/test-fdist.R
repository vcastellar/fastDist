# Tests for fdist() and its distance backends.
#
# Each metric is checked against an independent, pure-R reference
# implementation. Where the optional packages `proxy` and `parallelDist`
# are installed, a few cross-checks are added as well.

tol <- 1e-8

# Pure-R reference: applies `fun(x, y)` to every pair of rows (x from A, y
# from B) and returns the resulting m x n matrix.
ref_dist <- function(A, B, fun) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  out <- matrix(0, nrow(A), nrow(B))
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(nrow(B))) {
      out[i, j] <- fun(A[i, ], B[j, ])
    }
  }
  out
}

# Per-metric scalar reference functions.
ref_funs <- list(
  euclidean   = function(x, y) sqrt(sum((x - y)^2)),
  sqeuclidean = function(x, y) sum((x - y)^2),
  manhattan   = function(x, y) sum(abs(x - y)),
  supremum    = function(x, y) max(abs(x - y)),
  canberra    = function(x, y) {
    den <- abs(x) + abs(y)
    sum(ifelse(den > 0, abs(x - y) / den, 0))
  },
  correlation = function(x, y) 1 - stats::cor(x, y),
  cosine      = function(x, y) 1 - sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))),
  hamming     = function(x, y) sum(x != y)
)

ref_jaccard <- function(x, y) {
  a  <- sum(x == 1 & y == 1)
  bc <- sum((x == 1) != (y == 1))
  if (a + bc > 0) bc / (a + bc) else 0
}

ref_dice <- function(x, y) {
  a  <- sum(x == 1 & y == 1)
  bc <- sum((x == 1) != (y == 1))
  den <- 2 * a + bc
  if (den > 0) bc / den else 0
}

# ---------------------------------------------------------------------------
# Continuous-valued metrics, A vs B
# ---------------------------------------------------------------------------
test_that("continuous metrics match the pure-R reference (A vs B)", {
  set.seed(1)
  A <- matrix(rnorm(40), 8, 5)
  B <- matrix(rnorm(30), 6, 5)

  for (m in names(ref_funs)) {
    if (m == "hamming") next  # exercised separately
    expect_equal(
      fdist(A, B, method = m),
      ref_dist(A, B, ref_funs[[m]]),
      tolerance = tol,
      info = m
    )
  }
})

test_that("minkowski matches the reference for several values of p", {
  set.seed(2)
  A <- matrix(rnorm(40), 8, 5)
  B <- matrix(rnorm(30), 6, 5)

  for (p in c(1, 2, 3, 5)) {
    ref <- ref_dist(A, B, function(x, y) sum(abs(x - y)^p)^(1 / p))
    expect_equal(fdist(A, B, method = "minkowski", p = p), ref,
                 tolerance = tol, info = paste("p =", p))
  }
})

test_that("minkowski with p = 1/2 reproduces manhattan/euclidean", {
  set.seed(3)
  A <- matrix(rnorm(40), 8, 5)
  B <- matrix(rnorm(30), 6, 5)

  expect_equal(fdist(A, B, method = "minkowski", p = 1),
               fdist(A, B, method = "manhattan"), tolerance = tol)
  expect_equal(fdist(A, B, method = "minkowski", p = 2),
               fdist(A, B, method = "euclidean"), tolerance = tol)
})

test_that("sqeuclidean equals euclidean squared", {
  set.seed(4)
  A <- matrix(rnorm(40), 8, 5)
  B <- matrix(rnorm(30), 6, 5)

  expect_equal(fdist(A, B, method = "sqeuclidean"),
               fdist(A, B, method = "euclidean")^2, tolerance = tol)
})

# ---------------------------------------------------------------------------
# Symmetric case (B = NULL)
# ---------------------------------------------------------------------------
test_that("B = NULL reproduces the A vs A result with zero diagonal", {
  set.seed(5)
  A <- matrix(rnorm(40), 8, 5)

  for (m in c("euclidean", "sqeuclidean", "manhattan", "supremum",
              "canberra", "correlation", "cosine")) {
    res <- fdist(A, method = m)
    expect_equal(res, fdist(A, A, method = m), tolerance = tol, info = m)
    expect_equal(res, t(res), tolerance = tol, info = m)        # symmetric
    expect_equal(diag(res), rep(0, nrow(A)), tolerance = tol, info = m)
  }
})

# ---------------------------------------------------------------------------
# Hamming
# ---------------------------------------------------------------------------
test_that("hamming counts differing coordinates", {
  set.seed(6)
  A <- matrix(sample(0:3, 40, replace = TRUE), 8, 5)
  B <- matrix(sample(0:3, 30, replace = TRUE), 6, 5)

  expect_equal(fdist(A, B, method = "hamming"),
               ref_dist(A, B, ref_funs$hamming), tolerance = tol)
})

# ---------------------------------------------------------------------------
# Binary metrics
# ---------------------------------------------------------------------------
test_that("jaccard and dice match the reference on binary data", {
  set.seed(7)
  A <- matrix(sample(0:1, 40, replace = TRUE), 8, 5)
  B <- matrix(sample(0:1, 30, replace = TRUE), 6, 5)

  expect_equal(fdist(A, B, method = "jaccard"),
               ref_dist(A, B, ref_jaccard), tolerance = tol)
  expect_equal(fdist(A, B, method = "dice"),
               ref_dist(A, B, ref_dice), tolerance = tol)
})

test_that("jaccard/dice treat all-zero pairs as distance 0", {
  A <- matrix(0, 3, 4)
  expect_equal(fdist(A, method = "jaccard"), matrix(0, 3, 3))
  expect_equal(fdist(A, method = "dice"), matrix(0, 3, 3))
})

test_that("binary methods reject non-binary input", {
  bad <- matrix(c(0, 1, 2, 1), 2, 2)
  ok  <- matrix(c(0, 1, 1, 0), 2, 2)

  expect_error(fdist(bad, method = "jaccard"), "binary")
  expect_error(fdist(bad, method = "dice"), "binary")
  expect_error(fdist(ok, bad, method = "jaccard"), "binary")

  withNA <- matrix(c(0, 1, NA, 1), 2, 2)
  expect_error(fdist(withNA, method = "jaccard"), "NA")
})

# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------
test_that("unknown method raises an error", {
  A <- matrix(rnorm(20), 4, 5)
  expect_error(fdist(A, method = "not-a-method"), "not found")
})

test_that("data frames are accepted via as.matrix coercion", {
  df <- as.data.frame(matrix(rnorm(20), 4, 5))
  expect_equal(fdist(df, method = "euclidean"),
               fdist(as.matrix(df), method = "euclidean"), tolerance = tol)
})

# ---------------------------------------------------------------------------
# Cross-checks against proxy (optional dependency)
# ---------------------------------------------------------------------------
test_that("classic metrics agree with proxy::dist", {
  skip_if_not_installed("proxy")
  set.seed(8)
  A <- matrix(rnorm(40), 8, 5)
  B <- matrix(rnorm(30), 6, 5)

  mapping <- c(euclidean = "Euclidean", manhattan = "Manhattan",
               supremum = "Supremum", canberra = "Canberra",
               correlation = "correlation")

  for (m in names(mapping)) {
    expected <- as.matrix(proxy::dist(A, B, method = mapping[[m]]))
    expect_equal(unname(fdist(A, B, method = m)), unname(expected),
                 tolerance = 1e-6, info = m)
  }

  expected <- as.matrix(proxy::dist(A, B, method = "Minkowski", p = 3))
  expect_equal(unname(fdist(A, B, method = "minkowski", p = 3)),
               unname(expected), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Cross-check against parallelDist (optional dependency)
# ---------------------------------------------------------------------------
test_that("symmetric euclidean/manhattan agree with parallelDist", {
  skip_if_not_installed("parallelDist")
  set.seed(9)
  A <- matrix(rnorm(50), 10, 5)

  for (m in c("euclidean", "manhattan")) {
    expected <- as.matrix(parallelDist::parDist(A, method = m))
    expect_equal(unname(fdist(A, method = m)), unname(expected),
                 tolerance = 1e-6, info = m)
  }
})
