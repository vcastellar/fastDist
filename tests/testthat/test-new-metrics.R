# Validation of the Tier 1 distance backends against base-R references.

# Generic pairwise distance using a per-pair function f(x, y).
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

set.seed(42)
# general numeric data (can be negative)
A <- matrix(rnorm(5 * 4), 5, 4)
B <- matrix(rnorm(3 * 4), 3, 4)

# strictly positive data for distribution / abundance metrics
Ap <- matrix(runif(5 * 4, 0.1, 5), 5, 4)
Bp <- matrix(runif(3 * 4, 0.1, 5), 3, 4)

tol <- 1e-8

test_that("squared_euclidean matches sum of squared differences", {
  ref <- pairwise(A, B, function(x, y) sum((x - y)^2))
  expect_equal(fdist(A, B, method = "squared_euclidean"), ref, tolerance = tol)
  # equals the square of the euclidean distance
  expect_equal(fdist(A, B, method = "squared_euclidean"),
               fdist(A, B, method = "euclidean")^2, tolerance = tol)
})

test_that("bray_curtis matches reference", {
  bc <- function(x, y) {
    den <- sum(abs(x + y))
    if (den == 0) 0 else sum(abs(x - y)) / den
  }
  ref <- pairwise(Ap, Bp, bc)
  expect_equal(fdist(Ap, Bp, method = "bray_curtis"), ref, tolerance = tol)
})

test_that("hellinger matches reference", {
  hl <- function(x, y) sqrt(sum((sqrt(x) - sqrt(y))^2)) / sqrt(2)
  ref <- pairwise(Ap, Bp, hl)
  expect_equal(fdist(Ap, Bp, method = "hellinger"), ref, tolerance = tol)
})

test_that("chi_squared matches reference", {
  cs <- function(x, y) {
    s <- x + y
    keep <- s != 0
    0.5 * sum((x[keep] - y[keep])^2 / s[keep])
  }
  ref <- pairwise(Ap, Bp, cs)
  expect_equal(fdist(Ap, Bp, method = "chi_squared"), ref, tolerance = tol)
})

test_that("jensen_shannon matches reference and is zero on the diagonal", {
  js <- function(x, y) {
    x <- x / sum(x)
    y <- y / sum(y)
    m <- 0.5 * (x + y)
    kl <- function(p) sum(ifelse(p > 0, p * log(p / m), 0))
    sqrt(0.5 * kl(x) + 0.5 * kl(y))
  }
  ref <- pairwise(Ap, Bp, js)
  expect_equal(fdist(Ap, Bp, method = "jensen_shannon"), ref, tolerance = tol)
  d <- fdist(Ap, method = "jensen_shannon")
  expect_equal(diag(d), rep(0, nrow(Ap)), tolerance = tol)
})

test_that("haversine matches reference and requires two columns", {
  G1 <- matrix(c(40.4168, -3.7038,   # Madrid
                 41.3851,  2.1734,   # Barcelona
                 51.5074, -0.1278),  # London
               ncol = 2, byrow = TRUE)
  G2 <- matrix(c(48.8566, 2.3522),   # Paris
               ncol = 2, byrow = TRUE)
  hv <- function(x, y) {
    R <- 6371.0
    p <- pi / 180
    dlat <- (y[1] - x[1]) * p
    dlon <- (y[2] - x[2]) * p
    a <- sin(dlat / 2)^2 +
      cos(x[1] * p) * cos(y[1] * p) * sin(dlon / 2)^2
    2 * R * asin(sqrt(min(a, 1)))
  }
  ref <- pairwise(G1, G2, hv)
  expect_equal(fdist(G1, G2, method = "haversine"), ref, tolerance = 1e-6)
  # Madrid -> Paris is roughly 1050 km
  expect_true(abs(fdist(G1, G2, method = "haversine")[1, 1] - 1054) < 20)
  expect_error(fdist(A, B, method = "haversine"))
})

test_that("standardized_euclidean matches reference", {
  v <- apply(A, 2, var)
  se <- function(x, y) sqrt(sum((x - y)^2 / v))
  ref <- pairwise(A, B, se)
  expect_equal(fdist(A, B, method = "standardized_euclidean"), ref, tolerance = tol)
})

test_that("spearman matches 1 - rank correlation", {
  sp <- function(x, y) 1 - cor(x, y, method = "spearman")
  ref <- pairwise(A, B, sp)
  expect_equal(fdist(A, B, method = "spearman"), ref, tolerance = tol)
})

test_that("symmetric calls return symmetric matrices with zero diagonal", {
  for (m in c("squared_euclidean", "bray_curtis", "hellinger", "chi_squared",
              "jensen_shannon", "standardized_euclidean", "spearman")) {
    M <- if (m %in% c("bray_curtis", "hellinger", "chi_squared",
                      "jensen_shannon")) Ap else A
    d <- fdist(M, method = m)
    expect_equal(d, t(d), tolerance = tol, info = m)
    expect_equal(diag(d), rep(0, nrow(M)), tolerance = tol, info = m)
  }
})
