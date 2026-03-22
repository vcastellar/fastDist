test_that("crossdist_parallelDist matches fastDist on small matrices", {
  skip_if_not_installed("parallelDist")

  set.seed(42)
  A <- matrix(runif(12), nrow = 4)
  B <- matrix(runif(21), nrow = 7)

  expect_equal(
    crossdist_parallelDist(A, B, method = "euclidean", block_size = 3L),
    fdist(A, B, method = "euclidean"),
    tolerance = 1e-8
  )

  expect_equal(
    crossdist_parallelDist(A, B, method = "manhattan", block_size = 3L),
    fdist(A, B, method = "manhattan"),
    tolerance = 1e-8
  )

  expect_equal(
    crossdist_parallelDist(A, B, method = "minkowski", p = 3, block_size = 3L),
    fdist(A, B, method = "minkowski", p = 3),
    tolerance = 1e-8
  )
})

test_that("benchmark_parallelDist returns a summary for each method and B size", {
  skip_if_not_installed("parallelDist")
  skip_if_not_installed("microbenchmark")

  set.seed(99)
  A <- matrix(runif(24), nrow = 6)
  B <- matrix(runif(80), nrow = 20)

  benchmark_result <- benchmark_parallelDist(
    A = A,
    B = B,
    b_sizes = c(5L, 10L),
    methods = c("euclidean", "minkowski"),
    minkowski_p = 3,
    times = 1L,
    unit = "ms",
    block_size = 5L,
    check = TRUE
  )

  expect_true(is.list(benchmark_result))
  expect_named(
    benchmark_result,
    c("parameters", "summary", "benchmarks", "A", "B")
  )
  expect_equal(
    sort(unique(benchmark_result$summary$expr)),
    c("fastDist", "parallelDist")
  )
  expect_equal(length(unique(benchmark_result$summary$method)), 2L)
  expect_equal(length(unique(benchmark_result$summary$b_rows)), 2L)
})
