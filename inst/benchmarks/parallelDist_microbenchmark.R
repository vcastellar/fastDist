# Microbenchmark: fastDist vs parallelDist for A (1000 x 1000) against B (n x 1000)
#
# Example:
# install.packages(c("microbenchmark", "parallelDist"))
# devtools::load_all(".")
# source("inst/benchmarks/parallelDist_microbenchmark.R")

benchmark_result <- benchmark_parallelDist(
  a_nrow = 1000L,
  n_features = 1000L,
  # b_sizes = c(1000L, 5000L, 10000L, 20000L),
  b_sizes = c(1000L, 5000L, 10000L, 20000L),
  minkowski_p = 3,
  times = 1L,
  unit = "s",
  block_size = 1000L,
  check = FALSE
)

print(benchmark_result$summary)[, c("expr", "mean", "method", "b_rows", "a_rows", "n_features")]
