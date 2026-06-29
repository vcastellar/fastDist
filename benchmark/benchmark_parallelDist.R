# Reproducible benchmark: fastDist::fdist() vs parallelDist::parDist()
# -----------------------------------------------------------------------------
# Compares the time needed to compute the pairwise distances between the rows
# of two matrices A and B for the methods Euclidean, Manhattan and Minkowski.
#
# parallelDist::parDist() computes *square* distance matrices over the rows of
# a single matrix. To obtain the cross distances between A (m x k) and B
# (n x k), we stack rbind(A, B), compute the full square distance matrix and
# extract the A x B submatrix (rows 1:m, columns (m+1):(m+n)). This mirrors the
# block strategy described in the package README and lets us study how the time
# scales as the number of rows of B grows.
#
# Usage:
#   Rscript benchmark/benchmark_parallelDist.R
# or, from an interactive session:
#   source("benchmark/benchmark_parallelDist.R")
#   results <- run_benchmark()
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fastDist)
  library(parallelDist)
  library(microbenchmark)
})

# fastDist cross distances A x B for a given method.
fastdist_cross <- function(A, B, method, p = NULL) {
  if (method == "minkowski") {
    fdist(A, B, method = method, p = p)
  } else {
    fdist(A, B, method = method)
  }
}

# parallelDist cross distances A x B, extracted from the square matrix of
# rbind(A, B). `pd_method` is the parallelDist name of the metric.
pardist_cross <- function(A, B, pd_method, p = NULL) {
  m <- nrow(A)
  n <- nrow(B)
  M <- rbind(A, B)
  if (pd_method == "minkowski") {
    full <- as.matrix(parDist(M, method = "minkowski", p = p))
  } else {
    full <- as.matrix(parDist(M, method = pd_method))
  }
  full[seq_len(m), (m + 1):(m + n), drop = FALSE]
}

# methods to compare: name used by fdist() and the matching parallelDist name.
methods_tbl <- data.frame(
  fastDist     = c("euclidean", "manhattan", "minkowski"),
  parallelDist = c("euclidean", "manhattan", "minkowski"),
  stringsAsFactors = FALSE
)

#' Run the benchmark.
#'
#' @param a_rows Number of rows of the fixed matrix A.
#' @param b_rows Integer vector with the number of rows of B to sweep over.
#' @param n_features Number of columns (features) of A and B.
#' @param p Order used for the Minkowski metric.
#' @param times Number of repetitions passed to microbenchmark.
#' @return A data.frame with the mean time (seconds) per expr/method/size.
run_benchmark <- function(a_rows = 1000,
                          b_rows = c(1000, 5000, 10000, 20000),
                          n_features = 1000,
                          p = 2,
                          times = 5) {
  set.seed(123)
  A <- matrix(rnorm(a_rows * n_features), a_rows, n_features)

  results <- data.frame(
    expr = character(), method = character(), mean = numeric(),
    b_rows = integer(), a_rows = integer(), n_features = integer(),
    stringsAsFactors = FALSE
  )

  for (nb in b_rows) {
    B <- matrix(rnorm(nb * n_features), nb, n_features)

    for (i in seq_len(nrow(methods_tbl))) {
      fd_m <- methods_tbl$fastDist[i]
      pd_m <- methods_tbl$parallelDist[i]

      res <- microbenchmark(
        fastDist     = fastdist_cross(A, B, fd_m, p = p),
        parallelDist = pardist_cross(A, B, pd_m, p = p),
        times = times,
        unit = "s"
      )
      means <- aggregate(time ~ expr, data = res, FUN = mean)
      means$time <- means$time / 1e9 # nanoseconds -> seconds

      for (j in seq_len(nrow(means))) {
        results <- rbind(results, data.frame(
          expr = as.character(means$expr[j]),
          method = fd_m,
          mean = means$time[j],
          b_rows = nb,
          a_rows = a_rows,
          n_features = n_features,
          stringsAsFactors = FALSE
        ))
      }
    }
    cat(sprintf("done: b_rows = %d\n", nb))
  }

  results
}

# Optional plot of the results (requires ggplot2).
plot_benchmark <- function(results, file = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_benchmark()")
  }
  g <- ggplot2::ggplot(
    results,
    ggplot2::aes(b_rows, mean, color = expr)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      x = "rows of B", y = "mean time (s)", color = "package",
      title = "fastDist vs parallelDist"
    )
  if (!is.null(file)) {
    ggplot2::ggsave(file, g, width = 12, height = 5)
  }
  g
}

# Run when executed as a script (e.g. via Rscript).
if (sys.nframe() == 0) {
  results <- run_benchmark()
  print(results)
}
