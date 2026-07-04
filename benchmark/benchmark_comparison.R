# Reproducible benchmark: fastDist vs parallelDist, proxy and rdist
# -----------------------------------------------------------------------------
# Compares the time needed to compute the pairwise cross-distance matrix
# between the rows of two matrices A (m x k) and B (n x k) for every metric
# that fastDist shares with the packages parallelDist, proxy and rdist.
#
# How each package computes cross distances:
#   * fastDist::fdist(A, B, ...)  -- native cross distances.
#   * proxy::dist(A, B, ...)      -- native cross distances.
#   * rdist::cdist(A, B, ...)     -- native cross distances.
#   * parallelDist::parDist()     -- square matrices over one input only, so
#     we stack rbind(A, B), compute the full square matrix and extract the
#     A x B block (rows 1:m, columns (m+1):(m+n)). The extra work is part of
#     what a user would have to do in practice, so it is included in the
#     timing.
#
# Before timing, the results of every package are checked against fastDist
# (maximum absolute difference), so that we only compare implementations that
# agree on the answer.
#
# Missing packages and unsupported package/method pairs are skipped
# gracefully, so the script can be run with any subset of the four packages
# installed.
#
# Usage:
#   Rscript benchmark/benchmark_comparison.R
# or, from an interactive session:
#   source("benchmark/benchmark_comparison.R")
#   results <- run_benchmark()
#   plot_benchmark(results)      # requires ggplot2
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fastDist)
  library(microbenchmark)
})

# method names per package; NA marks an unsupported package/method pair.
# `binary` flags metrics that need 0/1 input data.
method_map <- data.frame(
  fastDist     = c("euclidean", "manhattan", "minkowski", "canberra",
                   "supremum", "hamming", "jaccard"),
  parallelDist = c("euclidean", "manhattan", "minkowski", "canberra",
                   "maximum", "hamming", "binary"),
  proxy        = c("Euclidean", "Manhattan", "Minkowski", "Canberra",
                   "Chebyshev", NA, "Jaccard"),
  rdist        = c("euclidean", "manhattan", "minkowski", "canberra",
                   "maximum", "hamming", "jaccard"),
  binary       = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# cross-distance wrappers -----------------------------------------------------

cross_fastDist <- function(A, B, method, p = 2) {
  if (method == "minkowski") {
    fdist(A, B, method = method, p = p)
  } else {
    fdist(A, B, method = method)
  }
}

cross_parallelDist <- function(A, B, method, p = 2) {
  m <- nrow(A)
  n <- nrow(B)
  M <- rbind(A, B)
  full <- if (method == "minkowski") {
    as.matrix(parallelDist::parDist(M, method = method, p = p))
  } else {
    as.matrix(parallelDist::parDist(M, method = method))
  }
  full[seq_len(m), (m + 1):(m + n), drop = FALSE]
}

cross_proxy <- function(A, B, method, p = 2) {
  res <- if (method == "Minkowski") {
    proxy::dist(A, B, method = method, p = p)
  } else {
    proxy::dist(A, B, method = method)
  }
  unclass(as.matrix(res))
}

cross_rdist <- function(A, B, method, p = 2) {
  if (method == "minkowski") {
    as.matrix(rdist::cdist(A, B, metric = method, p = p))
  } else {
    as.matrix(rdist::cdist(A, B, metric = method))
  }
}

cross_fun <- list(
  fastDist     = cross_fastDist,
  parallelDist = cross_parallelDist,
  proxy        = cross_proxy,
  rdist        = cross_rdist
)

# helpers ----------------------------------------------------------------------

available_packages <- function() {
  pkgs <- c("parallelDist", "proxy", "rdist")
  pkgs[vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
}

# Compute the cross distances for one package/method pair, returning NULL if
# the pair is unsupported or the computation fails.
try_cross <- function(pkg, A, B, method, p) {
  if (is.na(method)) return(NULL)
  tryCatch(
    cross_fun[[pkg]](A, B, method, p = p),
    error = function(e) {
      message(sprintf("  skipping %s / %s: %s", pkg, method, conditionMessage(e)))
      NULL
    }
  )
}

#' Run the benchmark.
#'
#' @param a_rows Number of rows of the fixed matrix A.
#' @param b_rows Integer vector with the number of rows of B to sweep over.
#' @param n_features Number of columns (features) of A and B.
#' @param p Order used for the Minkowski metric.
#' @param times Number of repetitions passed to microbenchmark.
#' @param methods Subset of method_map$fastDist to benchmark.
#' @param tol Tolerance used when checking that all packages agree.
#' @return A data.frame with columns package, method, mean (seconds),
#'   b_rows, a_rows, n_features and agrees (max abs difference vs fastDist).
run_benchmark <- function(a_rows = 1000,
                          b_rows = c(1000, 5000, 10000),
                          n_features = 1000,
                          p = 2,
                          times = 5,
                          methods = method_map$fastDist,
                          tol = 1e-6) {
  pkgs <- c("fastDist", available_packages())
  message("packages found: ", paste(pkgs, collapse = ", "))

  set.seed(123)
  A_num <- matrix(rnorm(a_rows * n_features), a_rows, n_features)
  A_bin <- matrix(rbinom(a_rows * n_features, 1, 0.5), a_rows, n_features)
  A_bin <- A_bin * 1.0  # ensure double storage

  results <- data.frame()

  for (nb in b_rows) {
    B_num <- matrix(rnorm(nb * n_features), nb, n_features)
    B_bin <- matrix(rbinom(nb * n_features, 1, 0.5), nb, n_features) * 1.0

    for (i in which(method_map$fastDist %in% methods)) {
      fd_method <- method_map$fastDist[i]
      is_binary <- method_map$binary[i]
      A <- if (is_binary) A_bin else A_num
      B <- if (is_binary) B_bin else B_num

      # correctness check against fastDist before timing anything
      ref <- cross_fastDist(A, B, fd_method, p = p)
      run_pkgs <- character(0)
      max_diff <- c(fastDist = 0)
      for (pkg in setdiff(pkgs, "fastDist")) {
        out <- try_cross(pkg, A, B, method_map[[pkg]][i], p)
        if (is.null(out)) next
        d <- max(abs(out - ref))
        max_diff[pkg] <- d
        if (d <= tol) {
          run_pkgs <- c(run_pkgs, pkg)
        } else {
          message(sprintf(
            "  %s / %s disagrees with fastDist (max abs diff = %.3g), skipping",
            pkg, fd_method, d
          ))
        }
      }

      exprs <- c(
        list(fastDist = bquote(cross_fastDist(A, B, .(fd_method), p = .(p)))),
        setNames(lapply(run_pkgs, function(pkg) {
          bquote(cross_fun[[.(pkg)]](A, B, .(method_map[[pkg]][i]), p = .(p)))
        }), run_pkgs)
      )

      mb <- microbenchmark(list = exprs, times = times, unit = "s")
      means <- aggregate(time ~ expr, data = mb, FUN = mean)
      means$time <- means$time / 1e9  # nanoseconds -> seconds

      results <- rbind(results, data.frame(
        package = as.character(means$expr),
        method = fd_method,
        mean = means$time,
        b_rows = nb,
        a_rows = a_rows,
        n_features = n_features,
        agrees = max_diff[as.character(means$expr)],
        stringsAsFactors = FALSE,
        row.names = NULL
      ))
      message(sprintf("done: method = %-10s b_rows = %d", fd_method, nb))
    }
  }

  results
}

# Relative speed of every package vs fastDist (values > 1 mean fastDist is
# faster).
speedup_table <- function(results) {
  base <- results[results$package == "fastDist",
                  c("method", "b_rows", "mean")]
  names(base)[names(base) == "mean"] <- "fastDist_mean"
  merged <- merge(results, base, by = c("method", "b_rows"))
  merged$speedup <- merged$mean / merged$fastDist_mean
  merged[order(merged$method, merged$b_rows, merged$package),
         c("method", "b_rows", "package", "mean", "speedup")]
}

# Optional plot of the results (requires ggplot2).
plot_benchmark <- function(results, file = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_benchmark()")
  }
  g <- ggplot2::ggplot(
    results,
    ggplot2::aes(b_rows, mean, color = package)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      x = "rows of B", y = "mean time (s)", color = "package",
      title = "Cross-distance computation: fastDist vs parallelDist, proxy, rdist"
    )
  if (!is.null(file)) {
    ggplot2::ggsave(file, g, width = 12, height = 7)
  }
  g
}

# Run when executed as a script (e.g. via Rscript).
if (sys.nframe() == 0) {
  results <- run_benchmark()
  print(results)
  cat("\nSpeed relative to fastDist (>1 means fastDist is faster):\n")
  print(speedup_table(results))
}
