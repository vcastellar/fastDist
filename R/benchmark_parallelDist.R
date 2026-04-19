#' Validate distance method names
#'
#' Checks that `method` belongs to the set of supported methods for
#' `parallelDist` benchmarking.
#'
#' @param method Character scalar with the distance method name.
#'
#' @return Invisibly returns `NULL`; throws an error if validation fails.
#' @keywords internal
#' @noRd
.validate_distance_method <- function(method) {
  valid_methods <- c("euclidean", "manhattan", "minkowski")
  if (!method %in% valid_methods) {
    stop(sprintf(
      "method must be one of %s",
      paste(shQuote(valid_methods), collapse = ", ")
    ))
  }
}

#' Ensure benchmark dependencies are available
#'
#' Verifies that required benchmarking packages are installed.
#'
#' @return Invisibly returns `NULL`; throws an error when dependencies are
#'   missing.
#' @keywords internal
#' @noRd
.require_benchmark_packages <- function() {
  missing_packages <- c()

  if (!requireNamespace("parallelDist", quietly = TRUE)) {
    missing_packages <- c(missing_packages, "parallelDist")
  }

  if (!requireNamespace("microbenchmark", quietly = TRUE)) {
    missing_packages <- c(missing_packages, "microbenchmark")
  }

  if (length(missing_packages) > 0) {
    stop(
      sprintf(
        "Missing required package(s): %s",
        paste(missing_packages, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

#' Compute Cross Distances with `parallelDist` in Blocks
#'
#' Computes pairwise distances between rows of `A` and rows of `B` using
#' `parallelDist::parDist()`, processing `B` by blocks to reduce memory
#' pressure.
#'
#' @param A Numeric matrix (or matrix-like object) of source observations.
#' @param B Numeric matrix (or matrix-like object) of target observations.
#' @param method Distance method. Supported values are `"euclidean"`,
#'   `"manhattan"`, and `"minkowski"`.
#' @param p Minkowski exponent. Required when `method = "minkowski"`.
#' @param threads Optional number of threads passed to
#'   `parallelDist::parDist()`.
#' @param block_size Number of rows of `B` processed per iteration. If `NULL`,
#'   defaults to `nrow(A)`.
#'
#' @return A numeric matrix of size `nrow(A) x nrow(B)` with cross distances.
#' @export
crossdist_parallelDist <- function(A,
                                   B,
                                   method = "euclidean",
                                   p = NULL,
                                   threads = NULL,
                                   block_size = NULL) {
  if (!requireNamespace("parallelDist", quietly = TRUE)) {
    stop("Package 'parallelDist' is required for crossdist_parallelDist().", call. = FALSE)
  }

  .validate_distance_method(method)

  A <- as.matrix(A)
  B <- as.matrix(B)

  if (ncol(A) != ncol(B)) {
    stop("A and B must have the same number of columns.", call. = FALSE)
  }

  if (is.null(block_size)) {
    block_size <- max(nrow(A), 1L)
  }

  block_size <- as.integer(block_size)
  if (is.na(block_size) || block_size < 1L) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }

  n_a <- nrow(A)
  n_b <- nrow(B)
  result <- matrix(NA_real_, nrow = n_a, ncol = n_b)

  for (start in seq.int(1L, n_b, by = block_size)) {
    end <- min(start + block_size - 1L, n_b)
    index <- start:end
    combined <- rbind(A, B[index, , drop = FALSE])

    call_args <- list(
      x = combined,
      method = method,
      threads = threads
    )

    if (identical(method, "minkowski")) {
      if (is.null(p)) {
        stop("Parameter 'p' is required when method = 'minkowski'.", call. = FALSE)
      }
      call_args$p <- p
    }

    block_dist <- do.call(parallelDist::parDist, call_args)
    block_matrix <- as.matrix(block_dist)
    result[, index] <- block_matrix[seq_len(n_a), n_a + seq_along(index), drop = FALSE]
  }

  result
}

#' Benchmark `fastDist` vs `parallelDist`
#'
#' Runs benchmarking experiments comparing [`fdist()`] and
#' [`crossdist_parallelDist()`] across multiple methods and `B` matrix sizes.
#' Optional correctness checks can verify numerical agreement before timing.
#'
#' @param A Optional numeric matrix of source observations. If `NULL`, a random
#'   matrix is generated.
#' @param B Optional numeric matrix of target observations. If `NULL`, a random
#'   matrix is generated.
#' @param a_nrow Number of rows to generate for `A` when `A` is `NULL`.
#' @param n_features Number of columns/features for generated matrices.
#' @param b_sizes Integer vector with the evaluated numbers of rows for `B`.
#' @param methods Character vector of methods to benchmark. Supported values are
#'   `"euclidean"`, `"manhattan"`, and `"minkowski"`.
#' @param minkowski_p Exponent used when benchmarking the Minkowski distance.
#' @param times Number of repetitions passed to `microbenchmark`.
#' @param seed Random seed used when generating `A` and `B`.
#' @param unit Time unit for benchmark results (passed to `microbenchmark`).
#' @param threads Optional number of threads for `parallelDist`.
#' @param block_size Optional block size used by [`crossdist_parallelDist()`].
#' @param check Logical; if `TRUE`, validate that both implementations produce
#'   equivalent results before timing.
#'
#' @return A list with benchmark parameters, a summary table, raw benchmark
#'   objects, and the matrices used in the benchmark run.
#' @export
benchmark_parallelDist <- function(A = NULL,
                                   B = NULL,
                                   a_nrow = 1000L,
                                   n_features = 1000L,
                                   b_sizes = c(1000L, 5000L, 10000L, 20000L),
                                   methods = c("euclidean", "manhattan", "minkowski"),
                                   minkowski_p = 3,
                                   times = 1L,
                                   seed = 123,
                                   unit = "s",
                                   threads = NULL,
                                   block_size = NULL,
                                   check = TRUE) {
  .require_benchmark_packages()

  methods <- unique(methods)
  invalid_methods <- setdiff(methods, c("euclidean", "manhattan", "minkowski"))
  if (length(invalid_methods) > 0) {
    stop(
      sprintf(
        "Unsupported method(s): %s",
        paste(shQuote(invalid_methods), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  b_sizes <- sort(unique(as.integer(b_sizes)))
  if (any(is.na(b_sizes)) || any(b_sizes < 1L)) {
    stop("b_sizes must contain positive integers.", call. = FALSE)
  }

  a_nrow <- as.integer(a_nrow)
  n_features <- as.integer(n_features)
  times <- as.integer(times)

  if (is.null(A)) {
    set.seed(seed)
    A <- matrix(runif(a_nrow * n_features), nrow = a_nrow, ncol = n_features)
  } else {
    A <- as.matrix(A)
    a_nrow <- nrow(A)
    n_features <- ncol(A)
  }

  max_b <- max(b_sizes)
  if (is.null(B)) {
    set.seed(seed + 1L)
    B <- matrix(runif(max_b * n_features), nrow = max_b, ncol = n_features)
  } else {
    B <- as.matrix(B)
    if (ncol(B) != n_features) {
      stop("B must have the same number of columns as A.", call. = FALSE)
    }
    if (nrow(B) < max_b) {
      stop("B must have at least max(b_sizes) rows.", call. = FALSE)
    }
  }

  if (is.null(block_size)) {
    block_size <- max(nrow(A), 1L)
  }

  benchmark_rows <- list()
  benchmark_objects <- list()
  benchmark_id <- 1L

  for (method in methods) {
    for (b_size in b_sizes) {
      current_B <- B[seq_len(b_size), , drop = FALSE]

      fastdist_call <- function() {
        if (identical(method, "minkowski")) {
          fdist(A, current_B, method = method, p = minkowski_p)
        } else {
          fdist(A, current_B, method = method)
        }
      }

      paralleldist_call <- function() {
        if (identical(method, "minkowski")) {
          crossdist_parallelDist(
            A,
            current_B,
            method = method,
            p = minkowski_p,
            threads = threads,
            block_size = block_size
          )
        } else {
          crossdist_parallelDist(
            A,
            current_B,
            method = method,
            threads = threads,
            block_size = block_size
          )
        }
      }

      if (isTRUE(check)) {
        fastdist_result <- fastdist_call()
        paralleldist_result <- paralleldist_call()

        if (!isTRUE(all.equal(fastdist_result, paralleldist_result, tolerance = 1e-8))) {
          stop(
            sprintf(
              "The results differ for method '%s' and b_size = %s.",
              method,
              b_size
            ),
            call. = FALSE
          )
        }
      }

      benchmark <- microbenchmark::microbenchmark(
        fastDist = fastdist_call(),
        parallelDist = paralleldist_call(),
        times = times,
        unit = unit,
        check = NULL
      )

      benchmark_summary <- summary(benchmark)
      benchmark_summary$method <- method
      benchmark_summary$b_rows <- b_size
      benchmark_summary$a_rows <- nrow(A)
      benchmark_summary$n_features <- ncol(A)
      benchmark_summary$minkowski_p <- if (identical(method, "minkowski")) minkowski_p else NA_real_
      benchmark_rows[[benchmark_id]] <- benchmark_summary
      benchmark_objects[[paste(method, b_size, sep = "__")]] <- benchmark
      benchmark_id <- benchmark_id + 1L
    }
  }

  summary_table <- do.call(rbind, benchmark_rows)
  rownames(summary_table) <- NULL

  list(
    parameters = list(
      a_rows = nrow(A),
      n_features = ncol(A),
      b_sizes = b_sizes,
      methods = methods,
      minkowski_p = minkowski_p,
      times = times,
      unit = unit,
      threads = threads,
      block_size = block_size,
      check = check
    ),
    summary = summary_table,
    benchmarks = benchmark_objects,
    A = A,
    B = B[seq_len(max_b), , drop = FALSE]
  )
}
