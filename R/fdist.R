#' Fast Pairwise Distance Computation Between Observations
#'
#' Computes pairwise distances between rows of `A` and `B` using distance
#' backends registered in `fdistregistry`.
#'
#' @param A A matrix (or object coercible to a matrix) containing source
#'   observations in rows.
#' @param B A matrix (or object coercible to a matrix) containing target
#'   observations in rows. If `NULL`, `A` is used (except for
#'   `method = "mahalanobis"`, which always uses only `A`).
#' @param method Character scalar with the distance method to use. Supported
#'   values are `"euclidean"`, `"manhattan"`, `"minkowski"`, `"correlation"`,
#'   `"cosine"`, `"canberra"`, `"supremum"`, `"squared_euclidean"`,
#'   `"bray_curtis"`, `"hellinger"`, `"chi_squared"`, `"jensen_shannon"`,
#'   `"haversine"`, `"standardized_euclidean"`, `"spearman"`, `"mahalanobis"`,
#'   `"hamming"`, `"jaccard"`, and `"gower"`.
#' @param ... Extra, method-specific parameters overriding the defaults
#'   stored in [fdistregistry] (see `fdistregistry$get_entry(method)$params`
#'   for the defaults of a given method). Recognized names are:
#'   \describe{
#'     \item{`p`}{(`method = "minkowski"`) exponent of the Minkowski metric
#'     (\eqn{p \ge 1} in the standard metric definition). Defaults to `2`.}
#'     \item{`radius`}{(`method = "haversine"`) sphere radius the result is
#'     expressed in (e.g. `6371` for kilometres or `3958.8` for miles).
#'     Defaults to `6371`.}
#'     \item{`base`}{(`method = "jensen_shannon"`) either the character
#'     shortcuts `"e"` / `"2"`, or a numeric logarithm base (\eqn{b > 0},
#'     \eqn{b \neq 1}), controlling the units of the underlying divergence
#'     (nats for `"e"`, bits for `"2"`). Defaults to `exp(1)` (natural log).}
#'     \item{`threshold`}{(`method = "jaccard"` or `"hamming"`) values
#'     strictly greater than `threshold` are treated as `1` (present)
#'     before comparing. Defaults to `0` for `"jaccard"`; for `"hamming"`
#'     it defaults to `NA`, which disables binarization and compares raw
#'     values for exact inequality instead.}
#'     \item{`regularize`}{(`method = "mahalanobis"`) ridge added to the
#'     diagonal of the covariance matrix before inversion, useful when it
#'     is singular or near-singular. Defaults to `0`.}
#'     \item{`weights`}{(`method = "standardized_euclidean"`) numeric
#'     vector with the per-feature scale used instead of the sample
#'     variance of `A` (its length must equal `ncol(A)`). Defaults to
#'     `NULL`, i.e. the sample variance of `A` is used.}
#'     \item{`cov`}{(`method = "mahalanobis"`) numeric matrix (`ncol(A)` x
#'     `ncol(A)`) with a precomputed covariance matrix used instead of the
#'     sample covariance of `A`. Defaults to `NULL`.}
#'   }
#'
#' @details
#' Let \eqn{x_i = (x_{i1}, \ldots, x_{ik})} be row `i` from `A` and
#' \eqn{y_j = (y_{j1}, \ldots, y_{jk})} be row `j` from `B`.
#'
#' Available `method` values are:
#'
#' \describe{
#'   \item{`"euclidean"`}{
#'   Standard \eqn{L_2} distance.
#'   \deqn{d(x_i, y_j) = \sqrt{\sum_{c=1}^{k} (x_{ic} - y_{jc})^2}.}
#'   }
#'
#'   \item{`"manhattan"`}{
#'   City-block (\eqn{L_1}) distance.
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k} |x_{ic} - y_{jc}|.}
#'   }
#'
#'   \item{`"minkowski"`}{
#'   Generalized \eqn{L_p} distance controlled by `p`.
#'   \deqn{d(x_i, y_j) = \left(\sum_{c=1}^{k} |x_{ic} - y_{jc}|^p\right)^{1/p}.}
#'   Special cases: `p = 1` gives Manhattan and `p = 2` gives Euclidean.
#'   }
#'
#'   \item{`"correlation"`}{
#'   Correlation distance, defined as one minus the Pearson correlation between
#'   centered row vectors.
#'   \deqn{d(x_i, y_j) = 1 - \frac{\sum_{c=1}^{k}(x_{ic}-\bar{x}_i)(y_{jc}-\bar{y}_j)}
#'   {\sqrt{\sum_{c=1}^{k}(x_{ic}-\bar{x}_i)^2}\sqrt{\sum_{c=1}^{k}(y_{jc}-\bar{y}_j)^2}}.}
#'   }
#'
#'   \item{`"cosine"`}{
#'   Cosine distance, defined as one minus cosine similarity.
#'   \deqn{d(x_i, y_j) = 1 - \frac{\sum_{c=1}^{k}x_{ic}y_{jc}}
#'   {\sqrt{\sum_{c=1}^{k}x_{ic}^2}\sqrt{\sum_{c=1}^{k}y_{jc}^2}}.}
#'   }
#'
#'   \item{`"canberra"`}{
#'   Canberra distance (feature-wise normalized absolute differences).
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k}\frac{|x_{ic}-y_{jc}|}{|x_{ic}|+|y_{jc}|},}
#'   where terms with zero denominator contribute `0`.
#'   }
#'
#'   \item{`"supremum"`}{
#'   Supremum (Chebyshev, \eqn{L_\infty}) distance.
#'   \deqn{d(x_i, y_j) = \max_{c=1,\ldots,k}|x_{ic}-y_{jc}|.}
#'   }
#'
#'   \item{`"squared_euclidean"`}{
#'   Squared Euclidean distance (Euclidean without the final square root),
#'   useful for clustering algorithms such as k-means.
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k} (x_{ic} - y_{jc})^2.}
#'   }
#'
#'   \item{`"bray_curtis"`}{
#'   Bray-Curtis dissimilarity, common for abundance/compositional data.
#'   \deqn{d(x_i, y_j) = \frac{\sum_{c=1}^{k}|x_{ic}-y_{jc}|}
#'   {\sum_{c=1}^{k}|x_{ic}+y_{jc}|},}
#'   where a zero denominator yields `0`.
#'   }
#'
#'   \item{`"hellinger"`}{
#'   Hellinger distance between non-negative vectors (e.g. distributions).
#'   \deqn{d(x_i, y_j) = \frac{1}{\sqrt{2}}
#'   \sqrt{\sum_{c=1}^{k}(\sqrt{x_{ic}}-\sqrt{y_{jc}})^2}.}
#'   }
#'
#'   \item{`"chi_squared"`}{
#'   Symmetric chi-squared distance for histograms / non-negative data.
#'   \deqn{d(x_i, y_j) = \frac{1}{2}\sum_{c=1}^{k}
#'   \frac{(x_{ic}-y_{jc})^2}{x_{ic}+y_{jc}},}
#'   where terms with zero denominator contribute `0`.
#'   }
#'
#'   \item{`"jensen_shannon"`}{
#'   Jensen-Shannon distance (square root of the Jensen-Shannon divergence,
#'   natural logarithm). Each row is normalized to sum to one beforehand. With
#'   \eqn{m = (x_i + y_j)/2},
#'   \deqn{d(x_i, y_j) = \sqrt{\tfrac{1}{2}\!\sum_c x_{ic}\log\frac{x_{ic}}{m_c}
#'   + \tfrac{1}{2}\!\sum_c y_{jc}\log\frac{y_{jc}}{m_c}}.}
#'   }
#'
#'   \item{`"haversine"`}{
#'   Great-circle distance in kilometres between points given as
#'   `(latitude, longitude)` in degrees (requires exactly two columns; Earth
#'   radius 6371 km).
#'   }
#'
#'   \item{`"standardized_euclidean"`}{
#'   Euclidean distance with each feature scaled by its sample variance
#'   \eqn{s_c^2} estimated from `A`.
#'   \deqn{d(x_i, y_j) = \sqrt{\sum_{c=1}^{k}\frac{(x_{ic}-y_{jc})^2}{s_c^2}},}
#'   where features with zero variance are dropped.
#'   }
#'
#'   \item{`"spearman"`}{
#'   Spearman correlation distance, i.e. one minus Spearman's rank correlation
#'   computed between the within-row average ranks of `x_i` and `y_j`.
#'   \deqn{d(x_i, y_j) = 1 - \rho_s(x_i, y_j).}
#'   }
#'
#'   \item{`"mahalanobis"`}{
#'   Mahalanobis distance using the sample covariance matrix \eqn{S} estimated
#'   from `A` (\eqn{B} is ignored).
#'   \deqn{d(x_i, x_j) = \sqrt{(x_i - x_j)^\top S^{-1}(x_i - x_j)}.}
#'   }
#'
#'   \item{`"hamming"`}{
#'   Hamming distance, counts the number of positions where values differ.
#'   Useful for categorical/binary data.
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k} \mathbb{1}(x_{ic} \neq y_{jc}).}
#'   }
#'
#'   \item{`"jaccard"`}{
#'   Jaccard distance for binary/set-based data (treating non-zero values as 1).
#'   \deqn{d(x_i, y_j) = 1 - \frac{|x_i \cap y_j|}{|x_i \cup y_j|},}
#'   where intersection and union are computed treating non-zero elements as 1.
#'   }
#'
#'   \item{`"gower"`}{
#'   Gower distance for mixed data types. Scales each feature by its range and
#'   averages the absolute differences.
#'   \deqn{d(x_i, y_j) = \frac{1}{k}\sum_{c=1}^{k}\frac{|x_{ic}-y_{jc}|}{R_c},}
#'   where \eqn{R_c} is the range (max - min) of feature \eqn{c} in `A`.
#'   }
#' }
#'
#' @return A numeric matrix where entry `[i, j]` is the distance between row
#'   `i` of `A` and row `j` of `B` (or row `j` of `A` when `B = NULL`).
#'
#' @examples
#' set.seed(1)
#' A <- matrix(rnorm(5 * 3), nrow = 5, ncol = 3)
#' B <- matrix(rnorm(4 * 3), nrow = 4, ncol = 3)
#'
#' # cross distances between the rows of A and the rows of B
#' fdist(A, B, method = "euclidean")
#'
#' # distances within the rows of A (symmetric, zero diagonal)
#' fdist(A, method = "manhattan")
#'
#' # Minkowski distance of order p = 3
#' fdist(A, B, method = "minkowski", p = 3)
#'
#' # Jensen-Shannon distance in bits (base 2) instead of nats
#' fdist(A, method = "jensen_shannon", base = "2")
#'
#' # Mahalanobis distance with a regularized (ridge) covariance matrix
#' fdist(A, method = "mahalanobis", regularize = 0.01)
#'
#' # binary data: Hamming and Jaccard distances
#' X <- matrix(rbinom(5 * 8, 1, 0.5), nrow = 5, ncol = 8)
#' fdist(X, method = "hamming")
#' fdist(X, method = "jaccard")
#'
#' # continuous data binarized at a custom threshold before comparing
#' fdist(A, method = "hamming", threshold = 0.5)
#'
#' # mixed-scale numeric data: Gower distance
#' fdist(A, method = "gower")
#'
#' # geographic coordinates (latitude, longitude) in degrees
#' cities <- rbind(c(40.4168, -3.7038),  # Madrid
#'                 c(41.3851,  2.1734),  # Barcelona
#'                 c(48.8566,  2.3522))  # Paris
#' fdist(cities, method = "haversine")
#'
#' # Haversine distance in miles instead of the default kilometres
#' fdist(cities, method = "haversine", radius = 3958.8)
#'
#' # all available methods
#' fdistregistry$get_entry_names()
#'
#' @seealso [fdistregistry] for the registry of available distance backends.
#' @export
fdist <- function(A, B = NULL, method, ...) {
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregistry"))
  }
  A <- as.matrix(A)
  entry <- fdistregistry$get_entry(method)
  params <- utils::modifyList(entry$params, list(...))

  if (identical(method, "jensen_shannon") && !is.null(params$base)) {
    params$base <- resolve_log_base(params$base)
  }

  if (method == "mahalanobis") {
    result <- do.call(entry$fun, c(list(A), params))
  } else {
    if (is.null(B)) {
      B <- A
    } else {
      B <- as.matrix(B)
    }
    result <- do.call(entry$fun, c(list(A, B), params))
  }

  result
}

# resolves the `base` argument of fdist() into a numeric logarithm base
resolve_log_base <- function(base) {
  if (is.character(base)) {
    base <- switch(base,
                    "e" = exp(1),
                    "2" = 2,
                    stop('base must be "e", "2", or a positive number other than 1'))
  }
  if (!is.numeric(base) || length(base) != 1L || is.na(base) ||
      base <= 0 || base == 1) {
    stop('base must be "e", "2", or a positive number other than 1')
  }
  base
}
