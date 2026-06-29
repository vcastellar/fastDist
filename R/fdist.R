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
#'   values are `"euclidean"`, `"sqeuclidean"`, `"manhattan"`, `"minkowski"`,
#'   `"correlation"`, `"cosine"`, `"canberra"`, `"supremum"`, `"hamming"`,
#'   `"jaccard"`, `"dice"`, and `"mahalanobis"`. The methods `"jaccard"` and
#'   `"dice"` require binary (0/1) input in `A` and `B`.
#' @param p Numeric scalar used only when `method = "minkowski"`. It is the
#'   exponent of the Minkowski metric (\eqn{p \ge 1} in the standard metric
#'   definition).
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
#'   \item{`"sqeuclidean"`}{
#'   Squared Euclidean distance (Euclidean without the final square root).
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k} (x_{ic} - y_{jc})^2.}
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
#'   \item{`"hamming"`}{
#'   Hamming distance: the number of coordinates in which the two rows differ.
#'   \deqn{d(x_i, y_j) = \sum_{c=1}^{k} \mathbf{1}(x_{ic} \neq y_{jc}).}
#'   }
#'
#'   \item{`"jaccard"`}{
#'   Jaccard distance for binary (0/1) vectors. With \eqn{a} the number of
#'   coordinates where both rows are 1, and \eqn{b + c} the number where exactly
#'   one row is 1,
#'   \deqn{d(x_i, y_j) = \frac{b + c}{a + b + c},}
#'   defined as `0` when \eqn{a + b + c = 0}. Requires binary input.
#'   }
#'
#'   \item{`"dice"`}{
#'   Dice distance for binary (0/1) vectors, using the same counts as Jaccard,
#'   \deqn{d(x_i, y_j) = \frac{b + c}{2a + b + c},}
#'   defined as `0` when \eqn{2a + b + c = 0}. Requires binary input.
#'   }
#'
#'   \item{`"mahalanobis"`}{
#'   Mahalanobis distance using the sample covariance matrix \eqn{S} estimated
#'   from `A` (\eqn{B} is ignored).
#'   \deqn{d(x_i, x_j) = \sqrt{(x_i - x_j)^\top S^{-1}(x_i - x_j)}.}
#'   }
#' }
#'
#' @return A numeric matrix where entry `[i, j]` is the distance between row
#'   `i` of `A` and row `j` of `B` (or row `j` of `A` when `B = NULL`).
#' @export
fdist <- function(A, B = NULL, method, p = NULL) {
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregistry"))
  }
  A <- as.matrix(A)
  if (method == "mahalanobis") {
    result <- fdistregistry$get_entry(method)$fun(A)
  } else {
    if (is.null(B)) {
      B <- A
    } else {
      B <- as.matrix(B)
    }
    if (method %in% .binary_methods) {
      .check_binary(A, "A")
      .check_binary(B, "B")
    }
    if (is.na(fdistregistry$get_entry(method)$p)) {
      result <- fdistregistry$get_entry(method)$fun(A, B)
    } else {
      result <- fdistregistry$get_entry(method)$fun(A, B, p)
    }
  }


  return(result)
}

# Distance methods that require binary (0/1) input.
.binary_methods <- c("jaccard", "dice")

# Validate that a matrix only contains binary values (0/1). Used by the
# binary distance methods. `name` is the argument name used in error messages.
.check_binary <- function(x, name) {
  if (anyNA(x)) {
    stop(paste0("'", name, "' must not contain NA values for binary methods"))
  }
  if (!all(x == 0 | x == 1)) {
    stop(paste0("'", name, "' must contain only binary values (0/1) for binary methods"))
  }
  invisible(TRUE)
}
