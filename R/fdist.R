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
#'   `"cosine"`, `"canberra"`, `"supremum"`, and `"mahalanobis"`.
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
    stop(paste(method, "not found in fdistregestry"))
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
    if (is.na(fdistregistry$get_entry(method)$p)) {
      result <- fdistregistry$get_entry(method)$fun(A, B)
    } else {
      result <- fdistregistry$get_entry(method)$fun(A, B, p)
    }
  }

  
  return(result)
}
