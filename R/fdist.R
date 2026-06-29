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
#'   `"haversine"`, `"standardized_euclidean"`, `"spearman"`, and
#'   `"mahalanobis"`.
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
    if (is.na(fdistregistry$get_entry(method)$p)) {
      result <- fdistregistry$get_entry(method)$fun(A, B)
    } else {
      result <- fdistregistry$get_entry(method)$fun(A, B, p)
    }
  }

  
  return(result)
}
