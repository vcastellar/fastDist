#' Internal Euclidean distance backend (.euclidean)
#'
#' @name .euclidean
#' @description
#' Computes the Euclidean distance between each pair of rows in `Ar` and `Br`.
#' @details
#' For two observations \eqn{x, y \in \mathbb{R}^k}, the Euclidean distance is
#' defined as:
#' \deqn{d(x, y) = \sqrt{\sum_{i=1}^{k}(x_i - y_i)^2}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @param Ar Numeric matrix of size `m x k` with observations in rows.
#' @param Br Numeric matrix of size `n x k` with observations in rows.
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .euclidean(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Manhattan distance backend (.manhattan)
#'
#' @name .manhattan
#' @description
#' Computes the Manhattan distance (\eqn{L_1} norm) between rows of `Ar` and `Br`.
#' @details
#' For two observations \eqn{x, y \in \mathbb{R}^k}, the Manhattan distance is
#' defined as:
#' \deqn{d(x, y) = \sum_{i=1}^{k}|x_i - y_i|}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .manhattan(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Minkowski distance backend (.minkowski)
#'
#' @name .minkowski
#' @description
#' Computes the Minkowski distance between rows of `Ar` and `Br` for a given `p`.
#' @details
#' For \eqn{p > 0} and observations \eqn{x, y \in \mathbb{R}^k}, the Minkowski
#' distance is defined as:
#' \deqn{d(x, y) = \left(\sum_{i=1}^{k}|x_i - y_i|^p\right)^{1/p}}
#' Special cases include \eqn{p=1} (Manhattan) and \eqn{p=2} (Euclidean).
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#' @param p Order of the Minkowski distance, with `p > 0`.
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .minkowski(Ar, Br, p)
#' @keywords internal
#' @noRd
NULL

#' Internal correlation distance backend (.correlation)
#'
#' @name .correlation
#' @description
#' Computes the correlation distance between rows of `Ar` and `Br`.
#' @details
#' Correlation distance is defined as one minus the Pearson correlation:
#' \deqn{d(x, y) = 1 - \mathrm{corr}(x, y)}
#' where
#' \deqn{\mathrm{corr}(x, y) =
#' \frac{\sum_{i=1}^{k}(x_i - \bar{x})(y_i - \bar{y})}
#' {\sqrt{\sum_{i=1}^{k}(x_i - \bar{x})^2}\sqrt{\sum_{i=1}^{k}(y_i - \bar{y})^2}}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .correlation(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal cosine distance backend (.cosine)
#'
#' @name .cosine
#' @description
#' Computes the cosine distance between rows of `Ar` and `Br`.
#' @details
#' Cosine distance is defined as:
#' \deqn{d(x, y) = 1 - \frac{x \cdot y}{\|x\|_2\|y\|_2}}
#' where \eqn{x \cdot y} is the dot product and \eqn{\|\cdot\|_2} is the
#' Euclidean norm.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .cosine(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Canberra distance backend (.canberra)
#'
#' @name .canberra
#' @description
#' Computes the Canberra distance between rows of `Ar` and `Br`.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}, the Canberra distance is defined
#' as:
#' \deqn{d(x, y) = \sum_{i=1}^{k}
#' \frac{|x_i - y_i|}{|x_i| + |y_i|}}
#' In the implementation, terms with a zero denominator are skipped.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .canberra(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal supremum distance backend (.supremum)
#'
#' @name .supremum
#' @description
#' Computes the supremum (Chebyshev, \eqn{L_\infty}) distance between rows of
#' `Ar` and `Br`.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}, the supremum distance is defined
#' as:
#' \deqn{d(x, y) = \max_{i=1,\ldots,k}|x_i - y_i|}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .supremum(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Mahalanobis distance backend (.mahalanobis)
#'
#' @name .mahalanobis
#' @description
#' Computes the Mahalanobis distance matrix between all rows in `Ar`.
#' @details
#' Let \eqn{S} be the sample covariance matrix of `Ar`, and \eqn{S^{-1}} its
#' inverse. For two observations \eqn{x, y \in \mathbb{R}^k}, Mahalanobis
#' distance is defined as:
#' \deqn{d(x, y) = \sqrt{(x - y)^\top S^{-1}(x - y)}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @param Ar Numeric matrix of size `m x k` with observations in rows.
#'
#' @return Numeric square `m x m` matrix with pairwise distances.
#' @usage .mahalanobis(Ar)
#' @keywords internal
#' @noRd
NULL
