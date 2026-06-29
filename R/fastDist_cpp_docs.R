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

#' Internal squared Euclidean distance backend (.sqeuclidean)
#'
#' @name .sqeuclidean
#' @description
#' Computes the squared Euclidean distance between rows of `Ar` and `Br`.
#' @details
#' For two observations \eqn{x, y \in \mathbb{R}^k}, the squared Euclidean
#' distance is defined as:
#' \deqn{d(x, y) = \sum_{i=1}^{k}(x_i - y_i)^2}
#' It is the Euclidean distance without the final square root, which is cheaper
#' to compute and is often sufficient (e.g. for k-means or nearest-neighbour
#' searches). This function is an internal backend implemented in C++ and
#' exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .sqeuclidean(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Hamming distance backend (.hamming)
#'
#' @name .hamming
#' @description
#' Computes the Hamming distance between rows of `Ar` and `Br`.
#' @details
#' For two observations \eqn{x, y \in \mathbb{R}^k}, the Hamming distance counts
#' the number of coordinates in which they differ:
#' \deqn{d(x, y) = \sum_{i=1}^{k}\mathbf{1}(x_i \neq y_i)}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .hamming(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Jaccard distance backend (.jaccard)
#'
#' @name .jaccard
#' @description
#' Computes the Jaccard distance between rows of `Ar` and `Br`, treating the
#' input as binary (0/1).
#' @details
#' Let \eqn{a} be the number of coordinates where both rows equal 1 and
#' \eqn{b + c} the number where exactly one row equals 1. The Jaccard distance
#' is:
#' \deqn{d(x, y) = \frac{b + c}{a + b + c}}
#' Pairs with an empty union (\eqn{a + b + c = 0}) are assigned a distance of 0.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .jaccard(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Dice distance backend (.dice)
#'
#' @name .dice
#' @description
#' Computes the Dice distance between rows of `Ar` and `Br`, treating the input
#' as binary (0/1).
#' @details
#' Using the same counts as the Jaccard distance, the Dice distance is:
#' \deqn{d(x, y) = \frac{b + c}{2a + b + c}}
#' Pairs with \eqn{2a + b + c = 0} are assigned a distance of 0.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .dice(Ar, Br)
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
