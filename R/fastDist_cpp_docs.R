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

#' Internal squared Euclidean distance backend (.squared_euclidean)
#'
#' @name .squared_euclidean
#' @description
#' Computes the squared Euclidean distance between rows of `Ar` and `Br`.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}:
#' \deqn{d(x, y) = \sum_{i=1}^{k}(x_i - y_i)^2}
#' This is the Euclidean distance without the final square root, computed via
#' the identity \eqn{(x - y)^2 = x^2 - 2xy + y^2}.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .squared_euclidean(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Bray-Curtis distance backend (.braycurtis)
#'
#' @name .braycurtis
#' @description
#' Computes the Bray-Curtis dissimilarity between rows of `Ar` and `Br`.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}:
#' \deqn{d(x, y) = \frac{\sum_{i=1}^{k}|x_i - y_i|}{\sum_{i=1}^{k}|x_i + y_i|}}
#' A zero denominator yields a distance of `0`.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .braycurtis(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Hellinger distance backend (.hellinger)
#'
#' @name .hellinger
#' @description
#' Computes the Hellinger distance between rows of `Ar` and `Br`. Inputs are
#' assumed to be non-negative.
#' @details
#' For non-negative observations \eqn{x, y \in \mathbb{R}^k}:
#' \deqn{d(x, y) = \frac{1}{\sqrt{2}}
#' \sqrt{\sum_{i=1}^{k}(\sqrt{x_i} - \sqrt{y_i})^2}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .hellinger(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal chi-squared distance backend (.chisquared)
#'
#' @name .chisquared
#' @description
#' Computes the symmetric chi-squared distance between rows of `Ar` and `Br`.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}:
#' \deqn{d(x, y) = \frac{1}{2}\sum_{i=1}^{k}\frac{(x_i - y_i)^2}{x_i + y_i}}
#' Terms with a zero denominator are skipped.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .chisquared(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Jensen-Shannon distance backend (.jensenshannon)
#'
#' @name .jensenshannon
#' @description
#' Computes the Jensen-Shannon distance between rows of `Ar` and `Br`. Each row
#' is normalized to sum to one before the computation.
#' @details
#' Let \eqn{m = (x + y)/2}. The Jensen-Shannon distance is the square root of the
#' Jensen-Shannon divergence (natural logarithm):
#' \deqn{d(x, y) = \sqrt{\tfrac{1}{2}\sum_i x_i\log\frac{x_i}{m_i}
#' + \tfrac{1}{2}\sum_i y_i\log\frac{y_i}{m_i}}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .jensenshannon(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Haversine distance backend (.haversine)
#'
#' @name .haversine
#' @description
#' Computes the great-circle (Haversine) distance in kilometres between rows of
#' `Ar` and `Br`. Both matrices must have exactly two columns holding latitude
#' and longitude in degrees.
#' @details
#' With latitudes/longitudes in radians and Earth radius \eqn{R = 6371} km:
#' \deqn{a = \sin^2\!\big(\tfrac{\Delta\phi}{2}\big)
#' + \cos\phi_1\cos\phi_2\sin^2\!\big(\tfrac{\Delta\lambda}{2}\big),\quad
#' d = 2R\,\mathrm{asin}(\sqrt{a})}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @param Ar Numeric matrix of size `m x 2` with `(latitude, longitude)` rows.
#' @param Br Numeric matrix of size `n x 2` with `(latitude, longitude)` rows.
#'
#' @return Numeric `m x n` matrix with pairwise distances in kilometres.
#' @usage .haversine(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal standardized Euclidean distance backend (.standardized_euclidean)
#'
#' @name .standardized_euclidean
#' @description
#' Computes the standardized Euclidean distance between rows of `Ar` and `Br`,
#' scaling each feature by its sample variance estimated from `Ar`.
#' @details
#' Let \eqn{s_c^2} be the sample variance of column `c` of `Ar`. Then:
#' \deqn{d(x, y) = \sqrt{\sum_{c=1}^{k}\frac{(x_c - y_c)^2}{s_c^2}}}
#' Features with zero variance are dropped.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .standardized_euclidean(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Spearman correlation distance backend (.spearman)
#'
#' @name .spearman
#' @description
#' Computes the Spearman correlation distance between rows of `Ar` and `Br`.
#' @details
#' Each row is replaced by its within-row average ranks, and the Pearson
#' correlation \eqn{\rho_s} is computed on those ranks:
#' \deqn{d(x, y) = 1 - \rho_s(x, y)}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .spearman(Ar, Br)
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

#' Internal Hamming distance backend (.hamming)
#'
#' @name .hamming
#' @description
#' Computes the Hamming distance between rows of `Ar` and `Br`, i.e. the
#' number of positions in which the two rows differ.
#' @details
#' For observations \eqn{x, y \in \mathbb{R}^k}:
#' \deqn{d(x, y) = \sum_{i=1}^{k}\mathbf{1}(x_i \neq y_i)}
#' Intended for binary or integer-coded categorical data.
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
#' Computes the Jaccard distance between rows of `Ar` and `Br`, treating
#' non-zero entries as set membership.
#' @details
#' With the intersection and union computed over the non-zero patterns of
#' \eqn{x} and \eqn{y}:
#' \deqn{d(x, y) = 1 - \frac{|x \cap y|}{|x \cup y|}}
#' An empty union yields a distance of `0`.
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .jaccard(Ar, Br)
#' @keywords internal
#' @noRd
NULL

#' Internal Gower distance backend (.gower)
#'
#' @name .gower
#' @description
#' Computes the Gower distance between rows of `Ar` and `Br`, scaling each
#' feature by its range estimated from `Ar`.
#' @details
#' Let \eqn{R_c} be the range (max - min) of column `c` of `Ar` (columns with
#' zero range use \eqn{R_c = 1}). Then:
#' \deqn{d(x, y) = \frac{1}{k}\sum_{c=1}^{k}\frac{|x_c - y_c|}{R_c}}
#' This function is an internal backend implemented in C++ and exposed via Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Numeric `m x n` matrix with pairwise distances.
#' @usage .gower(Ar, Br)
#' @keywords internal
#' @noRd
NULL
