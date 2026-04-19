#' Distancia euclidiana interna (.euclidean)
#'
#' @name .euclidean
#' @description
#' Calcula la distancia euclidiana entre cada par de filas de `Ar` y `Br`.
#' @details
#' Para dos observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia euclidiana
#' se define como:
#' \deqn{d(x, y) = \sqrt{\sum_{i=1}^{k}(x_i - y_i)^2}}
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @param Ar Matriz numérica de tamaño `m x k` con observaciones por fila.
#' @param Br Matriz numérica de tamaño `n x k` con observaciones por fila.
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .euclidean(Ar, Br)
#' @keywords internal
NULL

#' Distancia Manhattan interna (.manhattan)
#'
#' @name .manhattan
#' @description
#' Calcula la distancia Manhattan (norma \eqn{L_1}) entre filas de `Ar` y `Br`.
#' @details
#' Para dos observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia Manhattan
#' se define como:
#' \deqn{d(x, y) = \sum_{i=1}^{k}|x_i - y_i|}
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .manhattan(Ar, Br)
#' @keywords internal
NULL

#' Distancia Minkowski interna (.minkowski)
#'
#' @name .minkowski
#' @description
#' Calcula la distancia Minkowski entre filas de `Ar` y `Br` para un valor `p`.
#' @details
#' Para \eqn{p > 0} y observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia
#' Minkowski se define como:
#' \deqn{d(x, y) = \left(\sum_{i=1}^{k}|x_i - y_i|^p\right)^{1/p}}
#' Casos particulares: \eqn{p=1} (Manhattan) y \eqn{p=2} (Euclidiana).
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#' @param p Orden de la distancia Minkowski, con `p > 0`.
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .minkowski(Ar, Br, p)
#' @keywords internal
NULL

#' Distancia de correlación interna (.correlation)
#'
#' @name .correlation
#' @description
#' Calcula la distancia de correlación entre filas de `Ar` y `Br`.
#' @details
#' La distancia de correlación se define como uno menos la correlación de Pearson:
#' \deqn{d(x, y) = 1 - \mathrm{corr}(x, y)}
#' donde
#' \deqn{\mathrm{corr}(x, y) =
#' \frac{\sum_{i=1}^{k}(x_i - \bar{x})(y_i - \bar{y})}
#' {\sqrt{\sum_{i=1}^{k}(x_i - \bar{x})^2}\sqrt{\sum_{i=1}^{k}(y_i - \bar{y})^2}}}
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .correlation(Ar, Br)
#' @keywords internal
NULL

#' Distancia del coseno interna (.cosine)
#'
#' @name .cosine
#' @description
#' Calcula la distancia del coseno entre filas de `Ar` y `Br`.
#' @details
#' La distancia del coseno se define como:
#' \deqn{d(x, y) = 1 - \frac{x \cdot y}{\|x\|_2\|y\|_2}}
#' donde \eqn{x \cdot y} es el producto punto y \eqn{\|\cdot\|_2} es la norma
#' euclidiana.
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .cosine(Ar, Br)
#' @keywords internal
NULL

#' Distancia Canberra interna (.canberra)
#'
#' @name .canberra
#' @description
#' Calcula la distancia Canberra entre filas de `Ar` y `Br`.
#' @details
#' Para observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia Canberra se define
#' como:
#' \deqn{d(x, y) = \sum_{i=1}^{k}
#' \frac{|x_i - y_i|}{|x_i| + |y_i|}}
#' En la implementación, los términos con denominador cero no aportan al total.
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .canberra(Ar, Br)
#' @keywords internal
NULL

#' Distancia suprema interna (.supremum)
#'
#' @name .supremum
#' @description
#' Calcula la distancia suprema (Chebyshev, norma \eqn{L_\infty}) entre filas de
#' `Ar` y `Br`.
#' @details
#' Para observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia suprema se define
#' como:
#' \deqn{d(x, y) = \max_{i=1,\ldots,k}|x_i - y_i|}
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @inheritParams .euclidean
#'
#' @return Matriz numérica `m x n` con distancias por pares.
#' @usage .supremum(Ar, Br)
#' @keywords internal
NULL

#' Distancia Mahalanobis interna (.mahalanobis)
#'
#' @name .mahalanobis
#' @description
#' Calcula la matriz de distancias de Mahalanobis entre todas las filas de `Ar`.
#' @details
#' Sea \eqn{S} la matriz de covarianza muestral de `Ar`, y \eqn{S^{-1}} su inversa.
#' Para dos observaciones \eqn{x, y \in \mathbb{R}^k}, la distancia de Mahalanobis
#' se define como:
#' \deqn{d(x, y) = \sqrt{(x - y)^\top S^{-1}(x - y)}}
#' Esta función es un backend interno implementado en C++ y expuesto vía Rcpp.
#'
#' @param Ar Matriz numérica de tamaño `m x k` con observaciones por fila.
#'
#' @return Matriz numérica cuadrada `m x m` con distancias por pares.
#' @usage .mahalanobis(Ar)
#' @keywords internal
NULL
