#' Fast Distance Computation Between Observations
#'
#' Computes pairwise distances between rows of `A` and `B` using a distance
#' backend registered in `fdistregistry`. When `method = "mahalanobis"`, the
#' function dispatches to the method-specific signature that only requires `A`.
#'
#' @param A A matrix (or object coercible to a matrix) containing source
#'   observations in rows.
#' @param B A matrix (or object coercible to a matrix) containing target
#'   observations in rows. If `NULL`, `A` is used.
#' @param method The name of a distance method registered in `fdistregistry`.
#' @param p Optional extra parameter used by methods that require it (for
#'   example, Minkowski).
#'
#' @return A numeric matrix containing computed distances.
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
