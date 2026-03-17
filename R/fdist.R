#' Distancias rápidas entre observaciones
#'
#' Calcula distancias entre filas de `A` y `B` usando una implementación
#' registrada en `fdistregistry`. Si `method = "mahalanobis"`, el cálculo se
#' realiza con la firma específica de ese método.
#'
#' @param A Matriz o estructura coercible a matriz con las observaciones de
#'   origen (filas).
#' @param B Matriz o estructura coercible a matriz con las observaciones de
#'   destino (filas). Si es `NULL`, se usa `A`.
#' @param method Nombre del método de distancia registrado en `fdistregistry`.
#' @param p Parámetro adicional para métodos que lo requieren (por ejemplo,
#'   Minkowski).
#'
#' @return Una matriz numérica con las distancias calculadas.
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
