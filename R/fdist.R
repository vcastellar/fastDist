
#' Distancias por bloques en paralelo
#'
#' Calcula una matriz de distancias entre las filas de `A` y `B` repartiendo
#' las filas de `A` por bloques y ejecutando cada bloque en paralelo con
#' `parallel::parLapply()`. Internamente delega el cálculo en [fdist()].
#'
#' @param A Matriz o estructura coercible a matriz con las observaciones de
#'   origen (filas).
#' @param B Matriz o estructura coercible a matriz con las observaciones de
#'   destino (filas). Si es `NULL`, se usa `A`.
#' @param method Nombre del método de distancia registrado en `fdistregistry`.
#' @param ncores Número de núcleos para paralelizar. Si es `NULL`, se usa el
#'   valor configurado internamente por la función.
#' @param p Parámetro adicional para métodos que lo requieren (por ejemplo,
#'   Minkowski).
#'
#' @return Una matriz numérica con las distancias calculadas entre `A` y `B`.
#' @export
distf <- function(A, B = NULL, method, ncores = NULL, p = NULL) {
  library(parallel)

  A <- as.matrix(A)
  if (is.null(B)) {
    B <- A
  } else {
    B <- as.matrix(B)
  }
  
  nrowmax <- max(nrow(A), nrow(B))
  
  if (is.null(ncores)) {
    if (nrowmax > 100){
      numCores <- 1
      
    } else {
      numCores <- 1
    }
    
  }

  cl <- makeCluster(numCores, type = "PSOCK")
  clusterEvalQ(cl, library(fastDist))
  clusterExport(cl, varlist = c("B"))
  
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregestry"))
  }
  
  q <- nrow(A) %/% numCores
  r <- nrow(A) %% numCores
  
  my_list1 <- lapply(split(A,                    # Split matrix into list
                    c(rep(1:numCores, each = q ), 
                      rep(numCores,   each = r ))
                    ), matrix, ncol = ncol(A))
  
  
  
  as.matrix(do.call(rbind.data.frame, parLapply(cl, my_list1, fdist, method = method, B = B)))
}

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
