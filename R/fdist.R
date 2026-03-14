



distf <- function(A, B = NULL, method, ncores = NULL, p = NULL) {
  library(parallel)

  A <- as.matrix(A)
  if (is.null(B)) {
    B <- A
  } else {
    B <- as.matrix(B)
  }

  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregestry"))
  }

  numCores <- if (is.null(ncores)) {
    max(1L, min(detectCores(logical = FALSE), nrow(A)))
  } else {
    max(1L, min(as.integer(ncores), nrow(A)))
  }

  if (numCores == 1L) {
    return(fdist(A, B = B, method = method, p = p))
  }

  q <- nrow(A) %/% numCores
  r <- nrow(A) %% numCores

  groups <- c(rep(seq_len(numCores), each = q), rep(numCores, each = r))
  chunks <- lapply(split(A, groups), matrix, ncol = ncol(A))

  cl <- makeCluster(numCores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, library(fastDist))
  clusterExport(cl, varlist = c("B", "method", "p"), envir = environment())

  pieces <- parLapply(cl, chunks, fdist, B = B, method = method, p = p)
  do.call(rbind, pieces)
}

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
