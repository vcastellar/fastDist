library(RcppParallel)
RcppParallel::setThreadOptions(numThreads = 4)

fdist <- function(A, B = NULL, method, p = NULL) {
  if (!method %in% fdistregistry$get_entry_names()) {
    stop(paste(method, "not found in fdistregestry"))
  }
  A <- as.matrix(A)
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
  
  return(result)
}