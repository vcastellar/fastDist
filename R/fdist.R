fdist <- function(A,B, p = NULL, method) {
  if (is.na(fdistregistry$get_entry(method)$p)) {
    fdistregistry$get_entry(method)$fun(A, B)
  } else {
    fdistregistry$get_entry(method)$fun(A, B, p)
  }
}