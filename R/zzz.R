# Reinforce native symbol bindings for Rcpp wrappers at load time.
# This keeps generated R/RcppExports.R untouched while preventing
# NULL symbol-address errors in environments where bindings are stale.
.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  symbols <- c(
    "_fastDist_euclidean",
    "_fastDist_manhattan",
    "_fastDist_minkowski",
    "_fastDist_correlation",
    "_fastDist_cosine",
    "_fastDist_canberra",
    "_fastDist_supremum",
    "_fastDist_mahalanobis"
  )

  for (sym in symbols) {
    assign(sym, getNativeSymbolInfo(sym, PACKAGE = pkgname), envir = ns)
  }
}
