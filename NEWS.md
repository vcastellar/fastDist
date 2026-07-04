# fastDist 1.1.0

## New features

* New distance methods: `"hamming"` (categorical/binary data), `"jaccard"`
  (binary/set data) and `"gower"` (range-scaled mixed data), bringing the
  total to 19 metrics.
* `fdist(method = "minkowski")` no longer requires `p`: it defaults to the
  value stored in the registry (`p = 2`).
* New vignette *"Computing pairwise distances with fastDist"* covering all
  metrics, grouped by the kind of data they are intended for.
* New reproducible benchmark script `benchmark/benchmark_comparison.R`
  (package sources only) comparing `fdist()` with `parallelDist::parDist()`,
  `proxy::dist()` and `rdist::cdist()`, with a correctness check before
  timing.

## Bug fixes and internal changes

* `registry` is now declared in `Imports` and used via
  `registry::registry()` instead of `library(registry)` inside the package
  code; `LinkingTo` now only lists `Rcpp` and `RcppArmadillo`.
* `src/Makevars` and `src/Makevars.win` cleaned up: removed the unused
  `RcppParallel` linkage and the obsolete `CXX_STD = CXX11` setting; OpenMP
  is enabled through the standard `$(SHLIB_OPENMP_CXXFLAGS)` mechanism.
* `NAMESPACE` now exports `fdist()` and `fdistregistry` explicitly.
* `fdistregistry` is documented and `fdist()` gained runnable examples.

# fastDist 1.0

* First public version with the methods: Euclidean, Manhattan, Minkowski,
  correlation, cosine, Canberra, supremum, squared Euclidean, Bray-Curtis,
  Hellinger, chi-squared, Jensen-Shannon, Haversine, standardized Euclidean,
  Spearman and Mahalanobis.
