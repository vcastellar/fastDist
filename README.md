# fastDist

Fast computation of pairwise distances between the rows of two matrices.

`fastDist` provides a single entry point, `fdist()`, that dispatches to C++
backends built on [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)
and parallelised with OpenMP where available. Unlike `stats::dist()` and
packages built around it, `fdist()` computes **cross distances** directly:
given `A` (`m x k`) and `B` (`n x k`), it returns the `m x n` matrix of
distances between every row of `A` and every row of `B`, with no stacking or
slicing tricks. When `B` is omitted, only the upper triangle of the symmetric
result is computed.

## Installation

```r
# development version from GitHub
# install.packages("remotes")
remotes::install_github("vcastellar/fastDist")
```

## Usage

```r
library(fastDist)

set.seed(1)
A <- matrix(rnorm(1000 * 50), 1000, 50)
B <- matrix(rnorm( 500 * 50),  500, 50)

D <- fdist(A, B, method = "euclidean")   # 1000 x 500 cross distances
S <- fdist(A, method = "cosine")         # 1000 x 1000 symmetric matrix

# list all available methods
fdistregistry$get_entry_names()
```

## Supported distance methods

| method                    | intended data                       |
|---------------------------|-------------------------------------|
| `euclidean`               | general numeric                     |
| `manhattan`               | general numeric                     |
| `minkowski` (order `p`)   | general numeric                     |
| `supremum` (Chebyshev)    | general numeric                     |
| `squared_euclidean`       | general numeric (k-means, etc.)     |
| `standardized_euclidean`  | numeric, mixed units                |
| `correlation`             | profiles / expression data          |
| `cosine`                  | embeddings, text vectors            |
| `spearman`                | rank-based profiles                 |
| `canberra`                | non-negative, sparse counts         |
| `bray_curtis`             | abundances / compositional          |
| `hellinger`               | distributions                       |
| `chi_squared`             | histograms                          |
| `jensen_shannon`          | probability distributions           |
| `hamming`                 | binary / categorical codes          |
| `jaccard`                 | binary / sets                       |
| `gower`                   | mixed-scale numeric                 |
| `haversine`               | geographic (lat, lon) in degrees    |
| `mahalanobis`             | numeric, covariance-aware           |

See `?fdist` for the exact definition of each metric, and the vignette
(`vignette("fastDist")`) for a guided tour.

## Benchmarks

Two reproducible benchmark scripts live in the `benchmark/` directory of the
repository (they are not shipped with the installed package):

* `benchmark/benchmark_comparison.R` — compares `fdist()` against
  [`parallelDist`](https://cran.r-project.org/package=parallelDist),
  [`proxy`](https://cran.r-project.org/package=proxy) and
  [`rdist`](https://cran.r-project.org/package=rdist) for every metric the
  packages share (Euclidean, Manhattan, Minkowski, Canberra, supremum,
  Hamming, Jaccard). Results are validated against `fastDist` before timing,
  and a speed-up summary is printed at the end:

  ```r
  source("benchmark/benchmark_comparison.R")
  results <- run_benchmark(a_rows = 1000,
                           b_rows = c(1000, 5000, 10000),
                           n_features = 1000)
  speedup_table(results)
  plot_benchmark(results)   # requires ggplot2
  ```

* `benchmark/benchmark_parallelDist.R` — the original head-to-head
  comparison with `parallelDist` only.

As a reference, mean time in seconds to compute the cross distances between
`A` (`1000 x 1000`) and `B` (growing number of rows), measured with
`microbenchmark` on an Intel i5 10400. `parallelDist` computes square
matrices only, so its timings include the `rbind()` + slicing step needed to
obtain cross distances:

| expr          | method     | mean       | b_rows | a_rows | n_features |
|---------------|------------|------------|--------|--------|------------|
| fastDist      | euclidean  | 3.076106   | 1000   | 1000   | 1000       |
| parallelDist  | euclidean  | 5.131321   | 1000   | 1000   | 1000       |
| fastDist      | euclidean  | 17.333642  | 5000   | 1000   | 1000       |
| parallelDist  | euclidean  | 25.182693  | 5000   | 1000   | 1000       |
| fastDist      | euclidean  | 31.664112  | 10000  | 1000   | 1000       |
| parallelDist  | euclidean  | 51.903371  | 10000  | 1000   | 1000       |
| fastDist      | euclidean  | 79.663852  | 20000  | 1000   | 1000       |
| parallelDist  | euclidean  | 99.474717  | 20000  | 1000   | 1000       |
| fastDist      | manhattan  | 2.297173   | 1000   | 1000   | 1000       |
| parallelDist  | manhattan  | 4.679487   | 1000   | 1000   | 1000       |
| fastDist      | manhattan  | 17.064038  | 5000   | 1000   | 1000       |
| parallelDist  | manhattan  | 24.077118  | 5000   | 1000   | 1000       |
| fastDist      | manhattan  | 38.464093  | 10000  | 1000   | 1000       |
| parallelDist  | manhattan  | 54.814400  | 10000  | 1000   | 1000       |
| fastDist      | manhattan  | 89.921664  | 20000  | 1000   | 1000       |
| parallelDist  | manhattan  | 106.925199 | 20000  | 1000   | 1000       |
| fastDist      | minkowski  | 4.309440   | 1000   | 1000   | 1000       |
| parallelDist  | minkowski  | 9.192963   | 1000   | 1000   | 1000       |
| fastDist      | minkowski  | 37.301864  | 5000   | 1000   | 1000       |
| parallelDist  | minkowski  | 45.693993  | 5000   | 1000   | 1000       |
| fastDist      | minkowski  | 64.177987  | 10000  | 1000   | 1000       |
| parallelDist  | minkowski  | 94.068592  | 10000  | 1000   | 1000       |
| fastDist      | minkowski  | 139.660941 | 20000  | 1000   | 1000       |
| parallelDist  | minkowski  | 185.130829 | 20000  | 1000   | 1000       |

Timings depend on hardware, BLAS library and number of OpenMP threads; run
the scripts on your machine for meaningful numbers.

## License

GPL-3
