# fastDist
Fast computation of distances between the rows of two matrices

# benchmark
Computation time (seconds) of distances between the rows of a 10000x100 matrix. Comparison of different methods between the fastDist package and proxy. The test uses the microbenchmark package. Run on an Intel i5 10400 processor

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

## benchmark fastDist vs parallelDist

A reproducible benchmark was added to compare `fastDist` against `parallelDist`
when computing distances between the rows of `A` and `B` for the methods:

- Euclidean
- Manhattan
- Minkowski

The script is located at `inst/benchmarks/parallelDist_microbenchmark.R` and runs the case:

- `A`: `1000 x 1000`
- `B`: `1000`, `5000`, `10000` and `20000` rows

Suggested usage:

```r
install.packages(c("microbenchmark", "parallelDist"))
devtools::load_all(".")
source("inst/benchmarks/parallelDist_microbenchmark.R")
```

Internally, the comparison with `parallelDist` is solved in blocks over the rows
of `B`, because `parallelDist::parDist()` computes square distance matrices. This
way the cross submatrix `A x B` can be extracted in order to study how the time
scales as `B` grows.
