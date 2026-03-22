# fastDist
Cálculo rápido de distancias entre filas de dos matrices

# benchmark
Tiempo (segundos) de cálculo de distancias entre filas de una matriz de dimensiones 10000x100. Comparación de diferentes métodos entre el paquete fastDist y el proxy. Para el test, se usa el paquete microbenchmark. Ejecutado en un procesador pentium i5 10400


\begin{table}[ht]
\centering
\begin{tabular}{llrrrrr}
\hline
expr & method & mean & b\_rows & a\_rows & n\_features \\
\hline
fastDist      & euclidean & 3.076106  & 1000  & 1000 & 1000 \\
parallelDist  & euclidean & 5.131321  & 1000  & 1000 & 1000 \\
fastDist      & euclidean & 17.333642 & 5000  & 1000 & 1000 \\
parallelDist  & euclidean & 25.182693 & 5000  & 1000 & 1000 \\
fastDist      & euclidean & 31.664112 & 10000 & 1000 & 1000 \\
parallelDist  & euclidean & 51.903371 & 10000 & 1000 & 1000 \\
fastDist      & euclidean & 79.663852 & 20000 & 1000 & 1000 \\
parallelDist  & euclidean & 99.474717 & 20000 & 1000 & 1000 \\
fastDist      & manhattan & 2.297173  & 1000  & 1000 & 1000 \\
parallelDist  & manhattan & 4.679487  & 1000  & 1000 & 1000 \\
fastDist      & manhattan & 17.064038 & 5000  & 1000 & 1000 \\
parallelDist  & manhattan & 24.077118 & 5000  & 1000 & 1000 \\
fastDist      & manhattan & 38.464093 & 10000 & 1000 & 1000 \\
parallelDist  & manhattan & 54.814400 & 10000 & 1000 & 1000 \\
fastDist      & manhattan & 89.921664 & 20000 & 1000 & 1000 \\
parallelDist  & manhattan & 106.925199 & 20000 & 1000 & 1000 \\
fastDist      & minkowski & 4.309440  & 1000  & 1000 & 1000 \\
parallelDist  & minkowski & 9.192963  & 1000  & 1000 & 1000 \\
fastDist      & minkowski & 37.301864 & 5000  & 1000 & 1000 \\
parallelDist  & minkowski & 45.693993 & 5000  & 1000 & 1000 \\
fastDist      & minkowski & 64.177987 & 10000 & 1000 & 1000 \\
parallelDist  & minkowski & 94.068592 & 10000 & 1000 & 1000 \\
fastDist      & minkowski & 139.660941 & 20000 & 1000 & 1000 \\
parallelDist  & minkowski & 185.130829 & 20000 & 1000 & 1000 \\
\hline
\end{tabular}
\caption{Benchmark de tiempos medios para distintas métricas y tamaños}
\end{table}

## benchmark fastDist vs parallelDist

Se agregó un benchmark reproducible para comparar `fastDist` contra `parallelDist`
en el cálculo de distancias entre las filas de `A` y `B` para los métodos:

- Euclidean
- Manhattan
- Minkowski

El script está en `inst/benchmarks/parallelDist_microbenchmark.R` y ejecuta el caso:

- `A`: `1000 x 1000`
- `B`: `1000`, `5000`, `10000` y `20000` filas

Uso sugerido:

```r
install.packages(c("microbenchmark", "parallelDist"))
devtools::load_all(".")
source("inst/benchmarks/parallelDist_microbenchmark.R")
```

Internamente, la comparación con `parallelDist` se resuelve por bloques sobre las filas
de `B`, porque `parallelDist::parDist()` calcula matrices de distancia cuadradas. De esta
forma se puede extraer la submatriz cruzada `A x B` y estudiar cómo escala el tiempo a
medida que crece `B`.
