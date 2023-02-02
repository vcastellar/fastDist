# fastDist
Cálculo rápido de distancias entre filas de dos matrices

# benchmark
Tiempo (segundos) de cálculo de distancias entre filas de una matriz de dimensiones 10000x100. Comparación de diferentes métodos entre el paquete fastDist y el proxy. Para el test, se usa el paquete microbenchmark. Ejecutado en un procesador pentium i5 10400

|               expr|        min|         lq|       mean|     median|        uq|       max| neval|
|:---                |       ---:|       ---:|       ---:|       ---:|      ---:|      ---:|  ---:|
| fastDist_euclidean|   7.084282|   7.275046|   7.765501|   7.479555|  8.105057|  8.957104|    10|
|    proxy_euclidean|  44.228762|  45.879956|  47.817873|  46.242188| 50.503357| 53.047801|    10|
| fastDist_manhattan|  32.785770|  33.027793|  34.126146|  33.745409| 34.936131| 37.402658|    10|
|    proxy_manhattan|  45.070314|  46.908811|  48.802322|  47.368676| 52.432559| 53.713074|    10|
| fastDist_minkowsky| 187.029804| 188.430040| 191.739452| 190.020474|192.161082|202.709465|    10|
|    proxy_minkowsky| 427.731499| 429.940128| 437.436294| 432.872657|448.665750|454.400205|    10|

| Unit: | milliseconds         |          |          |          |          |          |          |       |
|-------|----------------------|----------|----------|----------|----------|----------|----------|-------|
|       | expr                 | min      | lq       | mean     | median   | uq       | max      | neval |
|       | fastDist_euclidean   | 63.3     | 64.6     | 69.0     | 66.7     | 72.3     | 81.9     | 30    |
|       | proxy_euclidean      | 389.2    | 396.3    | 405.9    | 404.6    | 414.1    | 441.3    | 30    |
|       | fastDist_manhattan   | 73.1     | 74.8     | 78.4     | 76.3     | 81.8     | 90.5     | 30    |
|       | proxy_manhattan      | 395.6    | 403.5    | 411.7    | 408.3    | 415.8    | 452.7    | 30    |
|       | fastDist_minkowsky   | 1714.7   | 1789.6   | 1839.4   | 1810.5   | 1860.5   | 2176.6   | 30    |
|       | proxy_minkowsky      | 4060.9   | 4144.8   | 4250.8   | 4211.8   | 4313.8   | 4717.8   | 30    |
|       | fastDist_canberra    | 114.9    | 117.3    | 128.5    | 121.6    | 129.8    | 188.0    | 30    |
|       | proxy_canberra       | 514.6    | 522.2    | 539.7    | 534.6    | 541.9    | 696.5    | 30    |
|       | fastDist_supremum    | 102.9    | 104.3    | 111.2    | 107.5    | 114.7    | 170.2    | 30    |
|       | proxy_supremum       | 395.2    | 401.1    | 413.0    | 408.1    | 419.7    | 483.7    | 30    |
|       | fastDist_mahalanobis | 2040.2   | 2078.8   | 2108.2   | 2100.7   | 2128.0   | 2244.6   | 30    |
|       | proxy_mahalanobis    | 107772.6 | 109274.8 | 111446.1 | 110321.4 | 112152.8 | 122567.6 | 30    |
