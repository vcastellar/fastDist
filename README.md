# fastDist
Cálculo rápido de distancias entre filas de dos matrices

# benchmark
Tiempo (segundos) de cálculo de distancias entre filas de una matriz de dimensiones 10000x100. Comparación de diferentes métodos entre el paquete fastDist y el proxy. Para el test, se usa el paquete microbenchmark. Ejecutado en un procesador pentium i5 10400

|               expr|        min|         lq|       mean|     median|        uq|       max| neval|
|---                |---        |---        |---        |---        |---       |---       |---   |
| fastDist_euclidean|   7.084282|   7.275046|   7.765501|   7.479555|  8.105057|  8.957104|    10|
|    proxy_euclidean|  44.228762|  45.879956|  47.817873|  46.242188| 50.503357| 53.047801|    10|
| fastDist_manhattan|  32.785770|  33.027793|  34.126146|  33.745409| 34.936131| 37.402658|    10|
|    proxy_manhattan|  45.070314|  46.908811|  48.802322|  47.368676| 52.432559| 53.713074|    10|
| fastDist_minkowsky| 187.029804| 188.430040| 191.739452| 190.020474|192.161082|202.709465|    10|
|    proxy_minkowsky| 427.731499| 429.940128| 437.436294| 432.872657|448.665750|454.400205|    10|
