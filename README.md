<<<<<<< HEAD
# fastDist
cálculo rápido de distancias entre filas de dos matrices

# benchmark
Distancias entre filas de una matriz de dimensiones 1000x100

|               expr|        min|         lq|       mean|     median|        uq|       max| neval|
|---                |---        |---        |---        |---        |---       |---       |---   |
| fastDist_euclidean|   57.83101|   63.26887|   69.61439|  6 7.79110|   73.5964|  107.1624|   100|
|    proxy_euclidean|  136.57597|  140.38191|  149.12514|  143.04433|  150.2885|  218.0038|   100|
| fastDist_manhattan|   90.16858|   92.71555|  100.71742|   95.24281|  104.2436|  167.1578|   100|
|    proxy_manhattan|  144.59291|  148.41961|  159.20060|  152.34862|  161.6213|  227.4017|   100|
| fastDist_minkowsky| 2279.70868| 2355.83950| 2460.24340| 2419.81075| 2496.2675| 3383.1511|   100|
|    proxy_minkowsky| 2741.66330| 2866.91688| 3002.01956| 2936.23840| 3074.8196| 4092.9548|   100|
