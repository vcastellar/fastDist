# fastDist
calculo rápido de distancias entre filas de dos matrices

# benchmark

| expr |     min |      lq |    mean |  median |      uq |    max |  neval |
|---       |---  |     --- |   ---    |---      |---    |---   |---    |
|fastDis::euclidean|   79.92418|   85.66241|   93.03568|   90.13945|   96.71315|  179.1589|   100|
|prox::euclidean|  196.35367|  200.07790|  212.72639|  204.49260|  215.31226|  325.3437|   100|
|fastDist::manhattan|  132.54304|  138.07719|  143.40045|  141.14112|  145.63924|  198.4304|   100|
|proxy::manhattan|  208.29112|  212.70243|  222.91464|  215.76621|  228.15421|  282.6630|   100|
|fastDist::minkowsky| 2428.73658| 2439.57058| 2534.84836| 2473.57942| 2574.30003| 3159.3491|   100|
|proxy::minkowsky| 2885.72169| 2916.22851| 3082.57815| 2965.26809| 3151.75746| 4602.3627|   100|
