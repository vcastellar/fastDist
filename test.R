library(proxy)
library(fastDist)

A <- matrix(rnorm(50), 5, 10)
B <- matrix(rnorm(500000), 5000, 100)

fastDist::euclidean(A, A)
proxy::dist(A, A, method = "Euclidean")

fastDist::manhattan(A, A)
proxy::dist(A, A, method = "Manhattan")

fastDist::minkowsky(A, A, p = 5)
proxy::dist(A, A, method = "Minkowski", p = 5)



system.time({
  fastDist::euclidean(B, B)
})
system.time({
  proxy::dist(B, B, method = "Euclidean")
})
system.time({
  fastDist::manhattan(B, B)
})
system.time({
  proxy::dist(B, B, method = "Manhattan")
})
system.time({
  fastDist::minkowsky(B, B, 5)
})
system.time({
  proxy::dist(B, B, method = "Minkowski", p = 5)
})

