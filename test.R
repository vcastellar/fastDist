library(proxy)
library(fastDist)
library(microbenchmark)
library(ggplot2)

A <- matrix(rnorm(50), 5, 10)
fastDist::euclidean(A, A)
proxy::dist(A, A, method = "Euclidean")

fastDist::manhattan(A, A)
proxy::dist(A, A, method = "Manhattan")

fastDist::minkowsky(A, A, p = 5)
proxy::dist(A, A, method = "Minkowski", p = 5)


# benchmark
#------------------------------------------------------------------------------
B <- matrix(rnorm(100000), 1000, 100)

rows <- seq(10, 100, 10) * 1e2
cols <- 100

res <- microbenchmark(
  fastDist_euclidean = euclidean(B, B),
  proxy_euclidean    = dist(B, B, method = "Euclidean"),
  fastDist_manhattan = manhattan(B, B),
  proxy_manhattan    = dist(B, B, method = "Manhattan"),
  fastDist_minkowsky = minkowsky(B, B, p = 5),
  proxy_minkowsky    = dist(B, B, method = "Minkowski", p = 5)
)
autoplot(res)

resultados = data.frame(method = character(),
                        package = character(),
                        nrow = integer(),
                        ncol = integer(),
                        t = numeric())

for (row in rows) {
  B <- matrix(rnorm(row * cols), row, cols)
  
  res <- system.time({
    fastDist::euclidean(B, B)
  })
  fila <- data.frame(method = "Euclidean", package = "fastDist", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  
  res <- system.time({
    proxy::dist(B, B, method = "Euclidean")
  })
  fila <- data.frame(method = "Euclidean", package = "proxy", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  
  res <- system.time({
    fastDist::manhattan(B, B)
  })
  fila <- data.frame(method = "Manhattan", package = "fastDist", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  
  
  res <- system.time({
    proxy::dist(B, B, method = "Manhattan")
  })
  fila <- data.frame(method = "Manhattan", package = "proxy", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  
  res <- system.time({
    fastDist::minkowsky(B, B, 5)
  })
  fila <- data.frame(method = "Minkowsky", package = "fastDist", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  
  res <- system.time({
    proxy::dist(B, B, method = "Minkowski", p = 5)
  })
  fila <- data.frame(method = "Minkowsky", package = "proxy", nrow = row, ncol = 100, t = res[3])
  resultados <- rbind(resultados, fila)
  print(resultados)
}

library(ggplot2)
png(filename = "test.png", width = 1200, height = 500)
ggplot(resultados, aes(nrow, t, color = package)) + 
  geom_point() + geom_line() +
  facet_wrap(. ~ method, scales = "free")
dev.off()
