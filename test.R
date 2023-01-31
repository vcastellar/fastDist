library(proxy)
library(fastDist)
library(microbenchmark)
library(ggplot2)
library(parallelDist)
library(parallel)
A <- matrix(1:15, 5, 3, byrow = T)
B <- matrix(1:12, 4, 3)

(X <- fdist(A, method = "euclidean"))
(Y <- distf(A, method = "euclidean"))
proxy::dist(A,  method = "Euclidean")

fdist(A, method = "manhattan")
fdist(A, B, method = "manhattan")
distf(A, B, method = "manhattan")


proxy::dist(A, method = "Manhattan")
proxy::dist(A, B, method = "Manhattan")

distf(A, method = "manhattan")

fdist(A, p = 5, method = "minkowski")
proxy::dist(A,  method = "Minkowski", p = 5)

fdist(A, B, p = 5, method = "minkowski")
proxy::dist(A, B, method = "Minkowski", p = 5)


data(mtcars)
fdist(mtcars, method = "euclidean")
proxy::dist(mtcars, method = "euclidean")

fdist(1:10, method = "euclidean")
proxy::dist(1:10, method = "euclidean")

proxy::dist(A, B,  method = "Minkowski", p = 5, pairwise = TRUE)
# benchmark
#------------------------------------------------------------------------------
B <- matrix(rnorm(1000000), 10000, 100)

rows <- seq(100, 100, 20) * 1e2
cols <- 100

res <- microbenchmark(
  fastDist_euclidean = fdist(B, B, method = "euclidean"),
  proxy_euclidean    = dist(B, B, method = "Euclidean"),
  fastDist_manhattan = fdist(B, B, method = "manhattan"),
  proxy_manhattan    = dist(B, B, method = "Manhattan"),
  fastDist_minkowsky = fdist(B, B, p = 5, method = "minkowski"),
  proxy_minkowsky    = dist(B, B, method = "Minkowski", p = 5),
  times = 10
)
autoplot(res)


#------------------------------------------------------------------------------
B <- matrix(rnorm(1000000), 10000, 100)

rows <- seq(100, 100, 20) * 1e2
cols <- 100

res <- microbenchmark(
  fastDist_euclidean = fdist(B, B, method = "euclidean"),
  fastDist_manhattan = fdist(B, B, method = "manhattan"),
  
  proxy_manhattan    = proxy::dist(B, method = "Manhattan"),
  proxy_euclidean    = proxy::dist(B, method = "Euclidean"),

  times = 10
)
autoplot(res)


resultados = data.frame(method = character(),
                        package = character(),
                        nrow = integer(),
                        ncol = integer(),
                        t = numeric())

for (row in rows) {
  B <- matrix(rnorm(row * cols), row, cols)
  
  res <- microbenchmark(fdist(B, B, method = "euclidean"), times = 10, unit = "seconds")
  fila <- data.frame(method = "Euclidean", package = "fastDist", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  res <- microbenchmark(proxy::dist(B, B, method = "Euclidean") ,times = 10, unit = "seconds")
  fila <- data.frame(method = "Euclidean", package = "proxy", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  res <- microbenchmark(fdist(B, B, method = "manhattan"), times = 10, unit = "seconds")
  fila <- data.frame(method = "Manhattan", package = "fastDist", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  res <- microbenchmark(proxy::dist(B, B, method = "Manhattan") ,times = 10, unit = "seconds")
  fila <- data.frame(method = "Manhattan", package = "proxy", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  res <- microbenchmark(fdist(B, B, method = "minkowski", p = 5), times = 10, unit = "seconds")
  fila <- data.frame(method = "Minkowski", package = "fastDist", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  res <- microbenchmark(proxy::dist(B, B, method = "Minkowski", p = 5) ,times = 10, unit = "seconds")
  fila <- data.frame(method = "Minkowski", package = "proxy", nrow = row, ncol = 100, t = summary(res)$mean)
  resultados <- rbind(resultados, fila)
  
  print(resultados)
}

library(ggplot2)
png(filename = "test.png", width = 1200, height = 500)
ggplot(resultados, aes(nrow, t, color = package)) + 
  geom_point() + geom_line() +
  facet_wrap(. ~ method, scales = "free")
dev.off()



