\name{manhattan}
\alias{manhattan}

\title{Función para el cálculo de distancias de Manhattan entre las filas de dos matrices}
\description{
  Función para el cálculo de distancias de Manhattan entre las filas de dos matrices
}
\usage{
fastPDist2(A, B)
}
\arguments{
  \item{A}{matriz numérica de dimensiones \code{(m, k)} }
  \item{B}{matriz numérica de dimensiones \code{(n, k)} } 
}
\value{
  Una matriz de dimensiones \code{(m, n)} donde el valor de \code{(i, j)} es la distancia de Manhattan  entre la fila   \code{i} de la matriz \code{A} y la fila \code{j} de la matriz \code{B} 
}
\details{
  La distancia de Manhattan  entre dos vectores $x$ e $y$ se define como
  \deqn{||x-y||_1 = (\sum_{i=1}^k |x_i - y_i| }
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  A <- matrix(rnorm(5000), 500, 10)
  B <- matrix(rnorm(1000), 100, 10)
  manhattan(A, B)
}
\author{Vicente Castellar Sebastiá}

