\name{fastPDist2}
\alias{fastPDist2}

\title{Función para el cálculo de distancias euclídeas entre las filas de dos matrices}
\description{
  Función para el cálculo de distancias euclídeas entre las filas de dos matrices
}
\usage{
fastPDist2(A, B)
}
\arguments{
  \item{A}{matriz numérica de dimensiones \code{(m, k)} }
  \item{B}{matriz numérica de dimensiones \code{(n, k)} } 
}
\value{
  Una matriz de dimensiones \code{(m, n)} donde el valor de \code{(i, j)} es la distancia euclídea entre la fila   \code{i} de la matriz \code{A} y la fila \code{j} de la matriz \code{B} 
}
\details{
  La distancia euclidea entre dos vectores $x$ e $y$ se define como
  \deqn{||x-y||_2 = (\sum_{i=1}^k (x_i - y_i)^2)^{1/2} }
  se usa la propiedad \eqn{(x - y)^2 = x^2 -2xy + y^2} para acelerar el cálculo
  
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  A <- matrix(rnorm(5000), 500, 10)
  B <- matrix(rnorm(1000), 100, 10)
  fastPdist2(A, B)
}
\author{Vicente Castellar Sebastiá}

