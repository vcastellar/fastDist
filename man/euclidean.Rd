\name{euclidean}
\alias{euclidean}

\title{Euclidean distances between the rows of two matrices}
\description{
  Function for calculating Euclidean distances between the rows of two matrices
}
\usage{
fastPDist2(A, B)
}
\arguments{
  \item{A}{numerical matrix with dimesions \code{(m, k)} }
  \item{B}{numerical matrix with dimesions \code{(n, k)} } 
}
\value{
  A numeric matrix with dimensions \code{(m, n)} where \code{(i, j)} is the euclidean distance between the row   \code{i} of mnatrix \code{A} and the row \code{j} of matrix \code{B} 
}
\details{
  The euclidean distance between two vectors $x$ e $y$ is
  \deqn{||x-y||_2 = (\sum_{i=1}^k (x_i - y_i)^2)^{1/2} }
  We use the property \eqn{(x - y)^2 = x^2 -2xy + y^2} for increase the speed
  
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

