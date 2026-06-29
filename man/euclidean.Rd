\name{euclidean}
\alias{euclidean}

\title{Euclidean distances between the rows of two matrices}
\description{
  Function for calculating Euclidean distances between the rows of two matrices
}
\arguments{
  \item{A}{numerical matrix with dimensions \code{(m, k)} }
  \item{B}{numerical matrix with dimensions \code{(n, k)} }
}
\value{
  A numeric matrix with dimensions \code{(m, n)} where \code{(i, j)} is the euclidean distance between the row \code{i} of matrix \code{A} and the row \code{j} of matrix \code{B}
}
\details{
  The euclidean distance between two vectors $x$ and $y$ is
  \deqn{||x-y||_2 = (\sum_{i=1}^k (x_i - y_i)^2)^{1/2} }
  We use the property \eqn{(x - y)^2 = x^2 -2xy + y^2} to increase the speed
  
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\author{Vicente Castellar Sebastiá}

