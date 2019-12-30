\name{mat.dirsum}
\alias{mat.dirsum}
\title{Forms the direct sum of a list of matrices}
\description{The direct sum is the partitioned matrices whose diagonal submatrices are 
             the matrices from which the direct sum is to be formed and whose off-diagonal 
             submatrices are conformable matrices of zeroes. The resulting 
             \code{\link{matrix}} is \eqn{m \times n}{m x n}, where \eqn{m}{m} is the sum of
             the numbers of rows of the contributing matrices and \eqn{n}{n} is the sum of
             their numbers of columns.}
\usage{mat.dirsum(matrices)}
\arguments{
 \item{matrices}{A list, each of whose component is a \code{\link{matrix}}.}
}
\value{An \eqn{m \times n}{m x n} \code{\link{matrix}}.}
\author{Chris Brien}
\seealso{\code{\link{mat.dirprod}}, \code{matmult}}
\examples{
       m1 <- matrix(1:4, nrow=2)
       m2 <- matrix(11:16, nrow=3)
       m3 <- diag(1, nrow=2, ncol=2)
       dsum <- mat.dirsum(list(m1, m2, m3))
}
\keyword{array}
