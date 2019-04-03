\name{mat.gau}
\alias{mat.gau}
\title{Forms an exponential correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order equal to the length of 
             \code{coordinates}. The \code{\link{matrix}} has diagonal 
             elements equal to one and the off-diagonal element in the ith row 
             and jth column equal to \eqn{\rho^k}{rho^(k)} where 
             \eqn{k = (coordinate[i]- coordinate[j])^2}{k = (coordinate[i]- coordinate[j])^2}.}
\usage{mat.gau(rho, coordinates)}
\arguments{
 \item{rho}{The correlation for points a distance of one apart.}
 \item{coordinates}{The coordinates of points whose correlation \code{\link{matrix}} 
                    is to be formed.}
}
\value{A correlation \code{\link{matrix}} whose elements depend on the power of the 
       absolute distance apart.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.banded}}, \code{\link{mat.ar1}}, 
         \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, \code{\link{mat.sar2}}, 
         \code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}},
         \code{\link{mat.exp}} 
}
\examples{
    corr <- mat.gau(coordinates=c(3:6, 9:12, 15:18), rho=0.1)
}
\keyword{array}
