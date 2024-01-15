\name{mat.ar1}
\alias{mat.ar1}
\title{Forms an ar1 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ar1 pattern. The \code{\link{matrix}} is banded and 
             has diagonal elements equal to one and the off-diagonal element in the ith row 
             and jth column equal to \eqn{\rho^k}{rho^(k)} where 
             \eqn{k = |i- j|}{k = abs(i- j)}.}
\usage{mat.ar1(rho, order)}
\arguments{
 \item{rho}{The correlation on the first off-diagonal.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an ar1 pattern.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, 
         \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
         \code{\link{mat.banded}}, \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, \code{\link{mat.sar2}}, 
         \code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.ar1(rho=0.4, order=4)
}
\keyword{array}
