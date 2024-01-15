\name{mat.cor}
\alias{mat.cor}
\title{Forms a correlation matrix in which all correlations have the same value.}
\description{Form the correlation \code{\link{matrix}} of order \code{order} in which  
             all correlations have the same value.}
\usage{mat.cor(rho, order)}
\arguments{
 \item{rho}{A \code{\link{numeric}} containing the single correlation value.}
 \item{order}{The order of the correlation \code{\link{matrix}} to be formed.}
}
\value{A correlation \code{\link{matrix}}.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.corg}}, \code{\link{mat.banded}}, \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.cor(rho = 0.4, order = 3)
}
\keyword{array}
