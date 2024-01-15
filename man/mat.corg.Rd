\name{mat.corg}
\alias{mat.corg}
\title{Forms a general correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} for which  
             all correlations potentially differ.}
\usage{mat.corg(rhos, order, byrow = FALSE)}
\arguments{
 \item{rhos}{A \code{\link{numeric}} containing the \code{p(p -1)/2} correlation values ordered 
             either by columns (if \code{byrow} is FALSE) or by rows (if \code{byrow} is TRUE).}
 \item{order}{The order of the correlation \code{\link{matrix}} to be formed.}
 \item{byrow}{A \code{\link{logical}}. If \code{FALSE} the lower-traingle of the matrix is filled 
             by columns, otherwise the the ower triangle is filled by rows.}
}
\value{A correlation \code{\link{matrix}}.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}},  \code{\link{mat.banded}}, \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.corg(rhos = c(0.4, 0.2, 0.1), order = 3)
}
\keyword{array}
