\name{mat.ma1}
\alias{mat.ma1}
\title{Forms an ma1 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ma1 pattern. The \code{\link{matrix}} is banded and 
             has diagonal elements equal to one and subdiagonal element equal to \cr 
             \code{-MAparameter / (1 + MAparameter*MAparameter)}.}
\usage{mat.ma1(MAparameter, order)}
\arguments{
 \item{MAparameter}{The moving average parameter, being the weight applied to the lag 1 random 
       pertubation.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an ma1 pattern.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, 
         \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
         \code{\link{mat.banded}}, \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, 
         \code{\link{mat.sar2}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.ma1(MAparameter=0.4, order=4)
}
\keyword{array}
