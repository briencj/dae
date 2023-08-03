\name{as.numfac}
\alias{as.numfac}
\title{Convert a factor to a numeric vector}
\description{Converts a \code{\link{factor}} to a numeric \code{vector} with approximately the 
             numeric values of its \code{levels}. Hence, the \code{levels} of the 
             \code{\link{factor}} must be numeric values, stored as characters. It uses the method 
             described in \code{\link{factor}}.  Use \code{\link{as.numeric}} to convert a 
             \code{\link{factor}} to a numeric \code{vector} with integers 1, 2, ... corresponding 
             to the positions in the list of levels. You can also use \code{\link{fac.recast}} to 
             recode the levels to numeric values. If a \code{numeric} is supplied, it is 
             left unchanged.}
\usage{as.numfac(factor)}
\arguments{
 \item{factor}{The \code{\link{factor}} to be converted.}
}
\value{A numeric \code{vector}. An \code{NA} will be stored for any value of the factor whose 
       level is not a number.}
\author{Chris Brien}
\seealso{\code{\link{as.numeric}}, \code{\link{fac.recast}} in package \pkg{dae}, \code{\link{factor}}.}
\examples{
## set up a factor and convert it to a numeric vector
a <- factor(rep(1:3, 4))
x <- as.numfac(a)
}
\keyword{factor}
\keyword{manip}