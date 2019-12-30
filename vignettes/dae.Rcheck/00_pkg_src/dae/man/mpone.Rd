\name{mpone}
\alias{mpone}
\title{Converts the first two levels of a factor into the numeric values
       -1 and +1}
\description{Converts the first two \code{levels} of a \code{\link{factor}} into the numeric
             values -1 and +1.}
\usage{mpone(factor)}
\arguments{
 \item{factor}{The \code{\link{factor}} to be converted.}
}
\value{A numeric \code{vector}.}
\section{Warning}{If the \code{\link{factor}} has more than two \code{levels} they will 
      be coerced to numeric values.}
\author{Chris Brien}
\seealso{\code{\link{mpone}} in package \pkg{dae}, \code{\link{factor}}, 
\code{\link{relevel}}.}
\examples{
## generate all combinations of two two-level factors
mp <- c("-", "+")
Frf3.trt <- fac.gen(list(A = mp, B = mp))

## add factor C, whose levels are the products of the levles of A and B
Frf3.trt$C <- factor(mpone(Frf3.trt$A)*mpone(Frf3.trt$B), labels = mp)
}
\keyword{factor}
\keyword{manip}