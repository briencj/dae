\name{fac.divide}
\alias{fac.divide}
\title{Divides a factor into several individual factors}
\description{Takes a \code{\link{factor}} and divides it into several individual 
 \code{\link{factor}s} as if the \code{levels} in the original \code{combined.factor} are numbered from one to its number of levels and correspond 
 to the numbering of the \code{levels} combinations of the individual 
 \code{\link{factor}s} when these are arranged in standard or Yates order.}
\usage{fac.divide(combined.factor, factor.names, order="standard")}
\arguments{
 \item{combined.factor}{A \code{\link{factor}} that is to be divided into the individual 
       \code{\link{factor}s} listed in \code{factor.names}.}
 \item{factor.names}{A \code{\link{list}} of \code{\link{factor}s} to be formed. The names in 
       the \code{\link{list}} are the names of the \code{\link{factor}s} and the component of 
       a name is either a) a single \code{numeric} value that is the number of 
       \code{levels}, b) a numeric \code{vector} that contains the \code{levels} 
       of the \code{\link{factor}}, or c) a character \code{vector} that contains the 
       \code{labels} of the \code{levels} of the \code{\link{factor}}.}
  \item{order}{Either \code{standard} or \code{yates}. The order in which the 
        \code{levels} combinations of the 
        \code{\link{factor}s} in \code{factor.names} are to be considered as numbered; 
        \code{standard} numbers them as if they are arranged in standard order, 
        that is with the first factor moving slowest and the last factor 
        moving fastest; \code{yates} numbers them 
        as if they are arranged in Yates order, that is with the first factor 
        moving fastest and last factor moving slowest.}
}
\value{A \code{\link{data.frame}} whose columns consist of the \code{\link{factor}s} listed in 
  \code{factor.names} and whose values have been computed from the combined 
  \code{\link{factor}}. All the \code{\link{factor}s} will be of the same length.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}} in package \pkg{dae}.}
\note{A single \code{\link{factor}} name may be supplied in the \code{\link{list}} in which case 
      a \code{\link{data.frame}} is produced that contains the single \code{\link{factor}} 
      computed from the numeric \code{vector}. This may be useful when calling 
      this \code{function} from others.}
\examples{
## generate a small completely randomized design for 6 treatments 
n <- 12
CRD.unit <- list(Unit = n)
treat <- factor(rep(1:4, each = 3))
CRD.lay <- designRandomize(allocated = treat, recipient = CRD.unit, seed=956)

## divide the treatments into two two-level factor A nd B
CRD.facs <- fac.divide(CRD.lay$treat, factor.names = list(A = 2, B = 2))
}
\keyword{factor}
\keyword{manip}