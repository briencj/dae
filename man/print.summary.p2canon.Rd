\name{print.summary.p2canon}
\alias{print.summary.p2canon}
\title{Prints the values in an \code{\link{summary.p2canon}} object}
\description{Prints a \code{summary.p2canon} object, which is also a 
             \code{data.frame}, in a pretty format.}
\usage{\method{print}{summary.p2canon}(x, ...)}
\arguments{
 \item{x}{A \code{summary.p2canon} object.}
 \item{...}{further arguments passed to \code{print}.}
}
\value{No value is returned.}


\author{Chris Brien}
\seealso{\code{\link{summary.p2canon}}}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain projectors using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and print summary
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summ <- summary(unit.trt.p2canon)
print(summ)
}
\keyword{design}
\keyword{projector}