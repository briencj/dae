\name{projs.combine.p2canon}
\alias{projs.combine.p2canon}
\title{Extract, from a p2canon object, the projectors that give the combined 
       canonical decomposition}
\description{Extracts, from a p2canon object obtained using 
             \code{\link{projs.2canon}}, the projectors that give the combined 
             canonical decomposition of two sets of projectors 
             (Brien and Bailey, 2009).}
\usage{projs.combine.p2canon(object)}
\arguments{
 \item{object}{A \code{list} of class \code{p2canon} produced by \code{projs.2canon}.}
}
\value{A \code{list}, each of whose components is a projector in the decomposition.}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{projs.2canon}}, \code{\link{proj2.eigen}}, \code{\link{proj2.combine}} in package \pkg{dae}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
UcombineT <- projs.combine.p2canon(unit.trt.p2canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
