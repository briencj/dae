\name{efficiency.criteria}
\alias{efficiency.criteria}
\title{Computes efficiency criteria from a set of efficiency factors}
\description{Computes efficiency criteria from a set of efficiency factors.}
\usage{efficiency.criteria(efficiencies)}
\arguments{
 \item{efficiencies}{A \code{numeric} containing a set of efficiency factors.}
}
\value{A \code{list} whose components are \code{aefficiency}, 
       \code{mefficiency}, \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
       \code{order} and \code{dforthog}.}
\details{The \code{aefficiency} criterion is the harmonic mean of the nonzero 
         efficiency factors. The \code{mefficiency} criterion is the  mean of 
         the nonzero efficiency factors. The \code{eefficiency} criterion is the 
         minimum of the nonzero efficiency factors. The \code{sefficiency} 
         criterion is the variance of the nonzero efficiency factors. The 
         \code{xefficiency} is the maximum of the efficiency factors. The 
         \code{order} is the order of balance and is the number of unique 
         nonzero efficiency factors. The \code{dforthog} is the number of 
         efficiency factors that are equal to one.}
\author{Chris Brien}
\seealso{\code{\link{proj2.efficiency}}, \code{\link{proj2.eigen}}, \code{\link{proj2.combine}} in package \pkg{dae}, 
\code{\link{eigen}}.

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

## save intrablock efficiencies
eff.inter <- proj2.efficiency(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])

## calculate efficiency criteria
efficiency.criteria(eff.inter)
}
\keyword{array}
\keyword{design}
\keyword{projector}
