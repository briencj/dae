\name{proj2.efficiency}
\alias{proj2.efficiency}
\title{Computes the canonical efficiency factors for the joint decomposition of 
       two projectors}
\description{Computes the canonical efficiency factors for the joint    
             decomposition of two projectors (James and Wilkinson, 1971).}
\usage{proj2.efficiency(Q1, Q2)}
\arguments{
 \item{Q1}{An object of class "\code{\link[=projector-class]{projector}}".}
 \item{Q2}{An object of class "\code{\link[=projector-class]{projector}}".}
}
\value{A \code{vector} containing the nonzero canonical efficiency factors.}
\details{
The nonzero canonical efficiency factors are the nonzero eigenvalues of 
Q1 \%*\% Q2 \%*\% Q1 (James and Wilkinson, 1971). An eigenvalue is regarded as 
zero if it is less than \code{daeTolerance}, which is initially set to 
\code{.Machine$double.eps ^ 0.5} (about 1.5E-08). 
The function \code{\link{set.daeTolerance}} can be used to change 
\code{daeTolerance}.}
\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{efficiency.criteria}}, \code{\link{proj2.eigen}}, \code{\link{proj2.combine}} in package \pkg{dae}, 
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
eff.intra <- proj2.efficiency(unit.struct$Q[["Block"]], trt.struct$Q[["trt"]])
}
\keyword{array}
\keyword{design}
\keyword{projector}
