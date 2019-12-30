\name{decomp.relate}
\alias{decomp.relate}
\title{Examines the relationship between the eigenvectors for two decompositions}
\description{Two decompositions produced by \code{\link{proj2.eigen}} are compared 
             by computing all pairs of crossproduct sums of eigenvectors from the 
             two decompositions. It is most useful when the calls to 
             \code{\link{proj2.eigen}} have the same Q1.}
\usage{decomp.relate(decomp1, decomp2)}
\arguments{
 \item{decomp1}{A \code{\link{list}} containing components efficiencies and eigenvectors such as 
                is produced by \code{\link{proj2.eigen}}.}
 \item{decomp2}{Another \code{\link{list}} containing components efficiencies and eigenvectors 
                such as is produced by \code{\link{proj2.eigen}}.}
}
\value{A \code{\link{matrix}} that is r1 x r2 where r1 and r2 are the numbers of efficiencies 
       of \code{decomp1} and \code{decomp2}, respectively. The \code{rownames} 
       and \code{columnnames} of the \code{\link{matrix}} are the values of the 
       efficiency factors from \code{decomp1} and \code{decomp2}, respectively.}

\details{
Each element of the r1 x r2 \code{\link{matrix}} is the sum of crossproducts of a pair of 
eigenvectors, one from each of the two decompositions. A sum is regarded as zero 
if it is less than \code{daeTolerance}, which is initially set to 
.Machine$double.eps ^ 0.5 (about 1.5E-08). The function \code{\link{set.daeTolerance}} 
can be used to change \code{daeTolerance}.}
\author{Chris Brien}
\seealso{\code{\link{proj2.eigen}}, \code{\link{proj2.combine}} in package \pkg{dae}, 
\code{\link{eigen}}.}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## obtain intra- and inter-block decompositions
decomp.inter <- proj2.eigen(unit.struct$Q[["Block"]], trt.struct$Q[["trt"]])
decomp.intra <- proj2.eigen(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])

## check that intra- and inter-block decompositions are orthogonal
decomp.relate(decomp.intra, decomp.inter) 
}
\keyword{array}
\keyword{design}
\keyword{projector}
