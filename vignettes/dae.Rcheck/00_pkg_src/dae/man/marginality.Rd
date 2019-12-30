\name{marginality}
\alias{marginality.pcanon}
\alias{marginality.pstructure}
\alias{marginality}
\title{Extracts the marginality matrix (matrices) from a \code{\link{pstructure.object}} or a \code{\link{pcanon.object}}.}
\description{Produces (i) a marginality \code{\link{matrix}} for the \code{formula} in a call to 
             \code{\link{pstructure.formula}} or (ii) a \code{list} containing the marginlity matrices, one for each 
             \code{formula} in the \code{formulae} argument of a call to 
             \code{\link{designAnatomy}}.

             A marginality matrix for a set of terms is a square \code{\link{matrix}} with 
             a row and a column for each ternon-aliased term. Its elements are zeroes and ones, 
             the entry in the ith row and jth column indicates whether or not the ith term is 
             marginal to the jth term i.e. the column space of the ith term is a subspace of 
             that for the jth term and so the source for the jth term will be orthogonal to 
             that for the ith term.
}
\usage{
\method{marginality}{pstructure}(object, ...)
\method{marginality}{pcanon}(object, ...)}
\arguments{
 \item{object}{A \code{\link{pstructure.object}} produced by \code{\link{pstructure.formula}} or  
               \code{\link{pcanon.object}} produced by \code{\link{designAnatomy}}.}
 \item{...}{Further arguments passed to or from other methods. Unused at present.}
}
\value{If \code{object} is a \code{\link{pstructure.object}} then a \code{matrix} containing 
       the marginality matrix for the terms obtained from the \code{formuula} in the call to 
       \code{\link{pstructure.formula}}. 

       If \code{object} is a \code{\link{pcanon.object}} then a \code{list} with a 
        component for each \code{formula}, each component having a marginality matrix that 
       corresponds to one of the formulae in the call to \code{\link{designAnatomy}}. The
       components of the \code{list} will have the same names as the componeents of the 
       \code{formulae} \code{list} and so will be unnamed if the components of the latter 
       \code{list} are unnamed.}
\author{Chris Brien}
\seealso{\code{\link{pstructure.formula}}, \code{\link{designAnatomy}}, \code{\link{summary.pcanon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, \code{\link{proj2.eigen}},  \cr 
         \code{\link{pstructure}} in package \pkg{dae}, 
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

##obtain pstructure.object and extract marginality matrix
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
unit.marg <- marginality(unit.struct)

##obtain combined decomposition and extract marginality matrices
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), data = PBIBD2.lay)
marg <- marginality(unit.trt.canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
