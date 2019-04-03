\name{efficiencies}
\alias{efficiencies}
\alias{efficiencies.pcanon}
\alias{efficiencies.p2canon}
\title{Extracts the canonical efficiency factors from a \code{\link{pcanon.object}} or a \code{\link{p2canon.object}}.}
\description{Produces a \code{list} containing the canonical efficiency factors 
             for the joint decomposition of two or more sets of projectors 
             (Brien and Bailey, 2009) obtained using \code{\link{designAnatomy}} or 
             \code{\link{projs.2canon}}.}
\usage{
\method{efficiencies}{pcanon}(object, which = "adjusted", ...)
\method{efficiencies}{p2canon}(object, which = "adjusted", ...)}
\arguments{
 \item{object}{A \code{\link{pcanon.object}} or an object of class \code{p2canon} produced by \code{\link{projs.2canon}}.}
 \item{which}{A character string, either \code{adjusted} or \code{pairwise}. For 
              \code{adjusted}, the canonical efficiency factor are adjusted for other 
              projectors from from the same set. For \code{pairwise}, they are the unadjusted 
              canonical efficiency factors between pairs of projectors consisting of one 
              projector from each of two sets.}
 \item{...}{Further arguments passed to or from other methods. Unused at present.}
}
\value{For a \code{pcanon.object}, a \code{list} with a component for each component of 
       \code{object}, except for the last component -- for more information about the components 
       see \code{\link{pcanon.object}} .

       For a \code{p2canon} object, a \code{list} with a component for each element of the \code{Q1} 
       argument from \code{\link{projs.2canon}}. Each component is \code{list}, each its components 
       corresponding to an element of the \code{Q2} argument from \code{\link{projs.2canon}}}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{designAnatomy}}, \code{\link{summary.pcanon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, \code{\link{proj2.eigen}},  \cr 
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

##obtain combined decomposition using designAnatomy and get the efficiencies
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), data = PBIBD2.lay)
efficiencies.pcanon(unit.trt.canon)

##obtain the projectors for each formula using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition projs.2canon and get the efficiencies
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
efficiencies.p2canon(unit.trt.p2canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
