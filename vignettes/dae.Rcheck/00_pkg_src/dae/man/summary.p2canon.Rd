\name{summary.p2canon}
\alias{summary.p2canon}
\alias{summary,p2canon-method}
\title{Summarize a canonical analysis of the relationships between two sets of projectors}
\description{Produces a summary of the efficiency criteria computed from the 
             canonical efficiency factors for the joint decomposition of two 
             sets of projectors (Brien and Bailey, 2009) obtained using 
             \code{\link{projs.2canon}}. It takes the form of a decomposition or skeleton 
             ANOVA table.}
\usage{
\method{summary}{p2canon}(object, which.criteria = c("aefficiency", "eefficiency", "order"), ...)}
\arguments{
 \item{object}{A \code{list} of class \code{p2canon} produced by 
               \code{projs.2canon}.}
 \item{which.criteria}{A character \code{vector} nominating the efficiency 
                       criteria to be included in the summary. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency},
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}.}
 \item{...}{further arguments affecting the summary produced.}
}
\value{An object of classes \code{summary.p2canon} and \code{data.frame}, whose 
       rows correspond to the pairs of projectors, one from the 
       \code{Q1} argument and the other from the \code{Q2} argument from 
       \code{\link{projs.2canon}}; only pairs with non-zero efficiency factors  
       are included. In addition, a line is included for each nonzero Residual 
       \code{Q1} projector.}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{projs.2canon}}, \code{\link{proj2.efficiency}}, 
         \code{\link{efficiency.criteria}}, \code{\link{proj2.combine}},  
         \code{\link{proj2.eigen}}, \code{\link{pstructure}}, \cr 
         \code{\link{print.summary.p2canon}}in package \pkg{dae}, 
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

##obtain projectors using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summary(unit.trt.p2canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
