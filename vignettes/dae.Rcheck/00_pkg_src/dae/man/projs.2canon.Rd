\name{projs.2canon}
\alias{projs.2canon}
\title{A canonical analysis of the relationships between two sets of projectors}
\description{Computes the canonical efficiency factors for the joint 
       decomposition of two structures or sets of mutually orthogonally 
       projectors (Brien and Bailey, 2009), orthogonalizing projectors in the 
       Q2 \code{list} to those earlier in the \code{list} of projectors with 
       which they are partially aliased. The results can be summarized in the 
       form of a skeleton ANOVA table.}
\usage{projs.2canon(Q1, Q2)}
\arguments{
 \item{Q1}{A \code{list} whose components are objects of class "\code{\link[=projector-class]{projector}}".}
 \item{Q2}{A \code{list} whose components are objectsof class "\code{\link[=projector-class]{projector}}".}
}
\value{A \code{\link{p2canon.object}}.}
\details{Two loops, one nested within the other, are performed. The first cycles 
         over the components of \code{Q1} and the nested loop cycles over the 
         components of \code{Q2}. The joint decomposition of the two projectors 
         in each cycle, one from \code{Q1} (say \code{Q1[[i]]}) and the other 
         from \code{Q2} (say \code{Q2[[j]]}) is obtained using 
         \code{\link{proj2.combine}}. 
         In particular, the nonzero canonical efficiency factors for the joint 
         decomposition of the two projectors is obtained. The nonzero canonical 
         efficiency factors are the nonzero eigenvalues of 
         \code{Q1[[i]] \%*\% Q2[[j]] \%*\% Q1[[i]]} (James and Wilkinson, 1971). 
         An eigenvalue is regarded as zero if it is less than 
         \code{daeTolerance}, which is initially set to 
         \code{.Machine$double.eps ^ 0.5} (about 1.5E-08). The function 
         \code{\link{set.daeTolerance}} can be used to change 
         \code{daeTolerance}.
         
         However, a warning occurs if any pair of Q2 projectors (say 
         \code{Q2[[j]]} and \code{Q2[[k]]}) do not have adjusted orthgonality 
         with respect to any Q1 projector (say \code{Q1[[i]]}), because they are 
         partially aliased. That is, if \code{Q2[[j]] \%*\% Q1[[i]] \%*\% Q2[[k]]} 
         is nonzero for any pair of different Q2 projectors and any 
         Q1 projector. When it is nonzero, the projector for the later term in 
         the list of projectors is orthogonalized to the projector that is 
         earlier in the list. A list o such projectors is returned in the 
         \code{aliasing} component of the \code{\link{p2canon.object}}.}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.

           James, A. T. and Wilkinson, G. N. (1971) Factorization of the residual 
           operator and canonical decomposition of nonorthogonal factors in the 
           analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{summary.p2canon}}, \code{\link{efficiencies.p2canon}}, 
         \code{\link{projs.combine.p2canon}}, \code{\link{pstructure}} , \cr
         \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, 
         \code{\link{proj2.eigen}}, \code{\link{efficiency.criteria}} 
         in package \pkg{dae}, \code{\link{eigen}}.

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
