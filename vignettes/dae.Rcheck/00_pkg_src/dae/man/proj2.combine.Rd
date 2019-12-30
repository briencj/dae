\name{proj2.combine}
\alias{proj2.combine}
\title{Compute the projection and Residual operators for two, possibly 
       nonorthogonal, projectors}
\description{The canonical relationship between a pair of projectors 
             is established by decomposing the range of Q1 into a part that 
             pertains to Q2 and a part that is orthogonal to Q2. It also 
             produces the nonzero canonical efficiency factors for the joint 
             decomposition of Q1 and Q and the corresponding eigenvectors of 
             Q1 (James and Wilkinson, 1971). Q1 and Q2 may be nonorthogonal.}
\usage{proj2.combine(Q1, Q2)}
\arguments{
 \item{Q1}{A symmetric \code{projector} whose range is to be decomposed.}
 \item{Q2}{A symmetric \code{projector} whose range in Q1 is required.}
}
\value{A \code{list} with the following components:
\enumerate{
\item{\bold{efficiencies:} a \code{vector} containing the nonzero canonical efficiency factors;}
\item{\bold{eigenvectors:} an n x r \code{\link{matrix}}, where n is the order of the projectors and 
            r is the number of nonzero canonical efficiency factors; it contains 
            the eigenvectors of Q1 corresponding to the nonzero canonical 
            efficiency factors.} 
\item{\bold{Qconf:} a \code{projector} onto the part of the range of Q1 with 
            which Q2 is confounded;}
\item{\bold{Qres:} a \code{projector} onto the part of the range of Q1 that is 
            orthogonal to the range of Q2.}}}
\details{
The nonzero canonical efficiency factors are the nonzero eigenvalues of 
\code{Q1 \%*\% Q2 \%*\% Q1} (James and Wilkinson, 1971). An eigenvalue is regarded 
as zero if it is less than \code{daeTolerance}, which is initially set to 
\code{.Machine$double.eps ^ 0.5} (about 1.5E-08). 
The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.

The eigenvectors are the eigenvectors of Q1 corresponding to the nonzero canonical 
efficiency factors. The eigenvectors for Q2 can be obtained by premultiplying 
those for Q1 by Q2.

Qres is computed using equation 4.10 from James and Wilkinson (1971), if the number of distinct
canonical efficiency factors is less than 10. If this fails to produce a projector or the number of distinct canonical efficiency factors is 10 or more, equation 5.3 of Payne and Tobias (1992) is used to obtain Qres. In this latter case, \code{Qres = Q1 - Q1 \%*\% ginv(Q2 \%*\% Q1 \%*\% Q2) \%*\% Q1}. Qconf is obtained by subtracting Qres from Q1.}

\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279--294.
            
            
            Payne, R. W. and R. D. Tobias (1992). General balance, combination of 
            information and the analysis of covariance. 
            \emph{Scandinavian Journal of Statistics}, \bold{19}, 3--23.}
\author{Chris Brien}
\seealso{\code{\link{proj2.eigen}}, \code{\link{proj2.efficiency}}, \code{\link{decomp.relate}} 
in package \pkg{dae}.

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

## obtain the projection operators for the interblock analysis
PBIBD2.Bops <- proj2.combine(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])
Q.B.T <- PBIBD2.Bops$Qconf
Q.B.res <- PBIBD2.Bops$Qres

## demonstrate their orthogonality
is.allzero(Q.B.T \%*\% Q.B.res)
}
\keyword{array}
\keyword{design}
\keyword{projector}
