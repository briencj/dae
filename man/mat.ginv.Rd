\name{mat.ginv}
\alias{mat.ginv}
\title{Computes the generalized inverse of a matrix}
\description{Computes the Moore-Penrose generalized inverse of a matrix.}
\usage{mat.ginv(x, tol = .Machine$double.eps ^ 0.5)}
\arguments{
 \item{x}{A \code{\link{matrix}} whose generalized inversed is to be computed.}
 \item{tol}{A \code{\link{numeric}} specifying the relative tolerance to 
            determine whether an eigenvalue of \code{x} is nonzero.}
}
\value{A \code{matrix}. An \code{NA} is returned if \code{\link{svd}} fails 
       during the compution of the generalized inverse.}
\author{Chris Brien}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)
## Compute the projector for a linear trend across Blocks
PBIBD2.lay <- within(PBIBD2.lay,
                     {
                       cBlock <- as.numfac(Block)
                       cBlock <- cBlock - mean(unique(cBlock))
                     })
X <- model.matrix(~ cBlock, data = PBIBD2.lay)
Q.cB <- projector((X \%*\% mat.ginv(t(X) \%*\% X) \%*\% t(X)))
}
\keyword{matrix}