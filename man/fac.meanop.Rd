\name{fac.meanop}
\alias{fac.meanop}
\title{computes the projection matrix that produces means}
\description{Computes the symmetric projection matrix that produces the means 
             corresponding to a (generalized) \code{\link{factor}}.}
\usage{fac.meanop(factor)}
\arguments{
 \item{factor}{The (generalized) \code{\link{factor}} whose means the projection matrix
               computes from an observation-length vector.}
}
\value{A \code{\link{projector}} containing the symmetric, projection matrix 
and its degrees of freedom.}
\details{
The design matrix \bold{X} for a (generalized) \code{\link{factor}} is formed with a 
column for each \code{level} of the (generalized) \code{\link{factor}}, this column 
being its indicator variable. The projection matrix is formed as 
\code{X \%*\% (1/diag(r) \%*\% t(X)}, where \code{r} is the \code{vector} of 
\code{levels} replications.

A generalized \code{\link{factor}} is a \code{\link{factor}} formed from the 
combinations of the \code{levels} of several original \code{\link{factor}s}. 
Generalized \code{\link{factor}s} can be formed using \code{\link{fac.combine}}.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}}, \code{\link{projector}}, \code{\link{degfree}}, 
         \code{\link{correct.degfree}}, \code{\link{fac.sumop}} in package \pkg{dae}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a generalized factor whose levels are the combinations of A and B
AB <- fac.combine(list(A,B))

## obtain the operator that computes the AB means from a vector of length 12
M.AB <- fac.meanop(AB)
}
\keyword{array}
\keyword{projector}
