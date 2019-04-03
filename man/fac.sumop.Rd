\name{fac.sumop}
\alias{fac.sumop}
\title{computes the summation matrix that produces sums corresponding to a (generalized) factor}
\description{Computes the matrix that produces the sums 
             corresponding to a (generalized) \code{\link{factor}}.}
\usage{fac.sumop(factor)}
\arguments{
 \item{factor}{The (generalized) \code{\link{factor}} whose sums the summation matrix
               computes from an observation-length vector.}
}
\value{A symmetric matrix.}
\details{
The design matrix \bold{X} for a (generalized) \code{\link{factor}} is formed with a 
column for each \code{level} of the (generalized) \code{\link{factor}}, this column 
being its indicator variable. The summation matrix is formed as 
\code{X \%*\% t(X)}.

A generalized \code{\link{factor}} is a \code{\link{factor}} formed from the combinations of 
the \code{levels} of several original \code{\link{factor}s}. Generalized \code{\link{factor}s} 
can be formed using \code{\link{fac.combine}}.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}}, \code{\link{fac.meanop}} in package \pkg{dae}.}

\examples{
## set up a two-level factoir and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a generlaized factor whose levels are the combinations of A and B
AB <- fac.combine(list(A,B))

## obtain the operator that computes the AB means from a vector of length 12
S.AB <- fac.sumop(AB)
}
\keyword{array}
\keyword{projector}
