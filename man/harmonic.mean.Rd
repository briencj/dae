\name{harmonic.mean}
\alias{harmonic.mean}
\title{Calcuates the harmonic mean.}
\description{A function to calcuate the harmonic mean of a set of nonzero numbers.}
\usage{harmonic.mean(x)}
\arguments{
 \item{x}{An object from whose elements the harmonic mean is to be computed.}
}
\details{
All the elements of \code{x} are tested as being less than \code{daeTolerance}, 
which is initially set to \code{.Machine$double.eps ^ 0.5} (about 1.5E-08). The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.}
\value{A numeric. Returns \code{Inf} if \code{x} contains a value close to zero}
\examples{
y <- c(seq(0.1,1,0.2))
harmonic.mean(y)
}
\keyword{manip}