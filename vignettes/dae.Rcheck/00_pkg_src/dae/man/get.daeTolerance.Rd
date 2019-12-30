\name{get.daeTolerance}
\alias{get.daeTolerance}
\title{Gets the value of {daeTolerance} for the package dae}
\description{A function that gets the \code{vector} of values such that, in \pkg{dae} 
             functions, values less than it are considered to be zero.}
\usage{get.daeTolerance()}
\value{The \code{vector} of two values for {daeTolerance}, one named \code{element.tol} 
       that is used for elements of matrices and a second named \code{element.eigen} 
       that is used for eigenvalues and quantities based on them, such as efficiency 
       factors.}
\author{Chris Brien}
\seealso{\code{\link{set.daeTolerance}}.}
\examples{
## get daeTolerance.
get.daeTolerance()
}
\keyword{manip}
\keyword{projector}