\name{set.daeTolerance}
\alias{set.daeTolerance}
\title{Sets the values of daeTolerance for the package dae}
\description{A function that sets the values such that, in \pkg{dae} functions, 
             values less than it are considered to be zero. The values are stored 
             in a \code{vector} named \code{daeTolerance} in the \code{daeEnv} 
             environment. The \code{vector} is of length two and, initially, both 
             values are set to \code{.Machine$double.eps ^ 0.5} (about 1.5E-08). 
             One value is named \code{element.tol} and is used for elements of 
             matrices; the second is named \code{element.eigen} and is used for 
             eigenvalues and quantities based on them, such as efficiency factors.}
\usage{set.daeTolerance(element.tol=NULL, eigen.tol=NULL)}
\arguments{
 \item{element.tol}{The value to to which the first element of the \code{daeTolerance} 
                    \code{vector} is to be set. If more than one value is 
                    supplied, only the first value is used. }
 \item{eigen.tol}{The value to to which the second element of the \code{daeTolerance} 
                    \code{vector} is to be set. If more than one value is 
                    supplied, only the first value is used. }
}
\value{The \code{vector} \code{daeTolerance} is returned invisibly.}
\author{Chris Brien}
\seealso{\code{\link{get.daeTolerance}}.}
\examples{
## set daeTolerance.
set.daeTolerance(1E-04, 1E-08)
}
\keyword{manip}
\keyword{projector}