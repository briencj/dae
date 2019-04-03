\name{print.projector}
\alias{print.projector}
\alias{print,projector-method}
\title{Print projectors}
\description{Print an object of class "\code{\link[=projector-class]{projector}}", 
             displaying the matrix and its degrees of freedom (rank).}
\usage{\method{print}{projector}(x, ...)}
\arguments{
 \item{x}{The object of class "\code{\link[=projector-class]{projector}}" to be printed.}
 \item{...}{Further arguments passed to or from other methods.}
}
\author{Chris Brien}
\seealso{\code{\link{print}}, \code{\link{print.default}}, \code{\link{show}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## print the object either using the Method function, the generic function or show
print.projector(proj.m)
print(proj.m)
proj.m
}
\keyword{array}
\keyword{projector}
