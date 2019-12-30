\name{projector}
\alias{projector}
\title{Create projectors}
\description{The class "\code{\link[=projector-class]{projector}}" is the 
             subclass of the class "\code{\link[=matrix-class]{matrix}}" 
             in which matrices are square, symmetric and idempotent. 
             
             The function \code{projector} tests whether a \code{\link{matrix}} 
             satisfies these criteria and if it does creates a 
             "\code{\link[=projector-class]{projector}}" object, computing the 
             projector's degrees of freedom and adding them to the object.}
\usage{projector(Q)}
\arguments{
 \item{Q}{The \code{\link{matrix}} to be made into a projector.}
}
\value{An object of Class "\code{\link[=projector-class]{projector}}" that 
       consists of a square, summetric,  idempotent \code{\link{matrix}} and 
       degrees of freedom (rank) of the \code{matrix.}}
\details{In checking that the \code{\link{matrix}} is square, symmetric and 
         idempotent, the equality of the \code{\link{matrix}} with either its 
         transpose or square is tested. In this, a difference in elements is 
         considered to be zero if it is less than \code{daeTolerance}, which is 
         initially set to \code{.Machine$double.eps ^ 0.5} (about 1.5E-08). 
         The function \code{\link{set.daeTolerance}} can 
         be used to change \code{daeTolerance}.}
\author{Chris Brien}
\seealso{\code{\link{degfree}}, \code{\link{correct.degfree}} in package \pkg{dae}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## check that it is a valid projector
is.projector(proj.m)
}
\keyword{array}
\keyword{projector}
