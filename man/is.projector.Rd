\name{is.projector}
\alias{is.projector}
\title{Tests whether an object is a valid object of class projector}
\description{Tests whether an  \code{object} is a valid object of class 
"\code{\link[=projector-class]{projector}}".}
\usage{is.projector(object)}
\arguments{
 \item{object}{The \code{\link{matrix}} to be made into a projector.}
}
\details{The function \code{is.projector} tests whether the object consists of a 
         \code{\link{matrix}} that is square, symmetric and idempotent. In checking 
         symmetry and idempotency, the equality of the matrix with either its 
         transpose or square is tested. In this, a difference in elements is 
         considered to be zero if it is less than \code{daeTolerance}, which is 
         initially set to \code{.Machine$double.eps ^ 0.5} (about 1.5E-08). 
         The function \code{\link{set.daeTolerance}} can 
         be used to change \code{daeTolerance}.}
\value{\code{TRUE} or \code{FALSE} depending on whether the object is a valid object of class 
       "\code{\link[=projector-class]{projector}}".}
\section{Warning}{The degrees of freedom are not checked. \code{\link{correct.degfree}} 
can be used to check them.}
\author{Chris Brien}
\seealso{\code{\link{projector}}, \code{\link{correct.degfree}} in package \pkg{dae}.

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
