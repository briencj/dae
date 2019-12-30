\name{correct.degfree}
\alias{correct.degfree}
\title{Check the degrees of freedom in an object of class projector}
\description{Check the degrees of freedom in an object of class "\code{\link[=projector-class]{projector}}".}
\usage{correct.degfree(object)}
\arguments{
 \item{object}{An object of class "\code{\link[=projector-class]{projector}}" whose degrees of freedom are to be checked.}
}
\value{TRUE or FALSE depending on whether the correct degrees of freedom have been
       stored in the object of class "\code{\link[=projector-class]{projector}}".}
\details{
The degrees of freedom of the projector are obtained as its number of nonzero 
eigenvalues. An eigenvalue is regarded as zero if it is less than 
\code{daeTolerance}, which is initially set to.Machine$double.eps ^ 0.5 (about 1.5E-08). 
The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.}
\author{Chris Brien}
\seealso{\code{\link{degfree}}, \code{\link{projector}} in package \pkg{dae}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom
degfree(proj.m) <- 1
    
## check degrees of freedom are correct
correct.degfree(proj.m)
}
\keyword{array}
\keyword{projector}
