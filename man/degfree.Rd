\name{degfree}
\alias{degfree}
\alias{degfree<-}
\title{Degrees of freedom extraction and replacement}
\description{Extracts the degrees of freedom from or replaces them in an object 
of class "\code{\link[=projector-class]{projector}}".}
\usage{degfree(object)

degfree(object) <- value}
\arguments{
 \item{object}{An object of class "\code{\link[=projector-class]{projector}}" 
               whose degrees of freedom are to be extracted or replaced.}
 \item{value}{An integer to which the degrees of freedom are to be set or 
              an object of class "\code{\link[=projector-class]{projector}}" or "matrix" from which the degrees 
              of freedom are to be calulated.}
}
\details{There is no checking of the correctness of the degrees of freedom, 
         either already stored or as a supplied integer value. This can be done using 
         \code{\link{correct.degfree}}.

When the degrees of freedom of the projector are to be calculated, they are obtained 
as the number of nonzero eigenvalues. An eigenvalue is regarded as zero if it is 
less than \code{daeTolerance}, which is initially set to \code{.Machine$double.eps ^ 0.5} 
(about 1.5E-08). The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.}
\value{An object of class "\code{\link[=projector-class]{projector}}" that consists 
       of a square, summetric, idempotent matrix and degrees of freedom (rank) of the matrix.}
\author{Chris Brien}
\seealso{\code{\link{correct.degfree}}, \code{\link{projector}} in package \pkg{dae}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## coerce to a projector
proj.m <- projector(m)

## extract its degrees of freedom
degfree(proj.m)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom and print the projector
degfree(proj.m) <- proj.m
print(proj.m)
}
\keyword{array}
\keyword{projector}
