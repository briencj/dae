\name{fac.vcmat}
\alias{fac.vcmat}
\title{forms the variance matrix for the variance component of a (generalized) factor}
\description{Form the variance matrix for a (generalized) factor whose effects 
  for its different levels are independently and identically distributed, with 
  their variance given by the variance component; elements of the matrix will 
  equal either zero or sigma2 and displays compound symmetry.}
\usage{fac.vcmat(factor, sigma2)}
\arguments{
 \item{factor}{The (generalized) \code{\link{factor}} for which the variance 
               matrix is required.}
 \item{sigma2}{The variance component, being the of the random effects for the 
               factor.}
}
\value{An n x n \code{\link{matrix}}, where n is the length of the 
       \code{\link{factor}}.}
\details{
The method is:
  a) form the n x n summation or relationship matrix whose elements are equal 
     to zero except for those elements whose corresponding elements in the 
     following two n x n matrices are equal: 1) each row contains the numeric 
     values corresponding to the observed levels of the factor, and 2) each 
     column contains the numeric values corresponding to the observed levels 
     of the factor, 
  b) multiply the summation matrix by sigma2.
}
\author{Chris Brien}
\seealso{\code{\link{fac.ar1mat}},  \code{\link{fac.meanop}}, 
         \code{\link{fac.sumop}} in package \pkg{dae}.}
\examples{
## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a 12 x 12 ar1 matrix corrresponding to B
vc.B <- fac.vcmat(B, 2)
}
\keyword{array}
