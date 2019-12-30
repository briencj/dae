\name{fac.ar1mat}
\alias{fac.ar1mat}
\title{forms the ar1 correlation matrix for a (generalized) factor}
\description{Form the correlation matrix for a (generalized) factor where the 
  correlation between the levels follows an autocorrelation of order 1 (ar1) 
  pattern.}
\usage{fac.ar1mat(factor, rho)}
\arguments{
 \item{factor}{The (generalized) \code{\link{factor}} for which the correlation between 
               its levels displays an ar1 pattern.}
 \item{rho}{The correlation parameter for the ar1 process.}
}
\value{An n x n \code{\link{matrix}}, where n is the length of the 
       \code{\link{factor}}.}
\details{
The method is:
  a) form an \code{n x n} matrix of all pairwise differences in the numeric values 
     corresponding to the observed levels of the factor by taking the 
     difference between the following two n x n matrices are equal: 1) each row 
     contains the numeric values corresponding to the observed levels of the 
     factor, and 2) each column contains the numeric values corresponding to 
     the observed levels of the factor, 
  b) replace each element of the pairwise difference matrix with rho raised to 
     the absolute value of the difference.
}
\author{Chris Brien}
\seealso{\code{\link{fac.vcmat}}, \code{\link{fac.meanop}}, 
         \code{\link{fac.sumop}} in package \pkg{dae}.}
\examples{
## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a 12 x 12 ar1 matrix corrresponding to B
ar1.B <- fac.ar1mat(B, 0.6)
}
\keyword{array}
