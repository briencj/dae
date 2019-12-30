\name{rmvnorm}
\alias{rmvnorm}
\title{generates a vector of random values from a multivariate normal distribution}
\description{Generates a vector of random values from an n-dimensional 
             multivariate normal distribution whose mean is given by the 
             n-vector \code{mean} and variance by the 
             n x n symmetric matrix \code{V}. It uses the method described by 
             Ripley (1987, p.98)}
\usage{rmvnorm(mean, V, method = 'eigenanalysis')}
\arguments{
 \item{mean}{The mean vector of the multivariate normal distribution from which 
               the random values are to be generated.}
 \item{V}{The variance matrix of the multivariate normal distribution from which 
               the random values are to be generated.}
 \item{method}{The method used to decompose the variance matrix in producing a 
               a matrix to transform the iid standard normal values. The two 
               methods available are \code{'eigenanalysis'} and \code{'choleski'}, 
               where only the first letter of each option is obligatory. 
               The default \code{method} is eigenanalysis, which is slower but 
               is likely to be more stable than Choleski decomposition.}
}
\value{A \code{\link{vector}} of length n, equal to the length of \code{mean}.}
\details{
The method is:
 a) uses either the eigenvalue or Choleski decomposition of the variance matrix, 
    \code{V}, to form the matrix that transforms an iid vector of values to a 
    vector with variance \code{V};
 b) generate a vector of length equal to \code{mean} of standard normal values;
 c) premultiply the vector of standard normal values by the transpose of the 
    upper triangular factor and, to the result, add \code{mean}.
}
\author{Chris Brien}
\references{Ripley, B. D. (1987) \emph{Stochastic simulation}. Wiley, New York.}
\seealso{\code{\link{fac.ar1mat}},  \code{\link{fac.vcmat}}, 
         in package \pkg{dae}, \code{\link{rnorm}}, and \code{\link{chol}}.}
\examples{
## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## generate random values from a multivariate normal for which 
#the mean is 20 for all variables and 
#the variance matrix has random effects for factor A, ar1 pattern for B and 
#residual random variation
mean <- rep(20, 12)
V <- fac.vcmat(A, 5) + fac.ar1mat(B, 0.6) + 2*mat.I(12)
y <- rmvnorm(mean, V)
}
\keyword{datagen}
