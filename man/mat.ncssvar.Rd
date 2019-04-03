\name{mat.ncssvar}
\alias{mat.ncssvar}
\title{Calculates the variance matrix of the random effects for a natural 
       cubic smoothing spline}
\description{Calculates the variance matrix of the random effects for a 
             natural cubic smoothing spline. It is the tri-diagonal matrix 
             \eqn{\bold{G}_s}{G_s} given by Verbyla et al., (1999) multiplied by 
             the variance component for the random spline effects.}
\usage{
mat.ncssvar(sigma2s = 1, knot.points, print = FALSE)}
\arguments{
 \item{sigma2s}{A \code{\link{numeric}} giving the value of the variance component
           for the random spline effects. The smoothing parameter is then the inverse 
           of the ratio of this component to the residual variance.}
 \item{knot.points}{A \code{\link{numeric}} giving the values of the knots point 
           used in fitting the spline. These must be orderd in increasing order.}
 \item{print}{A \code{\link{logical}} indicating whether to print the matrix.}
}
\value{A \code{\link{matrix}} containing the variances and covariances of the 
       random spline effects.}
\author{Chris Brien}
\seealso{\code{\link{Zncsspline}}.}
\references{Verbyla, A. P., Cullis, B. R., Kenward, M. G., and Welham, S. J. (1999). 
            The analysis of designed experiments and longitudinal data by using 
            smoothing splines (with discussion). \emph{Journal of the Royal 
            Statistical Society, Series C (Applied Statistics)}, \bold{48}, 269-311.}
\examples{
Gs <- mat.ncssvar(knot.points = 1:10)
}
\keyword{array}
\keyword{design}
