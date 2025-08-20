\name{Zncsspline}
\alias{Zncsspline}
\title{Calculates the design matrix for fitting the random component of a natural 
       cubic smoothing spline}
\description{Calculates the design matrix, \eqn{\bold{Z}_s}{Z_S}, of the random effects for a 
             natural cubic smoothing spline as described by Verbyla et al., (1999). 
             An initial design matrix, 
             \eqn{\bold{\Delta} \bold{\Delta}^{-1} \bold{\Delta}}{\bold{Delta} \bold{Delta}^{-1} \bold{Delta}}, 
              based on the knot points is computed. It can 
             then be post multiplied by a power of the tri-diagonal matrix 
             \eqn{\bold{G}_s}{G_S}, \eqn{\bold{G}_s}{G_S} being proportional to the 
             assumed variance matrix of the random spline effects. If the power is 
             set to 0.5, then the random spline effects based on the resulting design 
             matrix \eqn{\bold{Z}_s}{Z_S} are now independent with variance 
             \eqn{\sigma_s^2}{sigma_s^2}. The variance component that estimates 
             \eqn{\sigma_s^2}{sigma_s^2} will then be a variance ratio and the 
             smoothing parameter is the inverse of the ratio of this variance 
             component to the residual variance.}
\usage{
Zncsspline(knot.points, Gpower = 0, print = FALSE)}
\arguments{
 \item{knot.points}{A \code{\link{numeric}} giving the values of the knot points 
           to be used in fitting the spline. These must be orderd in increasing order.}
 \item{Gpower}{A \code{\link{numeric}} giving the power to take of the tri-diagonal matrix 
               \eqn{\bold{G}_s}{G_s}; the initial design matrix is postmultiplied  
               by the \eqn{\bold{G}_s}{G_s} raised to \code{Gpower} to produce 
               the design matrix \eqn{\bold{Z}_s}{Z_S} that is returned.}
 \item{print}{A \code{\link{logical}} indicating whether to print the 
              \eqn{\bold{\Delta}}{Delta} and \eqn{\bold{G}_s}{G_s} matrices.}
}
\value{A \code{\link{matrix}} that is the design matrix \eqn{\bold{Z}_s}{Z_S}.}
\author{Chris Brien}
\seealso{\code{\link{mat.ncssvar}}.}
\references{Verbyla, A. P., Cullis, B. R., Kenward, M. G., and Welham, S. J. (1999). 
            The analysis of designed experiments and longitudinal data by using 
            smoothing splines (with discussion). \emph{Journal of the Royal 
            Statistical Society, Series C (Applied Statistics)}, \bold{48}, 269-311.}
\examples{
Z <- Zncsspline(knot.points = 1:10, Gpower = 0.5)
}
\keyword{array}
\keyword{design}
