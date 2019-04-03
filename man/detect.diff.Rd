\name{detect.diff}
\alias{detect.diff}
\title{Computes the detectable difference for an experiment}
\description{Computes the delta that is detectable for specified replication, power, alpha.}
\usage{
detect.diff(rm=5, df.num=1, df.denom=10, sigma=1, alpha=0.05, power=0.8, 
            tol = 0.001, print=FALSE)
}
\arguments{
 \item{rm}{The number of observations used in computing a mean.}
 \item{df.num}{The degrees of freedom of the numerator of the F for testing
               the term involving the means.}
 \item{df.denom}{The degrees of freedom of the denominator of the F for
                 testing the term involving the means.}
 \item{sigma}{The population standard deviation.}
 \item{alpha}{The significance level to be used.}
 \item{power}{The minimum power to be achieved.}
 \item{tol}{The maximum difference tolerated between the power required and
            the power computed in determining the detectable difference.}
 \item{print}{\code{TRUE} or \code{FALSE} to have or not have a table of power 
              calculation details printed out.}
}
\value{A single \code{numeric} value containing the computed detectable difference.}
\author{Chris Brien}
\seealso{\code{\link{power.exp}}, \code{\link{no.reps}} in package \pkg{dae}.}
\examples{
## Compute the detectable difference for a randomized complete block design 
## with four treatments given power is 0.8 and alpha is 0.05. 
rm <- 5
detect.diff(rm = rm, df.num = 3, df.denom = 3 * (rm - 1),sigma = sqrt(20))
}
\keyword{design}