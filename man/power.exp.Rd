\name{power.exp}
\alias{power.exp}
\title{Computes the power for an experiment}
\description{Computes the power for an experiment.}
\usage{
power.exp(rm=5., df.num=1., df.denom=10., delta=1., sigma=1.,
          alpha=0.05, print=FALSE)}
\arguments{
 \item{rm}{The number of observations used in computing a mean.}
 \item{df.num}{The degrees of freedom of the numerator of the F for testing
               the term involving the means.}
 \item{df.denom}{The degrees of freedom of the denominator of the F for
                 testing the term involving the means.}
 \item{delta}{The true difference between a pair of means.}
 \item{sigma}{The population standard deviation.}
 \item{alpha}{The significance level to be used.}
 \item{print}{\code{TRUE} or \code{FALSE} to have or not have a table of power 
              calculation details printed out.}
}
\value{A single \code{numeric} value containing the computed power.}
\author{Chris Brien}
\seealso{\code{\link{no.reps}}, \code{\link{detect.diff}} in package \pkg{dae}.}
\examples{
## Compute power for a randomized complete block design with four treatments 
## and five blocks. 
rm <- 5
power.exp(rm = rm, df.num = 3, df.denom = 3 * (rm - 1), delta = 5,
          sigma = sqrt(20),print = TRUE)
}
\keyword{design}
