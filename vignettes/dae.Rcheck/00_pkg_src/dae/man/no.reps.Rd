\name{no.reps}
\alias{no.reps}
\title{Computes the number of replicates for an experiment}
\description{Computes the number of pure replicates required in an experiment
             to achieve a specified power.}
\usage{
no.reps(multiple=1., df.num=1.,
        df.denom=expression((df.num + 1.) * (r - 1.)), delta=1.,
        sigma=1., alpha=0.05, power=0.8, tol=0.1, print=FALSE)
}
\arguments{
 \item{multiple}{The multiplier, m, which when multiplied by the number of
                 pure replicates of a treatment, r, gives the number of
                 observations rm used in computing means for some, not
                 necessarily proper, subset of the treatment factors;
                 m is the replication arising from other treatment factors.
                 However, for single treatment factor experiments the
                 subset can only be the treatment factor and m = 1.}
 \item{df.num}{The degrees of freedom of the numerator of the F for testing
               the term involving the treatment factor subset.}
 \item{df.denom}{The degrees of freedom of the denominator of the F for
                 testing the term involving the treatment factor subset.}
 \item{delta}{The true difference between a pair of means for some, not
              necessarily proper, subset of the treatment factors.}
 \item{sigma}{The population standard deviation.}
 \item{alpha}{The significance level to be used.}
 \item{power}{The minimum power to be achieved.}
 \item{tol}{The maximum difference tolerated between the power required and
            the power computed in determining the number of replicates.}
 \item{print}{\code{TRUE} or \code{FALSE} to have or not have a table of power 
              calculation details printed out.}
}
\value{A list containing \code{nreps}, a single \code{numeric} value containing the computed number of pure replicates, and \code{power}, a single \code{numeric} value containing the power for the computed number of pure replicates.}
\author{Chris Brien}
\seealso{\code{\link{power.exp}}, \code{\link{detect.diff}} in package \pkg{dae}.}
\examples{
## Compute the number of replicates (blocks) required for a randomized 
## complete block design with four treatments. 
no.reps(multiple = 1, df.num = 3,
        df.denom = expression(df.num * (r - 1)), delta = 5,
	        sigma = sqrt(20), print = TRUE)
}
\keyword{design}