\name{tukey.1df}
\alias{tukey.1df}
\title{Performs Tukey's one-degree-of-freedom-test-for-nonadditivity}
\description{Performs Tukey's one-degree-of-freedom-test-for-nonadditivity 
             on a set of residuals from an analysis of variance.}
\usage{tukey.1df(aov.obj, data, error.term="Within")}
\arguments{
 \item{aov.obj}{An \code{aov} object or \code{aovlist} object created from a call to 
\code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function whose residuals are 
        to be tested for nonadditivity. Only required when the \code{Error} 
        function used in call to \code{aov}, so that an  \code{aovlist} object 
        is created.}
 \item{data}{A \code{data.frame} containing the original response variable and 
       factors used in the call to \code{\link{aov}}.}
}
\value{A \code{list} containing Tukey.SS, Tukey.F, Tukey.p, Devn.SSq being the SSq 
        for the 1df test, F value for test and the p-value for the test.}
\note{In computing the test quantities fitted values must be obtained. 
      If \code{error.term} is specified, fitted values will be the sum of 
      effects extracted from terms from the \code{Error} function, but only down 
      to that specified by \code{error.term}.The order of terms is as given in the 
      ANOVA table. If \code{error.term} is unspecified, all effects for terms 
      external to any \code{Error} terms are extracted and summed.
      
      Extracted effects will only be for terms external to any \code{Error} function. 
      If you want effects for terms in the \code{Error} function to be included, 
      put them both inside and outside the \code{Error} function so they are 
      occur twice.}
\author{Chris Brien}
\seealso{\code{\link{fitted.errors}}, \code{\link{resid.errors}} in package \pkg{dae}.}
\examples{
## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## Obtain the quantities for Tukey's test
tukey.1df(RCBDPen.aov, RCBDPen.dat, error.term = "Blend:Flask")
}
\keyword{models}
\keyword{htest}