\name{resid.errors}
\alias{resid.errors}
\title{Extract the residuals for a fitted model}
\description{An alias for the generic function \code{\link{residuals}}. When it is 
     available, the method \code{\link{residuals.aovlist}} extracts residuals, which is provided 
     in the package \pkg{dae} to cover \code{aovlist} objects.}
\usage{resid.errors(...)}
\arguments{
 \item{...}{Arguments passed to \code{\link{residuals.aovlist}}.}
}
\value{A numeric \code{vector} containing the residuals.}
\note{See \code{\link{residuals.aovlist}} for specific information about the 
      residuals when an \code{Error} function is used in the call to the 
      \code{\link{aov}} function.}
\author{Chris Brien}
\seealso{\code{\link{fitted.errors}}, \code{\link{residuals.aovlist}}, 
\code{\link{tukey.1df}} in package \pkg{dae}.}
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

## two equivalent ways of extracting the residuals
res  <- residuals.aovlist(RCBDPen.aov)
res <- residuals(RCBDPen.aov, error.term = "Blend:Flask")
res <- resid.errors(RCBDPen.aov)
}
\keyword{models}
\keyword{htest}