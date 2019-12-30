\name{residuals.aovlist}
\alias{residuals.aovlist}
\alias{residuals}
\title{Extract the residuals from an aovlist object}
\description{Extracts the residuals from \code{error.term} or, if \code{error.term} 
             is not specified, the last \code{error.term} in the analysis. It is a 
             method for the generic function \code{\link{residuals}}.}
\usage{\method{residuals}{aovlist}(object, error.term=NULL, ...)}
\arguments{
 \item{object}{An \code{aovlist} object created from a call to \code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function for which the 
                   residuals are to be extracted. If \code{error.term} is 
                   \code{NULL} the residuals are extracted from the last 
                   \code{Error} term.}
 \item{...}{Further arguments passed to or from other methods.}
}
\value{A numeric \code{vector} containing the residuals.}
\author{Chris Brien}
\seealso{\code{\link{fitted.errors}}, \code{\link{resid.errors}}, 
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
}
\keyword{methods}
\keyword{models}
\keyword{htest}