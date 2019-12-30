\name{fitted.aovlist}
\alias{fitted.aovlist}
\alias{fitted}
\title{Extract the fitted values for a fitted model from an aovlist object}
\description{Extracts the fitted values as the sum of the effects 
     for all the fitted terms in the model, stopping at \code{error.term} 
     if this is specified. It is a method for the generic function 
     \code{\link{fitted}}.}
\usage{\method{fitted}{aovlist}(object, error.term=NULL, ...)}
\arguments{
 \item{object}{An \code{aovlist} object created from a call to \code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function down to which effects 
      are extracted for adding to the fitted values. The order of terms is as 
      given in the ANOVA table. If \code{error.term} is \code{NULL} effects are 
      extracted from all \code{Error} terms.}
 \item{...}{Further arguments passed to or from other methods.}
}
\value{A numeric \code{vector} of fitted values.}
\note{Fitted values will be the sum of effects for terms from the model, but only 
      for terms external to any \code{Error} function. If you want effects for 
      terms in the \code{Error} function to be included, put them both inside 
      and outside the \code{Error} function so they are occur twice.}
\author{Chris Brien}
\seealso{\code{\link{fitted.errors}}, \code{\link{resid.errors}}, \code{\link{tukey.1df}} 
in package \pkg{dae}.}
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

## two equivalent ways of extracting the fitted values
fit  <- fitted.aovlist(RCBDPen.aov)
fit <- fitted(RCBDPen.aov, error.term = "Blend:Flask")
}
\keyword{methods}
\keyword{models}
\keyword{htest}