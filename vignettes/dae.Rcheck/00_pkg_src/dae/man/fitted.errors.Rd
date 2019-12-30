\name{fitted.errors}
\alias{fitted.errors}
\title{Extract the fitted values for a fitted model}
\description{An alias for the generic function \code{\link{fitted}}. When it is 
     available, the method \code{\link{fitted.aovlist}} extracts the fitted values, which is provided 
     in the \pkg{dae} package to cover \code{aovlist} objects.}
\usage{\method{fitted}{errors}(object, error.term=NULL, ...)}
\arguments{
 \item{object}{An \code{aovlist} object created from a call to \code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function down to which effects 
      are extracted for adding to the fitted values. The order of terms is as 
      given in the ANOVA table. If \code{error.term} is \code{NULL} effects are 
      extracted from all \code{Error} terms.}
 \item{...}{Further arguments passed to or from other methods.}
}
\value{A numeric vector of fitted values.}
\section{Warning}{See \code{\link{fitted.aovlist}} for specific information about fitted 
      values when an \code{Error} function is used in the call to the 
      \code{\link{aov}} function.}
\author{Chris Brien}
\seealso{\code{\link{fitted.aovlist}}, \code{\link{resid.errors}}, \code{\link{tukey.1df}} 
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

## three equivalent ways of extracting the fitted values
fit  <- fitted.aovlist(RCBDPen.aov)
fit <- fitted(RCBDPen.aov, error.term = "Blend:Flask")
fit <- fitted.errors(RCBDPen.aov, error.term = "Blend:Flask")
}
\keyword{models}
\keyword{htest}