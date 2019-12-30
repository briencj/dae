\name{yates.effects}
\alias{yates.effects}
\title{Extract Yates effects}
\description{Extracts Yates effects from an \code{aov} object or \code{aovlist} object.}
\usage{yates.effects(aov.obj, error.term="Within", data=NULL)}
\arguments{
 \item{aov.obj}{An \code{aov} object or \code{aovlist} object created from a call to 
\code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function from which the Yates 
                   effects are estimated. Only required when \code{Error} used 
                   in call to \code{aov}.}
 \item{data}{A \code{data.frame} in which the variables specified in the \code{aov.obj}
          will be found. If missing, the variables are searched for in
          the standard way.}
}
\value{A \code{vector} of the Yates effects.}
\details{Yates effects are specific to \eqn{2^k}{2^k} experiments, where Yates
      effects are conventionally defined as the difference between the upper
      and lower levels of a factor. We follow the convention used in
      Box, Hunter and Hunter (1978) for scaling of higher order interactions:
      all the Yates effects are on the same scale, and represent the average
      difference due to the interaction between two different levels.
      Effects are estimated only from the error term supplied to the
      \code{error.term} argument.}
\author{Chris Brien}
\seealso{\code{\link{qqyeffects}} in package \pkg{dae}, \code{\link{aov}}.}
\examples{
## analysis of 2^4 factorial experiment from Table 10.6 of Box, Hunter and 
## Hunter (1978) Statistics for Experimenters. New York, Wiley.
## use ?Fac4Proc.dat for data set details
data(Fac4Proc.dat)
Fac4Proc.aov <- aov(Conv ~ Catal * Temp * Press * Conc + Error(Runs),
                                                            Fac4Proc.dat)
round(yates.effects(Fac4Proc.aov, error.term="Runs", data=Fac4Proc.dat), 2)
}
\keyword{design}
\keyword{htest}