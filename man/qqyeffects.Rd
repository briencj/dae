\name{qqyeffects}
\alias{qqyeffects}
\title{Half or full normal plot of Yates effects}
\description{Produces a half or full normal plot of the Yates effects from a
       \eqn{2^k}{2^k} factorial experiment.}
\usage{
qqyeffects(aov.obj, error.term="Within", data=NULL, pch=16, 
           full=FALSE, ...)}
\arguments{
 \item{aov.obj}{An \code{aov} object or \code{aovlist}object created from a call to 
\code{\link{aov}}.}
 \item{error.term}{The term from the \code{Error} function from which the Yates 
                   effects are estimated. Only required when \code{Error} used 
                   in call to \code{aov}.}
 \item{data}{A \code{data.frame} in which the variables specified in the \code{aov.obj}
          will be found. If missing, the variables are searched for in
          the standard way.}
 \item{pch}{The number of a plotting symbol to be drawn when plotting points
 (use \code{help(points)} for details).}
 \item{full}{whether a full or half normal plot is to be produced. The
       default is for a half-normal plot; \code{full=TRUE} produces a full normal
       plot.}
 \item{\dots}{Further graphical parameters may be specified (use
       \code{help(par)} for possibilities.}
}
\value{Returns, invisibly, a list with components x and y, giving coordinates
       of the plotted points.}
\details{A half or full normal plot of the Yates effects is produced.
      You will be able to interactively select effects to be labelled (click 
      reasonably close to the point and on the side where you want the label placed).
      \bold{Right click on the graph and select Stop when you have finished
      labelling effects.} A regression line fitted to the unselected effects
      and constrained to go through the origin is plotted. Also, a list of the
      labelled effects, if any, are printed to standard ouptut.}
\author{Chris Brien}
\seealso{\code{\link{yates.effects}} in package \pkg{dae}, \code{\link{qqnorm}}.}
\examples{
## analysis of 2^4 factorial experiment from Table 10.6 of Box, Hunter and 
## Hunter (1978) Statistics for Experimenters. New York, Wiley.
## use ?Fac4Proc.dat for data set details
data(Fac4Proc.dat)
Fac4Proc.aov <- aov(Conv ~ Catal * Temp * Press * Conc + Error(Runs),
                                                            Fac4Proc.dat)
qqyeffects(Fac4Proc.aov, error.term="Runs", data=Fac4Proc.dat)
}
\keyword{iplot}
\keyword{hplot}
\keyword{design}
\keyword{htest}