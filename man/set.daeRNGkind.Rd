\name{set.daeRNGkind}
\alias{set.daeRNGkind}
\title{Sets the values of daeRNGkind for the package dae in the daeEnv environment}
\description{A function that sets the \code{character} value \code{daeRNGkind} that 
             specifies the \code{kind} of Random Number generator to use in \code{dae}.
             The value is stored in a \code{character} named \code{daeRNGkind} in the 
             \code{daeEnv} environment. It is initially set to "Mersenne-Twister" and 
             can be changed using \code{\link{get.daeRNGkind}}. For details of the 
             different Random Number Generators available in \code{R}, see the \code{R} 
             help for \code{RNGkind}.}
\usage{set.daeRNGkind(kind = "Mersenne-Twister")}
\arguments{
 \item{kind}{A \code{character} to which \code{daeRNGkind} is to be set.}
}
\value{The value of \code{daeRNGkind} is returned invisibly.}
\author{Chris Brien}
\seealso{\code{\link{get.daeRNGkind}}.}
\examples{
## set daeRNGkind to L'Ecuyer-CMRG.
set.daeRNGkind("L'Ecuyer-CMRG")
}
\keyword{manip}
