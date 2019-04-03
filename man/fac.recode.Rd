\name{fac.recode}
\alias{fac.recode}
\title{Recodes factor \code{levels} using values in a vector. The values in the vector do 
             not have to be unique.}
\description{Recodes the \code{levels} and values of a factor using each value in the 
       \code{newlevels} vector  to replace the corresponding value in the vector of 
       \code{levels} of the \code{\link{factor}}.}
\usage{fac.recode(factor, newlevels, ...)}
\arguments{
 \item{factor}{The \code{\link{factor}} to be recoded.}
 \item{newlevels}{A \code{vector }of length \code{levels(factor)} containing values to use in 
                  the recoding.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}}.}
\author{Chris Brien}
\seealso{\code{\link{as.numfac}} and \code{\link{mpone}} in package \pkg{dae}, 
\code{\link{factor}}, \code{\link{relevel}}.}
\examples{
## set up a factor with labels
a <- factor(rep(1:4, 4), labels=c("A","B","C","D"))
 
## recode "A" and "D" to 1 and "B" and "C" to 2
b <- fac.recode(a, c(1,2,2,1), labels = c("a","b"))
}
\keyword{factor}
\keyword{manip}