\name{fac.nested}
\alias{fac.nested}
\title{creates a factor, the nested factor, whose values are generated within those of the factor 
   nesting.fac}
\description{Creates a nested \code{\link{factor}} whose \code{levels} are generated 
  within those of the factor \code{nesting.fac}. All elements of \code{nesting.fac} 
  having the same level are numbered from 1 to the number of different elements 
  having that level.}
\usage{fac.nested(nesting.fac, nested.levs=NA, nested.labs=NA, ...)}
\arguments{
 \item{nesting.fac}{The \code{\link{factor}} within each of whose \code{levels} the 
       created \code{\link{factor}} is to be generated.}
 \item{nested.levs}{Optional \code{vector} of levels for the \code{\link{factor}}. Any data 
       value that does not match a value in \code{levels} will be \code{NA} in 
       the \code{\link{factor}}. The default value of \code{nested.levs} is the  list of 
       numbers from 1 to the maximum replication of the \code{levels} of 
       \code{nesting.fac}, represented as characters.}
 \item{nested.labs}{Optional \code{vector} of values to use as \code{labels} for the 
       \code{levels} of the \code{\link{factor}}. The default is \code{as.character(nested.levs)}.} 
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}} that is a character \code{vector} with class attribute 
       "\code{\link[=factor-class]{factor}}" and a \code{levels} attribute which 
       determines what character strings may be included in the \code{vector}. It has 
       a different level for of the values of the nesting.fac with the same level.}
\note{The levels of \code{nesting.fac} do not have to be equally replicated.}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}}, \code{\link{fac.multinested}} in package \pkg{dae}, \code{\link{factor}}.}
\examples{
## set up factor A
A <- factor(c(1, 1, 1, 2, 2))

## create nested factor
B <- fac.nested(A)
}
\keyword{factor}
\keyword{manip}