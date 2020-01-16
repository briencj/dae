\name{fac.uselogical}
\alias{fac.uselogical}
\title{Forms a two-level \code{\link{factor}} from a \code{\link{logical}} object.}
\description{Forms a two-level \code{\link{factor}} from a \code{\link{logical}} object. 
       It can be used to recode a \code{\link{factor}} when the resulting \code{\link{factor}} 
       is to have only two \code{levels}.}
\usage{fac.uselogical(x, levels = c(TRUE, FALSE), labels = c("yes", "no"), ...)}
\arguments{
 \item{x}{A \code{\link{logical}} vector with values \code{TRUE} or \code{FALSE}. If the vector is 
                not a \code{logical}, \code{\link{as.logical}} will be used in an attempt to 
                coerce it to \code{logical}.}
 \item{levels}{A \code{vector} of length two with values \code{TRUE} or \code{T} and 
               \code{FALSE} or \code{F}, in either order depending on which of \code{TRUE} 
               or \code{FALSE} is to be the first level.}
 \item{labels}{A \code{vector} of length two with values to be used as labels for the first 
               and second levels, respectively.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}}.}
\author{Chris Brien}
\seealso{\code{\link{fac.recode}}, \code{\link{as.numfac}} and \code{\link{mpone}} in package \pkg{dae}, 
\code{\link{factor}}, \code{\link{relevel}}.}
\examples{
## set up a factor with labels
Treats <- factor(rep(1:4, 4), labels=c("A","B","C","D"))
 
## recode "A" and "D" to "a" and "B" and "C" to "b"
B <- fac.uselogical(Treats \%in\% c("A", "D"), labels = c("a","b"))
B <- fac.uselogical(Treats \%in\% c("A", "D"), labels = c(-1,1))

## suppose level A in factor a is a control treatment
## set up a factor Control to discriminate between control and treated
Control <- fac.uselogical(Treats == "A")
}
\keyword{factor}
\keyword{manip}