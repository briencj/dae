\name{fac.split}
\alias{fac.split}
\title{Splits a factor whose levels consist of several delimited strings into several factors}
\description{Splits a \code{\link{factor}}, whose \code{levels} consist of strings delimited by a 
             separator character, into several \code{\link{factor}s}. It uses the function 
             \code{\link{strsplit}}, with \code{fixed = TRUE} to split the \code{levels}.}
\usage{fac.split(combined.factor, factor.names, sep=",", ...)}
\arguments{
 \item{combined.factor}{A \code{\link{factor}} to be split into several \code{\link{factor}s}.}
 \item{factor.names}{A \code{\link{list}} of names for factors and associated \code{levels}, if required.
   The names of the componeonts of the \code{\link{list}} are used for the names of the new factors. 
   Each component of the \code{\link{list}} should either be \code{NULL} or a vector of levels for 
   the new \code{\link{factor}}. If a component is \code{NULL} then the unique values  for the 
   supplied factor are used as the \code{levels}, which are sorted into alphabetical order. 
   If a either a \code{\link{numeric}} or a \code{\link{character}} vector is supplied for a 
   component, then these are supplied as the \code{levels} of the new factor.}
 \item{sep}{A \code{character} string that separates the levels in the \code{combined.factor}.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating 
            the new \code{\link{factor}}.}
}
\value{A \code{\link{data.frame}} containing the new \code{\link{factor}s}.}
\author{Chris Brien}
\seealso{\code{\link{fac.divide}}, \code{\link{fac.uncombine}}, \code{\link{fac.combine}} in package \pkg{dae} and 
          \code{\link{strsplit}}.}
\examples{
## Forma combined factor to split
data(Oats.dat)
tmp <- within(Oats.dat, Trts <- fac.combine(list(Variety, Nitrogen), combine.levels = TRUE))

##Variety levels sorted into alphabetical order
trts.dat <- fac.split(combined.factor = tmp$Trts, 
                      factor.names = list(Variety = NULL, Nitrogen = NULL))

##Variety levels order from Oats.dat retained
trts.dat <- fac.split(combined.factor = tmp$Trts, 
                      factor.names = list(Variety = levels(tmp$Variety), Nitrogen = NULL))
}
\keyword{factor}
\keyword{manip}