\name{fac.uncombine}
\alias{fac.uncombine}
\title{Cleaves a single factor, each of whose levels has delimited strings, into several 
       factors using the separated strings.}
\description{Cleaves a single \code{\link{factor}} into several factors whose \code{levels},  
the levels of the original \code{\link{factor}} consisting of several delimited strings that 
can be separated to form the levels of the \code{new.factors}. That is, it reverses the process 
of combining factors that \code{\link{fac.combine}} performs.}
\usage{fac.uncombine(factor, new.factors, sep=",", ...)}
\arguments{
 \item{factor}{A \code{\link{factor}} or \code{\link{character}} that has values that are strings 
   deleimited by the delimiter specified by \code{sep}.}
 \item{new.factors}{A \code{\link{list}}, whose component names are the names of the 
   new \code{\link{factor}s} to be formed. If a component is not \code{NULL}, 
   then they are used as the levels of the corresponding \code{factor}.}
 \item{sep}{A \code{character} string that separates the levels of the \code{new.factors} in the 
   levels \code{factor}.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating 
            the new \code{\link{factor}}.}
}
\value{A \code{\link{data.frame}} whose columns consist of the \code{\link{factor}s} listed in 
  \code{new.factors} and whose values have been computed from the values of the combined 
  \code{\link{factor}}.}
\author{Chris Brien}
\seealso{\code{\link{fac.split}}, \code{\link{fac.combine}}, \code{\link{fac.divide}} in package \pkg{dae} and 
          \code{\link{strsplit}}.}
\examples{
## set up two factors and combine them
facs <- fac.gen(list(A = letters[1:3], B = 1:2), each = 4)
facs$AB <- with(facs, fac.combine(list(A, B), combine.levels = TRUE))

## now reverse the proces and uncombine the two factors
new.facs <- fac.uncombine(factor = facs$AB, 
                          new.factors = list(A = letters[1:3], B = NULL), 
                          sep = ",")
new.facs <- fac.uncombine(factor = facs$AB, 
                          new.factors = list(A = NULL, B = NULL), 
                          sep = ",")
}
\keyword{factor}
\keyword{manip}