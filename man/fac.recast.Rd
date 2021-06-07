\name{fac.recast}
\alias{fac.recast}
\title{Recasts a factor by modifying the values in the factor vector and/or the levels attribute, possibly combining some levels into a single level.}
\description{A \code{\link{factor}} is comprised of a vector of values and a \code{\link{levels}} attribute. 
       This function can modify these separately or jointly. The \code{newlevels} argument recasts 
       both the values of a \code{\link{factor}} vector and the \code{\link{levels}} attribute, using each 
       value in the \code{newlevels} vector to replace the corresponding value in both \code{\link{factor}} 
       vector and the \code{\link{levels}} attribute. The \code{\link{factor}}, possibly 
       with the new levels, can have its \code{levels} attribute reordered and/or new 
       \code{labels} associated with the \code{levels} using the \code{levels.order} and \code{newlabels} 
       arguments.}
\usage{fac.recast(factor, newlevels = NULL, levels.order = NULL, newlabels = NULL, ...)}
\arguments{
 \item{factor}{The \code{\link{factor}} to be recast.}
 \item{newlevels}{A \code{vector }of length \code{levels(factor)} that changes both the values in the 
                  \code{factor} vector and its \code{\link{levels}} attribute. The values in the 
                  \code{newlevels} vector need not be unique, but there must be as many values as 
                  there are levels in the supplied \code{factor}. The levels in the vector of the 
                  supplied \code{factor} that have the same value in \code{newlevels} will be combined 
                  in the recast \code{\link{factor}}. The values in the new \code{\link{levels}} 
                  attribute can be re-oredered using \code{levels.order}.}
 \item{levels.order}{A \code{\link{vector}} that specifies the order of the levels in the 
                  \code{\link{levels}} attribute of the recast \code{\link{factor}}. 
                  If \code{newlevels} is \code{NULL}, must be of 
                  length \code{levels(factor)} and contain the old levels in the new order for 
                  the recast \code{\link{factor}}. If \code{newlevels} is not \code{NULL}, the vector 
                  must be of \code{length(unique(newlevels))} and contain the unique values in 
                  \code{newlevels} in the new order for the recast \code{\link{factor}}. 
                  The values in the \code{\link{factor}} vector whose \code{\link{levels}} are being 
                  re-ordered will be unchanged. If \code{levels.order} is \code{NULL}, then the 
                  current \code{\link{levels}} attribute of \code{factor} is used.}
 \item{newlabels}{A \code{\link{vector}} of length \code{levels(factor)} if \code{newlevels} is 
                  \code{NULL}, and of length \cr\code{unique(newlevels)} if it is not \code{NULL}. 
                  It should contain the values to be used as \code{labels} in the recast \code{\link{factor}}. 
                  Effectively, this changes the values in the \code{\link{factor}} vector to those 
                  given in \code{newlabels} and the \code{levels} attribute to \code{newlabels}.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}}.}
\author{Chris Brien}
\seealso{\code{\link{fac.uselogical}, \link{as.numfac}} and \code{\link{mpone}} in package \pkg{dae}, 
\code{\link{factor}}, \code{\link{relevel}}.}
\examples{
## set up a factor with labels
Treats <- factor(rep(1:4, 4), labels=letters[1:4])
 
## recast to reduce the levels: "a" and "d" to 1 and "b" and "c" to 2, i.e. from 4 to 2 levels
A <- fac.recast(Treats, newlevels = c(1,2,2,1), labels = letters[1:2])
A <- fac.recast(Treats, newlevels = letters[c(1,2,2,1)])

#reduce the levels from 4 to 2, with re-ordering the levels vector without changing the values 
#of the new recast factor vector
A <- fac.recast(Treats, newlevels = letters[c(1,2,2,1)], levels.order = letters[2:1])  

#reassign the values in the factor vector without re-ordering the levels attribute
A <- fac.recast(Treats, newlevels = letters[4:1])  

#reassign the values in the factor vector, with re-ordering the levels attribute
A <- fac.recast(Treats, newlabels = letters[4:1])

#reorder the levels attribute with changing the values in the factor vector
A <- fac.recast(Treats, levels.order = letters[4:1])  

#reorder the values in the factor vector without changing the levels attribute
A <- fac.recast(Treats, newlevels = 4:1, newlabels = levels(Treats))

}
\keyword{factor}
\keyword{manip}