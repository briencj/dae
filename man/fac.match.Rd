\name{fac.match}
\alias{fac.match}
\title{Match, for each combination of a set of columns in \code{x}, the row that has the same combination in \code{table}}
\description{Match, for each combination of a set of columns in \code{x}, 
             the rows that has the same combination in \code{table}. 
             The argument \code{multiples.allow} controls what happens when there are 
             multple matches in \code{table} of a combination in \code{x}. }
\usage{fac.match(x, table, col.names, nomatch = NA_integer_, multiples.allow = FALSE)}
\arguments{
 \item{x}{an R object, normally a \code{data.frame}, possibly a 
          \code{matrix}.}
 \item{table}{an R object, normally a \code{data.frame}, possibly a 
       \code{matrix}.}
 \item{col.names}{A \code{character vector} giving the columns in \code{x} 
                  and \code{table} that are to be matched.}
 \item{nomatch}{The value to be returned in the case when no match is found. 
                Note that it is coerced to integer.}
 \item{multiples.allow}{A \code{logical} indicating whether multiple matches of a 
                        combination in \code{x} to those in \code{table}  is allowed. 
                        If \code{multiples.allow} is \code{FALSE}, an error is generated. 
                        If \code{multiples.allow} is \code{TRUE}, the first occuence in 
                        \code{table} is matched. This function can be viewed as a 
                        generalization to multiple vectors of the \code{\link{match}} 
                        function that applies to single vectors.}
}
\value{A \code{\link{vector}} of length equal to \code{x} that gives the 
       rows in \code{table} that match the combinations of 
       \code{col.names} in \code{x}. The order of the rows is the same as 
       the order of the combintions in \code{x}. The value returned if a combination is
       unmatched is specified in the \code{nomatch} argument.}
\author{Chris Brien}
\seealso{\code{\link{match}}}
\examples{\dontrun{
#A single unmatched combination
kdata <- data.frame(Expt="D197-5", 
                    Row=8, 
                    Column=20, stringsAsFactors=FALSE)
index <- fac.match(kdata, D197.dat, c("Expt", "Row", "Column"))

# A matched and an unmatched combination
kdata <- data.frame(Expt=c("D197-5", "D197-4"), 
                    Row=c(8, 10), 
                    Column=c(20, 8), stringsAsFactors=FALSE)
index <- fac.match(kdata, D197.dat, c("Expt", "Row", "Column"))
}}
\keyword{design}
\keyword{factor}
