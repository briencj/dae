\name{rep.data.frame}
\alias{rep.data.frame}
\title{Replicate the rows of a data.frame by repeating each row consecutively and/or repeating all rows as a group}
\description{Replicate the rows of a \code{data.frame} by repeating each row consecutively and/or repeating all rows as a group.}
\usage{
\method{rep}{data.frame}(x, times=1, each=1, ...)}
\arguments{
 \item{x}{A \code{\link{data.frame}} whose rows are to be repeated.}
 \item{times}{The number of times to repeat the whole set of rows, after the rows have been 
              replicated consecutively \code{each} times.}
 \item{each}{The number of times to replicate consecutively each row in the \code{data.frame}.}
 \item{...}{Further arguments passed to or from other methods. Unused at present.}
}
\value{A \code{\link{data.frame}} with replicated rows.}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}} in package \pkg{dae} and \code{rep}}
\examples{
rep(fac.gen(list(a = 2, b = 2)), times=2, each=2)
}
\keyword{datagen}
