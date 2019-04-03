\name{fac.combine}
\alias{fac.combine}
\title{Combines several factors into one}
\description{Combines several \code{\link{factor}s} into one whose \code{levels} are the 
combinations of the used \code{levels} of the individual \code{\link{factor}s}.}
\usage{fac.combine(factors, order="standard", combine.levels=FALSE, sep=",", ...)}
\arguments{
 \item{factors}{A \code{\link{list}} of \code{\link{factor}s} all of the same length.}
 \item{order}{Either \code{standard} or \code{yates}. The order in which the 
   \code{levels} combinations of the \code{\link{factor}s} are to be considered as 
   numbered when forming the \code{levels} of the combined \code{\link{factor}}; 
   \code{standard} numbers them as if they are arranged in 
   standard order, that is with the levels of the first factor moving slowest and 
   those of the last factor moving fastest; \code{yates} numbers them as if they 
   are arranged in Yates order, that is with the levels of the first factor moving 
   fastest and those of the last factor moving slowest.}
 \item{combine.levels}{A \code{logical} specifying whether the \code{levels} 
   labels of the new \code{\link{factor}} are to be combined from those of the 
   \code{\link{factor}s} being combined. The default is to 
   use the integers from 1 to the product of the numbers of combinations 
   of used \code{levels} of the individual \code{\link{factor}s}, numbering the levels 
   according to \code{order}.}
 \item{sep}{A \code{character} string to separate the levels when 
\code{combine.levels = TRUE}.}
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating 
            the new \code{\link{factor}}.}
}
\value{A \code{\link{factor}} whose \code{levels} are formed form the observed 
       combinations of the \code{levels} of the individual \code{\link{factor}s}.}
\author{Chris Brien}
\seealso{\code{\link{fac.divide}} in package \pkg{dae}.}
\examples{
## set up two factors
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## obtain six-level factor corresponding to the combinations of A and B
AB <- fac.combine(list(A,B))
}
\keyword{factor}
\keyword{manip}