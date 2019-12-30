\name{fac.gen}
\alias{fac.gen}
\title{Generate all combinations of several factors and, optionally, replicate them}
\description{Generate all combinations of several factors and, optionally, replicate them.}
\usage{fac.gen(generate, each=1, times=1, order="standard")}
\arguments{
 \item{generate}{A \code{\link{list}} of named objects and numbers that specify the 
                     \code{\link{factor}s} whose \code{levels} are to be generated and 
                     the pattern in these \code{levels}.If a component of the 
                     \code{\link{list}} is named, then the component 
                     should be either a) a single \code{numeric} value that is the 
                     number of \code{levels}, b) a numeric \code{vector} that contains the 
                     \code{levels} of the \code{\link{factor}}, or c) a character \code{vector} that 
                     contains the \code{labels} of the \code{levels} of the \code{\link{factor}}.}
 \item{each}{The number of times to replicate consecutively the elements of the 
                     \code{levels} generated according to pattern specified by the 
                     \code{generate} argument.}
 \item{times}{The number of times to repeat the whole generated pattern of 
                     \code{levels} generated according to pattern specified by the 
                     \code{generate} argument.}
 \item{order}{Either \code{standard} or \code{yates}. The order in which the 
                     speed of cycling through the \code{levels} is to move;
                     combinations of the \code{\link{factor}s} are to be considered as numbered; 
                     \code{standard} cycles through the levels of the first 
                     factor slowest and the last factor moving fastest; 
                     \code{yates} cycles through the levels of the first factor 
                     fastest and last factor moving slowest.}
}
\value{A \code{\link{data.frame}} of generated \code{levels} with columns 
       corresponding to the code{factors} in the \code{generate} list.}
\details{The \code{levels} of each \code{\link{factor}} are generated in a hierarchical 
      pattern where the \code{levels} of one \code{\link{factor}} are held constant 
      while those of the adjacent \code{\link{factor}} are cycled through the complete 
      set once. If a number is supplied instead of a name, the pattern is 
      generated as if a \code{\link{factor}} with that number of \code{levels} had been 
      supplied in the same position as the number. However, no \code{levels} are 
      stored for this unamed \code{\link{factor}}.}
\section{Warning}{Avoid using factor names F and T as these might be confused
                  with FALSE and TRUE.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}} in package \pkg{dae}}
\examples{
## generate a 2^3 factorial experiment with levels - and +, and 
## in Yates order
mp <- c("-", "+")
fnames <- list(Catal = mp, Temp = mp, Press = mp, Conc = mp)
Fac4Proc.Treats <- fac.gen(generate = fnames, order="yates")

## Generate the factors A, B and D. The basic pattern has 4 repetitions
## of the levels of D for each A and B combination and 3 repetitions of 
## the pattern of the B and D combinations for each level of A. This basic 
## pattern has each combination repeated twice, and the whole of this 
## is repeated twice. It generates 864 A, B and D combinations.
gen <- list(A = 3, 3, B = c(0,100,200), 4, D = c("0","1"))
fac.gen(gen, times=2, each=2)
}
\keyword{design}
\keyword{factor}
\keyword{datagen}
