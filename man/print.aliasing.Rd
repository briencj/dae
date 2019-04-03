\name{print.aliasing}
\alias{print.aliasing}
\title{Print an aliasing data.frame}
\description{Prints an aliasing \code{\link{data.frame}}.}
\usage{\method{print}{aliasing}(x, which.criteria = c("aefficiency","eefficiency","order"),  ...)}
\arguments{
 \item{x}{The \code{\link{data.frame}} that is also of class \code{aliasing} and is to be printed.}
 \item{which.criteria}{A character \code{vector} nominating the efficiency 
                       criteria to be included in the summary of aliasing 
                       between terms. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency}, 
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}. 
                       If \code{none}, no criteria are printed.}
 \item{...}{Further arguments passed to the \code{print} method for \code{data.frame}.}
}
\author{Chris Brien}
\seealso{\code{\link{print}}, \code{\link{print.default}}, \code{\link{show}}.}
\examples{
## Generate a data.frame with 3 factors length 12
pseudo.lay <- data.frame(pl = factor(1:12),
                         ab = factor(rep(1:4, times=3)),
                         a = factor(rep(1:2, times=6)))


## create a pstructure object
trt.struct <- pstructure(~ ab+a, data = pseudo.lay)

## print the object either using the Method function, the generic function or show
print.aliasing(trt.struct$aliasing)
print(trt.struct$aliasing, which.criteria = "none")
trt.struct$aliasing
}
\keyword{array}
\keyword{projector}
