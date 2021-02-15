\name{print.pstructure}
\alias{print.pstructure}
\title{Prints a pstructure.object}
\description{Prints a \code{\link{pstructure.object}}, which is of class \code{pstructure}. 
       The df, terms and sources are coerced into a \code{data.frame} and printed; 
       the marginality matrix is printed separately.}
\usage{\method{print}{pstructure}(x, which = "all", ...)}
\arguments{
 \item{x}{The \code{\link{pstructure.object}}, which is of class \code{pstructure} and is to be printed.}
 \item{which}{A character \code{vector} nominating the components of the 
              \code{\link{pstructure.object}} to print. Must be \code{all} 
              or some combination of \code{projectors}, \code{marginality}, 
              and \code{aliasing}.}
 \item{...}{Further arguments passed to or from other methods.}
}
\author{Chris Brien}
\seealso{\code{\link{print}}, \code{\link{print.default}}, \code{\link{show}}.}
\examples{
## Generate a data.frame with 4 factors, each with three levels, in standard order
ABCD.lay <- fac.gen(list(A = 3, B = 3, C = 3, D = 3))

## create a pstructure object based on the formula ((A*B)/C)*D
ABCD.struct <- pstructure.formula(~ ((A*B)/C)*D, data =ABCD.lay)

## print the object either using the Method function, the generic function or show
print.pstructure(ABCD.struct)
print(ABCD.struct)
ABCD.struct
}
\keyword{array}
\keyword{projector}
