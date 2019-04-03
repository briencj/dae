\name{as.data.frame.pstructure}
\alias{as.data.frame.pstructure}
\title{Coerces a pstructure.object to a data.frame.}
\description{Coerces a \code{\link{pstructure.object}}, which is of class \code{pstructure}, 
       to a \code{\link{data.frame}}. One can choose whether or not to include the \code{marginality} 
       matrix in the data.frame. The \code{aliasing} component is excluded.}
\usage{
\method{as.data.frame}{pstructure}(x, row.names = NULL, optional = FALSE, ..., 
              omit.marginality = FALSE)}
\arguments{
   \item{x}{The \code{\link{pstructure.object}}, which is of class \code{pstructure} and is to be coerced.}
 \item{row.names}{NULL or a \code{\link{character}} vector giving the row names for 
          the data frame. Missing values are not allowed.}
 \item{optional}{A \code{\link{logical}} passed to \code{as.data.frame}. If TRUE, 
           setting row names and converting column names (to syntactic names: 
           see \code{make.names)} is optional. Note that all of R's base 
           package \code{as.data.frame()} methods use \code{optional} only for 
           column names treatment, basically with the meaning of 

           \code{data.frame(*, check.names = !optional)}.}
 \item{...}{Further arguments passed to or from other methods.}
 \item{omit.marginality}{A \code{\link{logical}}, which, if \code{TRUE}, 
           results in the \code{marginality} matrix being omitted from the 
           \code{data.frame}.}
}
\value{A \code{\link{data.frame}} with as many rows as there are non-aliased terms 
       in the \code{\link{pstructure.object}}. The columns are \code{df}, \code{terms}, 
       \code{sources} and, if \code{omit.marginality} is \code{FALSE}, the columns of 
       the generated \code{levels} with columns of the \code{marginality} matrix 
       that is stored in the \code{marginality} component of the object.}
\author{Chris Brien}
\seealso{\code{\link{as.data.frame}}.}
\examples{
## Generate a data.frame with 4 factors, each with three levels, in standard order
ABCD.lay <- fac.gen(list(A = 3, B = 3, C = 3, D = 3))

## create a pstructure object based on the formula ((A*B)/C)*D
ABCD.struct <- pstructure.formula(~ ((A*B)/C)*D, data =ABCD.lay)

## print the object either using the Method function or the generic function 
ABCS.dat <- as.data.frame.pstructure(ABCD.struct)
as.data.frame(ABCD.struct)
}
\keyword{array}
\keyword{projector}
