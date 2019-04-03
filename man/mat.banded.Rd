\name{mat.banded}
\alias{mat.banded}
\title{Form a banded matrix from a vector of values}
\description{Takes the first value in \code{x}  and places it down the diagonal of the 
             \code{\link{matrix}}. Takes the second value in \code{x}  and places it 
             down the first subdiagonal, both below and above the diagonal of the 
             \code{\link{matrix}}.  The third value is placed in the second subdiagonal 
             and so on, until the bands for which there are elements in \code{x} have 
             been filled. All other elements in the \code{\link{matrix}} will be zero.}
\usage{mat.banded(x, nrow, ncol)}
\arguments{
 \item{x}{A \code{\link{numeric}} containing the values for each band from 1 to the 
          length of \code{x}.}
 \item{nrow}{The number of rows in the banded \code{\link{matrix}} being formed.}
 \item{ncol}{The number of columns in the banded \code{\link{matrix}} being formed.}
}
\value{An \eqn{nrow \times ncol}{nrow \times ncol} \code{\link{matrix}}.}
\author{Chris Brien}
\seealso{\code{matmult}, \code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, 
         \code{\link{mat.sar2}}, \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
         \code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}
         \code{\link{mat.I}}, \code{\link{mat.J}}}
\examples{
      m <- mat.banded(c(1,0.6,0.5), 5,5)
      m <- mat.banded(c(1,0.6,0.5), 3,4)
      m <- mat.banded(c(1,0.6,0.5), 4,3)

}
\keyword{array}
