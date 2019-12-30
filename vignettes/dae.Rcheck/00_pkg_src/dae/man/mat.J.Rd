\name{mat.J}
\alias{mat.J}
\title{Forms a square matrix of ones}
\description{Form the square \code{\link{matrix}} of ones of order \code{order}.}
\usage{mat.J(order)}
\arguments{
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A square \code{\link{matrix}} all of whose elements are one.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.ar1}}}
\examples{
    col.J <- mat.J(order=4)
}
\keyword{array}
