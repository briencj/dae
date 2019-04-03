\name{mat.I}
\alias{mat.I}
\title{Forms a unit matrix}
\description{Form the unit or identity \code{\link{matrix}} of order \code{order}.}
\usage{mat.I(order)}
\arguments{
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A square \code{\link{matrix}} whose diagonal elements are one and its off-diagonal 
       are zero.}
\author{Chris Brien}
\seealso{\code{\link{mat.J}}, \code{\link{mat.ar1}}}
\examples{
    col.I <- mat.I(order=4)
}
\keyword{array}
