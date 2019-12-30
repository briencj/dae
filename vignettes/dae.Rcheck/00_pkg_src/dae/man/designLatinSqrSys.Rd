\name{designLatinSqrSys}
\alias{designLatinSqrSys}
\title{Generate a systematic plan for a Latin Square design}
\description{Generates a systematic plan for a Latin Square design using the method of cycling the integers 1 to the number of treatments. The start of the cycle for each row, or the first column, can be specified as a vector of integers.}
\usage{designLatinSqrSys(order, start = NULL)}
\arguments{
 \item{order}{The number of treatments.}
 \item{start}{A \code{\link{numeric}} containing \code{order} unique values between one and \code{order}. These are interpreted as the value for the fist column for each row. If \code{NULL}, \code{1:order} is used.}
}
\value{A  \code{\link{numeric}} containing \code{order} x \code{order} integers between 1 and \code{order} such that, when the \code{\link{numeric}} is considered as a square matrix of size \code{order}, each integer occurs once and only once in each row and column of the matrix.}
\seealso{\code{\link{designRandomize}}, \code{\link{designPlot}}, \code{\link{designAnatomy}} in package \pkg{dae}.}
\examples{
   matrix(designLatinSqrSys(5, start = c(seq(1, 5, 2), seq(2, 5, 2))), nrow=5)
   designLatinSqrSys(3)
}
\keyword{array}
\keyword{design}
