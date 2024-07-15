\name{as.numfac}
\alias{as.numfac}
\title{Convert a factor to a numeric vector, possibly centering or scaling the values}
\description{Converts a \code{\link{factor}} to a numeric \code{vector} with approximately the 
             numeric values of its \code{levels}. Hence, the \code{levels} of the 
             \code{\link{factor}} must be numeric values, stored as characters. It uses the method 
             described in \code{\link{factor}}. Use \code{\link{as.numeric}} to convert a 
             \code{\link{factor}} to a numeric \code{vector} with integers 1, 2, ... corresponding 
             to the positions in the list of levels. The numeric values can be centred and/or scaled. 
             You can also use \code{\link{fac.recast}} to recode the levels to numeric 
             values. If a \code{numeric} is supplied and both \code{center} and \code{scale} are
             \code{FALSE}, it is left unchanged; otherwise, it will be centred and scaled 
             according to the settings of \code{center} and \code{scale}.}
\usage{as.numfac(factor, center = FALSE, scale = FALSE)}
\arguments{
 \item{factor}{The \code{\link{factor}} to be converted.}
 \item{center}{Either a \code{\link{logical}} value or a \code{\link{numeric}}-alike value, 
               where \emph{numeric-alike} means that \code{\link{as.numeric}} will be applied 
               successfully if \code{\link{is.numeric}} is not \code{TRUE}.}
 \item{scale}{Either a \code{\link{logical}} value, a \code{\link{numeric}}-alike vector,   
               where \emph{numeric-alike} means that \code{\link{as.numeric}} will be applied 
               successfully if \code{\link{is.numeric}} is not \code{TRUE}.}
}
\details{The value of \code{center} specifies how the centring is done. If it is \code{TRUE},
         the mean of the unique values of the \code{\link{factor}} are subtracted, after the 
         \code{\link{factor}} is converted to a \code{\link{numeric}}. If \code{center} is 
         \code{\link{numeric}}, it is subtracted from \code{\link{factor}}'s converted 
         \code{\link{numeric}} values.  If \code{center} is \code{FALSE} no scaling is done.

          The value of \code{scale} specifies how scaling is performed, after any centring is 
          done. If \code{scale} is \code{TRUE}, the unique converted values of the 
          \code{\link{factor}} are divided by (i) the standard deviaton when the values have 
          been centred and (ii) the root-mean-square error if they have not been centred; 
          the root-mean-square is given by \eqn{\sqrt{\Sigma(x^2)/(n-1)}}{sqrt(sum(x^2)/(n-1))}, 
          where \code{x} contains the unique converted \code{\link{factor}} values and \emph{n} is the number 
          of unique values. If \code{scale} is \code{FALSE} no scaling is done.
}
\value{A numeric \code{vector}. An \code{NA} will be stored for any value of the factor whose 
       level is not a number.}
\author{Chris Brien}
\seealso{\code{\link{as.numeric}}, \code{\link{fac.recast}} in package \pkg{dae}, 
         \code{\link{factor}}, \code{\link{scale}}.}
\examples{
## set up a factor and convert it to a numeric vector
a <- factor(rep(1:6, 4))
x <- as.numfac(a)
x <- as.numfac(a, center = TRUE)
x <- as.numfac(a, center = TRUE, scale = TRUE)
}
\keyword{factor}
\keyword{manip}