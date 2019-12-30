\name{mat.sar}
\alias{mat.sar}
\title{Forms an sar correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the sar pattern. The resulting \code{\link{matrix}}
             is banded. 
}
\usage{mat.sar(SARparameter, order)}
\arguments{
 \item{SARparameter}{A \code{\link{numeric}} containing the single value of the 
              parameter from which the correlations are calculated.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an sar pattern.}
\details{The values of the correlations in the correlation matrix, \code{corr} say, are calculated 
         from the SARparameter, gamma as follows. The values in 
             \itemize{
             \item the diagonal of \code{corr} (\code{k = 1}) are one; 
             \item the first subdiagonal band (\code{k = 2}) of \code{corr} are equal to 
                   \code{gamma/(1 + (gamma * gamma / 4))}; 
             \item the subsequent subdiagonal bands, (\code{k = 3:order}), of \code{corr} are equal to \cr
                      \code{gamma * corr[k-1] - (gamma * gamma/4) * corr[k-2].}}
}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.banded}}, \code{\link{mat.exp}}, 
\code{\link{mat.gau}}, \code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, 
\code{\link{mat.sar2}}, \code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.sar(SARparameter = -0.4, order = 4)
}
\keyword{array}
