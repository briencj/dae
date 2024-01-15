\name{mat.ar3}
\alias{mat.ar3}
\title{Forms an ar3 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ar3 pattern. The resulting \code{\link{matrix}} 
             is banded.}
\usage{mat.ar3(ARparameters, order)}
\arguments{
 \item{ARparameters}{A \code{\link{numeric}} containing the three autoregressive parameter values 
            of the process, being the weights given to the lag 1, lag 2 and lag 3 
            response values.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an ar3 pattern.}
\details{The correlations in the correlation matrix, \code{corr} say, are calculated 
         from the autoregressive parameters, \code{ARparameters}. \cr
         Let \code{omega = 1 - ARparameters[2] - ARparameters[3] * (ARparameters[1] + ARparameters[3])}. 
             Then the values in 
             \itemize{
             \item the diagonal of \code{corr} (\code{k = 1}) are one; 
             \item the first subdiagonal band (\code{k = 2}) of \code{corr} are equal to \cr 
                   \code{(ARparameters[1] + ARparameters[2]*ARparameters[3]) / omega}; 
             \item the second subdiagonal band (\code{k = 3}) of \code{corr} are equal to \cr
                   \code{(ARparameters[1] * (ARparameters[1] + ARparameters[3]) + } \cr 
                                 \code{ARparameters[2] * (1 - ARparameters[2])) / omega}; 
             \item the subsequent subdiagonal bands, (\code{k = 4:order}), of \code{corr} are equal to \cr
                      \code{ARparameters[1]*corr[k-1] + ARparameters[2]*corr[k-2] + ARparameters[3]*corr[k-3]}.}}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, \code{\link{mat.banded}}, 
\code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.ar3(ARparameters = c(0.4, 0.2, 0.1), order = 4)
}
\keyword{array}
