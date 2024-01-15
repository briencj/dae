\name{mat.ar2}
\alias{mat.ar2}
\title{Forms an ar2 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ar2 pattern. The resulting \code{\link{matrix}} 
             is banded.}
\usage{mat.ar2(ARparameters, order)}
\arguments{
 \item{ARparameters}{A \code{\link{numeric}} containing the two autoregressive parameter values 
            of the process, being the weights given to the lag 1 and lag 2 response values.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an ar2 pattern.}
\details{The correlations in the correlation matrix, \code{corr} say, are calculated 
         from the autoregressive parameters, \code{ARparameters}.
             The values in 
             \itemize{
             \item the diagonal (\code{k = 1}) of \code{corr} are one; 
             \item the first subdiagonal band (\code{k = 2}) of \code{corr} 
                   are equal to \cr 
                   \code{ARparameters[1]/(1-ARparameters[2])}; 
             \item in subsequent disgonal bands, (\code{k = 3:order}), 
             of \code{corr} are \cr 
             \code{ARparameters[1]*corr[k-1] + ARparameters[2]*corr[k-2]}.}}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, 
\code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.banded}}, \code{\link{mat.ar1}}, \code{\link{mat.ar3}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.ar2(ARparameters = c(0.4, 0.2), order = 4)
}
\keyword{array}
