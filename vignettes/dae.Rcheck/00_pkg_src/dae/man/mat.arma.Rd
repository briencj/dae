\name{mat.arma}
\alias{mat.arma}
\title{Forms an arma correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the arma pattern. The resulting \code{\link{matrix}} 
             is banded.}
\usage{mat.arma(ARparameter, MAparameter, order)}
\arguments{
 \item{ARparameter}{A \code{\link{numeric}} value for the autoregressive parameter
            of the process, being the weight given to the lag 1 response values.}
 \item{MAparameter}{A \code{\link{numeric}} value for the moving average parameter
            of the process, being the weight given to the lag 1 random variable.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an arma pattern.}
\details{The correlations in the correlation matrix, \code{corr} say, are calculated 
         from the correlation parameters, \code{ARparameters}.
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
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.banded}}, \code{\link{mat.ar1}}, \code{\link{mat.ar3}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}}
\examples{
    corr <- mat.arma(ARparameter = 0.4, MAparameter = -0.2, order = 4)
}
\keyword{array}
