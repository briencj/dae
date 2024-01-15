\name{mat.ma2}
\alias{mat.ma2}
\title{Forms an ma2 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ma2 pattern. The resulting \code{\link{matrix}} 
             is banded.}
\usage{mat.ma2(MAparameters, order)}
\arguments{
 \item{MAparameters}{A \code{\link{numeric}} containing the two moving average parameter values 
            of the process, being the weights given to the lag 1 and lag 2 random pertubations.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an ma2 pattern.}
\details{The correlations in the correlation matrix, \code{corr} say, are calculated 
         from the moving average parameters, \code{MAparameters}.
             The values in 
             \itemize{
             \item the diagonal (\code{k = 1}) of \code{corr} are one; 
             \item the first subdiagonal band (\code{k = 2}) of \code{corr} 
                   are equal to \cr 
                   \code{-MAparameters[1]*(1 - MAparameters[2]) / div}; 
             \item the second subdiagonal bande (\code{k = 3}) of \code{corr} 
                   are equal to \code{-MAparameters[2] / div}; 
             \item in subsequent disgonal bands, (\code{k = 4:order}), 
             of \code{corr} are zero,}
         where \code{div = 1 + MMAparameters[1]*MAparameter[1] + MAparameters[2]*MAparameters[2]}.}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, 
\code{\link{mat.exp}}, \code{\link{mat.gau}}, 
\code{\link{mat.banded}}, \code{\link{mat.ar1}}, \code{\link{mat.ar3}}, \code{\link{mat.sar2}}, 
\code{\link{mat.ma1}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.ma2(MAparameters = c(0.4, -0.2), order = 4)
}
\keyword{array}
