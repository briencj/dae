\name{mat.sar2}
\alias{mat.sar2}
\title{Forms an sar2 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the sar2 pattern, a pattern used in crop competition 
             models. The resulting \code{\link{matrix}}
             is banded and is a constrained AR3 matrix. 
}
\usage{mat.sar2(gamma, order, print = NULL)}
\arguments{
 \item{gamma}{A \code{\link{numeric}} containing the two values of \code{gamma}, being 
              parameters linked with spatial dependence and competition.}
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
 \item{print}{A \code{\link{character}} giving the object to be printed. Currently, 
              only the claculated values of the \code{ar3parameters} can be printed. 
              If NULL, nothing is printed.}
}
\value{A banded correlation \code{\link{matrix}} whose elements follow an sar2 pattern.}
\details{The values of the AR3 parameters, phi, are calculated from the gammas as follows: \cr
         \code{phi[1] = gamma[1] + 2 * gamma[2]}; \code{phi[2] = -gamma[2] * (2*gamma[2] + gamma[1])}; \cr
         \code{phi[3] = gamma[1] * gamma[2] * gamma[2]}.

         Then the correlations in the correlation matrix, \code{corr} say, are calculated 
         from the correlation parameters, phi.
         Let \code{omega = 1 - phi[2] - phi[3] * (phi[1] + phi[3])}. 
             Then the values in 
             \itemize{
             \item the diagonal of \code{corr} (\code{k = 1}) are one; 
             \item the first subdiagonal band (\code{k = 2}) of \code{corr} are equal to 
                   \code{(phi[1] + phi[2]*phi[3]) / omega}; 
             \item the second subdiagonal band (\code{k = 3}) of \code{corr} are equal to \cr
                   \code{(phi[1] * (phi[1] + phi[3]) + phi[2] * (1 - phi[2])) / omega}; 
             \item the subsequent subdiagonal bands, (\code{k = 4:order}), of \code{corr} are equal to \cr
                      \code{phi[1]*corr[k-1] + phi[2]*corr[k-2] + phi[3]*corr[k-3]}.}}
\author{Chris Brien}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.cor}}, \code{\link{mat.corg}}, 
\code{\link{mat.banded}}, \code{\link{mat.exp}}, 
\code{\link{mat.gau}}, \code{\link{mat.ar1}}, \code{\link{mat.ar2}}, \code{\link{mat.ar3}}, \code{\link{mat.sar}}, 
\code{\link{mat.ma1}}, \code{\link{mat.ma2}}, \code{\link{mat.arma}}}
\examples{
    corr <- mat.sar2(gamma = c(-0.4, 0.2), order = 4)
    corr <- mat.sar2(gamma = c(-0.4, 0.2), order = 4, print = "ar3")
}
\keyword{array}
