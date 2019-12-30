\name{mat.Vpred}
\alias{mat.Vpred}
\title{Calculates the variances of a set of predicted effects from a mixed model}
\description{For \code{n} observations, \code{w} effects to be predicted, 
             \code{f} nuiscance fixed effects and \code{r} nuisance random effects, 
             the variances of a set of predicted effects is calculated using  
             the incidence matrix for the effects to be predicted and, optionally, 
             a variance matrix of the effects, an incidence matrix for the 
             nuisance fixed factors and covariates, the variance matrix of the nuisance 
             random effects in the mixed model and the residual variance matrix.
             
             This function has been superseded by \code{\link{mat.Vpredicts}}, which 
             allows the use of both matrices and \code{\link{formula}}e.}
\usage{
mat.Vpred(W, Gg = 0, X = matrix(1, nrow = nrow(W), ncol = 1), Vu = 0, R, eliminate)}
\arguments{
 \item{W}{The \code{n x w} incidence \code{\link{matrix}} for the \code{w} effects 
          to be predicted.}
 \item{Gg}{The \code{w x w} variance \code{\link{matrix}} of the \code{w} effects 
           to be predicted. 
           If the effects to be predicted are fixed, set to 0.}
 \item{X}{The \code{n x f} incidence \code{\link{matrix}} for the \code{f} nuisance 
          fixed factors and covariates. The default is a column vector of ones.}
 \item{Vu}{The \code{n x r} variance \code{\link{matrix}} of the \code{r} nuisance 
           random effects. If there are none, set to zero.}
 \item{R}{The residual variance \code{\link{matrix}}.}
 \item{eliminate}{The \code{n x n} \code{\link{projector}} onto the subspace 
          corresponding to the effects to be eliminated from the information matrix 
          prior to inverting it to form the variance \code{\link{matrix}} of the 
          predicted effects. It is only appropriate to use this option when the 
          effects to be predicted are fixed.}
}
\value{A \code{w x w} \code{\link{matrix}} containing the variances and covariances of the 
       predicted effects.}
\details{Firstly the information matrix is calculated as \cr
         \code{A <- t(W) \%*\% Vinv \%*\% W + ginv(Gg) - A\%*\%ginv(t(X)\%*\%Vinv\%*\%X)\%*\%t(A)}, 
         where \code{Vinv <- ginv(Vu + R)}, \code{A = t(W) \%*\% Vinv \%*\% X} and ginv(B) is the unique Moore-Penrose inverse of B formed using the eigendecomposition of B.
         
         If \code{eliminate} is set and the effects to be predicted are fixed then the reduced information matrix is calculated as \code{A <- (I - eliminate) Vinv (I - eliminate)}.
         
         Finally, the variance of the predicted effects is calculated: \code{Vpred <- ginv(A)}.}
\author{Chris Brien}
\seealso{\code{\link{designAmeasures}}, \code{\link{mat.Vpredicts}}.}
\references{Smith, A. B., D. G. Butler, C. R. Cavanagh and B. R. Cullis (2015). 
            Multi-phase variety trials using both composite and individual replicate 
            samples: a model-based design approach. 
            \emph{Journal of Agricultural Science}, \bold{153}, 1017-1029.}
\examples{
## Reduced example from Smith et al. (2015)
## Generate two-phase design
mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                            levels = c("Y","W","G","M","D","E"))
start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
rownames(start.design) <- NULL

## Set up matrices
n <- nrow(start.design)
W <- model.matrix(~ -1+ Variety, start.design)
ng <- ncol(W)
Gg<- diag(1, ng)
Vu <- with(start.design, fac.vcmat(Mrep, 0.3) + 
                         fac.vcmat(fac.combine(list(Mrep, Mday)), 0.2) + 
                         fac.vcmat(Frep, 0.1) + 
                         fac.vcmat(fac.combine(list(Frep, Fplot)), 0.2))
R <- diag(1, n)
  
## Calculate the variance matrix of the predicted random Variety effects
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
designAmeasures(Vp)

## Calculate the variance matrix of the predicted fixed Variety effects, 
## elminating the grand mean
Vp.reduc <- mat.Vpred(W = W, Gg = 0, Vu = Vu, R = R, 
                      eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
designAmeasures(Vp.reduc)
}
\keyword{array}
\keyword{design}
