\name{mat.Vpredicts}
\alias{mat.Vpredicts}
\title{Calculates the variances of a set of predicted effects from a mixed model, based  on supplied matrices or formulae.}
\description{For \code{n} observations, \code{w} effects to be predicted, 
             \code{f} nuiscance fixed effects, \code{r} nuisance random effects and 
             \code{n} residuals, 
             the variances of a set of predicted effects is calculated using the   
             incidence matrix for the effects to be predicted and, optionally, 
             a variance matrix of these effects, an incidence matrix for the 
             nuisance fixed factors and covariates, the variance matrix of the 
             nuisance random effects and the residual variance matrix. 
             The \code{\link[=matrix]{matrices}} can be supplied directly or 
             using \code{\link{formula}}e 
             and a  \code{\link{matrix}} specifying the variances of 
             the nuisance random effects. The difference between 
             \code{mat.Vpredicts} and \code{\link{mat.Vpred}} is that the 
             former has different names for equivalent arguments and the latter 
             does not allow for the use of \code{\link{formula}}e.}
\usage{
mat.Vpredicts(target, Gt = 0, fixed = ~ 1, random, G, R, design, 
              eliminate, keep.order = TRUE, result = "variance.matrix")}
\arguments{
 \item{target}{The \code{n x w} incidence \code{\link{matrix}} for the \code{w} 
          effects targetted for prediction, or a \code{\link{formula}} from which 
          the \code{\link{matrix}} can be generated.}
 \item{Gt}{The value of the variance component for the targetted effects or the 
           \code{w x w} variance \code{\link{matrix}} of the \code{w} targetted 
           effects. If the targetted effects are fixed, set \code{Gt} to 0.}
 \item{fixed}{The \code{n x f} incidence \code{\link{matrix}} for the \code{f} 
          nuisance fixed effects and covariates, or a \code{\link{formula}} from 
          which the \code{\link{matrix}} can be generated. The default is a 
          \code{\link{formula}} for an intercept-only model.}
 \item{random}{A \code{\link{formula}} or a \code{\link{matrix}}. If a 
          \code{\link{formula}}, it specifies the random effects from which 
          the \code{\link{matrix}} for the contribution of the random effects 
          to the variance matrix can be generated. If it is a matrix, it must be an 
          \code{n x n} \code{\link{matrix}} and will be passed on to form the variance 
          matrix of the observations. The default is 0, which 
          implies that there are no random effects.}
 \item{G}{This term only needs to be set if \code{random} is a 
           \code{\link{formula}}. Then it is set to a \code{\link{list}}, in which 
           each component is either a single value or a \code{\link{matrix}}; 
           there needs to be a component for each term in the expanded 
           \code{\link{formula}}, with the order of the terms and components 
           matching. If it is a single value, a diagonal matrix of dimension 
           equal to the product of the numbers of levels of the factors in 
           its term. If it is a matrix, its dimension must be equal to the 
           product of the numbers of levels of the factors in its term.}
 \item{R}{The \code{n x n} residual variance \code{\link{matrix}}. If \code{R} is 
          not set in the call, then it is set to the identity \code{\link{matrix}}.}
 \item{design}{A \code{\link{data.frame}} containing the design to be used in an 
          experiment from which predictions are to be obtained. It is not required  
          when the only \code{\link{formula}} specified is an intercept-only 
          \code{\link{formula}}.}
 \item{eliminate}{The \code{n x n} \code{\link{projector}} onto the subspace 
          corresponding to the effects to be eliminated from the information matrix 
          prior to inverting it to form the variance \code{\link{matrix}} of the 
          predicted effects. It is only appropriate to use this option when the 
          effects to be predicted are fixed.}
 \item{keep.order}{A \code{\link{logical}} indicating whether the terms should keep 
                   their position in the expanded \code{formula} projector, or 
                   reordered so that main effects precede two-factor 
                   interactions, which precede three-factor interactions and 
                   so on.}
 \item{result}{A \code{\link{character}} indicating which matrix is to be returned: \code{variance.matrix} or \code{information.matrix}.}
}
\value{A \code{w x w} \code{\link{matrix}} containing the variances and covariances of the 
       predicted effects or the information matrix for the effects, depending on the setting of \code{result}. The matrix has its rank as an attribute.}
\details{The mixed model for which the predictions are to be obtained is of the form 
         \eqn{\bold{Y} = \bold{X\beta} + \bold{Ww} + \bold{Zu} + \bold{e}}{Y = X\beta + Ww + Zu + e}, 
        where   \eqn{\bold{W}}{W} is the incidence matrix for the \code{target} predicted 
        effects  \eqn{\bold{w}}{w},  \eqn{\bold{X}}{X} is the is incidence matrix for the 
        \code{fixed} nuisance effects \eqn{\bold{\beta}}{\beta},  \eqn{\bold{Z}}{Z} is the 
        is incidence matrix for the \code{random} nuisance effects  \eqn{\bold{u}}{u}, 
         \eqn{\bold{e}}{e} are the residuals; the  \eqn{\bold{u}}{u} are assumed 
        to have variance matrix  \eqn{\bold{G}}{G} so that their contribution to the 
        variance matrix for \eqn{\bold{Y}}{Y} is  
        \eqn{\bold{Vu} = \bold{ZGZ}^T}{Vu = Z \%*\% G \%*\% t(Z)} and \eqn{\bold{e}}{e} 
        is assumed to have variance matrix \eqn{\bold{R}}{R}. 
        If the \code{target} effects are random then the variance matrix for 
        \eqn{\bold{w}}{w} is \eqn{\bold{G}_t}{Gt} so that their 
        contribution to the variance matrix for \eqn{\bold{Y}}{Y} is 
        \eqn{\bold{WG}_t\bold{W}^T}{W \%*\% Gt \%*\% t(W)}.   
         
         As described in Hooks et al. (2009, Equation 19), the information matrix is 
         calculated as \cr
         \code{A <- t(W) \%*\% Vinv \%*\% W + ginv(Gg) - A\%*\%ginv(t(X)\%*\%Vinv\%*\%X)\%*\%t(A)}, 
         where \code{Vinv <- ginv(Vu + R)}, \code{A = t(W) \%*\% Vinv \%*\% X} 
         and ginv(B) is the unique Moore-Penrose inverse of B formed using the 
         eigendecomposition of B.
         
         Then, if \code{eliminate} is set and the effects to be predicted are 
         fixed then the reduced information matrix is calculated as 
         \code{A <- (I - eliminate) Vinv (I - eliminate)}.
         
         Finally, if \code{result} is set to \code{variance.matrix},  the variance of the predicted effects is calculated: 
         \code{Vpred <- ginv(A)} and returned; otherwise the information matrix A is returned. The rank of the matrix to be returned is obtain via a singular value decomposition of the information matrix, it being the number of nonzero eigenvalues. An eigenvalue is regarded as zero if it is less than 
\code{daeTolerance}, which is initially set to.Machine$double.eps ^ 0.5 (about 1.5E-08). The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.}
\author{Chris Brien}
\seealso{\code{\link{designAmeasures}}, \code{\link{mat.random}}, \code{\link{mat.Vpred}}.}
\references{Hooks, T., Marx, D., Kachman, S., and Pedersen, J. (2009). Optimality 
            criteria for models with random effects. \emph{Revista Colombiana de 
            Estadistica}, \bold{32}, 17-31. 
   
            Smith, A. B., D. G. Butler, C. R. Cavanagh and B. R. Cullis (2015). 
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

## Set gammas
terms <- c("Variety", "Frep", "Frep:Fplot", "Mrep", "Mrep:Mday", "Mrep:Mday:Mord")
gammas <- c(1, 0.1, 0.2, 0.3, 0.2, 1)
names(gammas) <- terms

## Specify matrices to calculate the variance matrix of the predicted fixed Variety effects 
W <- model.matrix(~ -1 + Variety, start.design)
Vu <- with(start.design, fac.vcmat(Mrep, gammas["Mrep"]) + 
                         fac.vcmat(fac.combine(list(Mrep,Mday)), gammas["Mrep:Mday"]) + 
                         fac.vcmat(Frep, gammas["Frep"]) + 
                         fac.vcmat(fac.combine(list(Frep,Fplot)), gammas["Frep:Fplot"]))
R <- diag(1, nrow(start.design))
  
## Calculate variance matrix
Vp <- mat.Vpredicts(target = W, random=Vu, R=R, design = start.design)

## Calculate the variance matrix of the predicted random Variety effects using formulae
Vp <- mat.Vpredicts(target = ~ -1 + Variety, Gt = 1, 
                    fixed = ~ 1, 
                    random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                    G = as.list(gammas[c(4,5,2,3)]), 
                    R = R, design = start.design)
designAmeasures(Vp)

## Calculate the variance matrix of the predicted fixed Variety effects, 
## elminating the grand mean
n <- nrow(start.design)
Vp.reduc <- mat.Vpredicts(target = ~ -1 + Variety, 
                          random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                          G = as.list(gammas[c(4,5,2,3)]), 
                          eliminate = projector(matrix(1, nrow = n, ncol = n)/n), 
                          design = start.design)
designAmeasures(Vp.reduc)
}
\keyword{array}
\keyword{design}
