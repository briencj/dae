\name{mat.random}
\alias{mat.random}
\title{Calculates the variance matrix for the random effects from a mixed model, based  on a supplied formula or a matrix.}
\description{For \code{n} observations, compute the variance matrix of the random effects. 
             The \code{\link{matrix}} can be specified using a \code{\link{formula}} 
             for the \code{random} effects and a \code{\link{list}} of values of the  
             variance components for the terms specified in the \code{random} \code{\link{formula}}. 
             If a  \code{\link{matrix}} specifying the variances of 
             the nuisance random effects is supplied then it is returned as the value 
             from the function.}
\usage{
mat.random(random, G, design, keep.order = TRUE)}
\arguments{
 \item{random}{A \code{\link{formula}} or a \code{\link{matrix}}. If a 
          \code{\link{formula}}, it specifies the random effects from which 
          the \code{\link{matrix}} for the contribution of the random effects 
          to the variance matrix can be generated. If it is a \code{\link{matrix}}, 
          it must be an \code{n x n} \code{\link{matrix}} and will be passed 
          through as the required variance matrix for the random effects. 
          The default is 0, which implies that there are no random effects.}
 \item{G}{This term only needs to be set if \code{random} is a 
           \code{\link{formula}}. Then it is set to a \code{\link{list}}, in which 
           each component is either a single value or a \code{\link{matrix}}; 
           there needs to be a component for each term in the expanded 
           \code{\link{formula}}, with the order of the terms and components 
           matching. If it is a single value, a diagonal matrix of dimension 
           equal to the product of the numbers of levels of the factors in 
           its term. If it is a matrix, its dimension must be equal to the 
           product of the numbers of levels of the factors in its term.}
 \item{design}{A \code{\link{data.frame}} containing the design to be used in an 
          experiment and for which the variane matrix is required. It is not required  
          when the only \code{\link{formula}} specified is an intercept-only 
          \code{\link{formula}}.}
 \item{keep.order}{A \code{\link{logical}} indicating whether the terms should keep 
                   their position in the expanded \code{formula} projector, or 
                   reordered so that main effects precede two-factor 
                   interactions, which precede three-factor interactions and 
                   so on.}
}
\value{A \code{n x n} \code{\link{matrix}} containing the variance matrix for the random effects.}
\details{If \eqn{\bold{Z}_i}{Zi} is the is incidence matrix for the \code{random} nuisance effects 
         in  \eqn{\bold{u}_i}{ui} for a term in \code{random} and \eqn{\bold{u}_i}{ui} has 
         variance matrix  \eqn{\bold{G}_i}{Gi} so that the contribution of the random effectst to 
         the variance matrix for \eqn{\bold{Y}}{Y} is  
        \eqn{\bold{V}_u =  \Sigma (\bold{Z}_i\bold{G}_i(\bold{Z}_i)^T)}{Vu = sum(Zi \%*\% Gi\%*\% t(Zi))}.   
}
\author{Chris Brien}
\seealso{\code{\link{mat.Vpredicts}}.}
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
Vu <- with(start.design, fac.vcmat(Mrep, gammas["Mrep"]) + 
                         fac.vcmat(fac.combine(list(Mrep,Mday)), gammas["Mrep:Mday"]) + 
                         fac.vcmat(Frep, gammas["Frep"]) + 
                         fac.vcmat(fac.combine(list(Frep,Fplot)), gammas["Frep:Fplot"]))

## Calculate the variance matrix of the predicted random Variety effects using formulae
Vu <- mat.random(random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                 G = as.list(gammas[c(4,5,2,3)]), 
                 design = start.design)
}
\keyword{array}
\keyword{design}
