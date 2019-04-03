\name{designAmeasures}
\alias{designAmeasures}
\title{Calculates the average variance of pairwise differences from the variance matrix 
       for predictions that can be obtained using mat.Vpredicts}
\description{Calculates the average variance of pairwise differences between 
             (of elementary contrasts of) preductions, possibly for different subgroups 
             of the predictions, from the variance matrix for the predictions. If groups 
             are specified then the A-optimality measures are calculated for the differences 
             between predictions within each group and for those between predictions from 
             different groups. If groupsizes are specified, but groups are not, the
             predictions will be sequentially broken into groups of the size specified by 
             the elements of groupsizes. The groups can be named.
}
\usage{designAmeasures(Vpred, groupsizes = NULL, groups = NULL)}
\arguments{
 \item{Vpred}{The variance \code{\link{matrix}} of the predictions.}
 \item{groupsizes}{A \code{\link{numeric}} containing group sizes. The sum of the 
                   elements of \code{groupsizes} must be less than or equal to the 
                   order of \code{Vpred}. If \code{groupsizes} is a named 
                   \code{vector}, the names are used to label the groups. 
                   If \code{NULL}, either \code{groups} is used or 
                   the average for all pairwise differences is obtained.}
 \item{groups}{A \code{\link{list}}, each element of which specifies a subroup of the 
               predictions over whose pairwise differences the variances are 
               to be averaged. If there is more than one group, the variances of 
               all between and within group pairwise differences are averaged. 
               If the elements of \code{groups} are named , the names are used to 
               label the groups. If \code{groups} is \code{NULL}, either 
               \code{groupsizes} is used or the average for all pairwise differences 
               is obtained.}
}
\value{A \code{\link{matrix}} containing the within and between group A-optimality measures.}
\author{Chris Brien}
\seealso{\code{\link{mat.Vpred}}, \code{\link{designAnatomy}}.}
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
  
## Calculate the varaince matrix of the predicted random Variety effects
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
  
## Calculate A-optimality measure
designAmeasures(Vp)
designAmeasures(Vp, groups=list(fldUndup = c(1:4), fldDup = c(5,6)))
grpsizes <- c(4,2)
names(grpsizes) <- c("fldUndup", "fldDup")
designAmeasures(Vp, groupsizes = grpsizes)
designAmeasures(Vp, groupsizes = c(4))
designAmeasures(Vp, groups=list(c(1,4),c(5,6)))

## Calculate the variance matrix of the predicted fixed Variety effects, elminating the grand mean
Vp.reduc <- mat.Vpred(W = W, Gg = 0, Vu = Vu, R = R, 
                      eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
## Calculate A-optimality measure
designAmeasures(Vp.reduc)


}
\keyword{design}
