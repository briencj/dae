\name{fac.layout}
\alias{fac.layout}
\title{Generate a randomized layout for an experiment}
\description{Warning: \code{fac.layout} will be deprecated in future versions, its 
             synonym \code{\link{designRandomize}} being preferred. 

             Given a systematic design, generate a layout for an experiment 
             consisting of randomized \code{\link{factor}s} that are randomized 
             to the unrandomized \code{\link{factor}s}, taking into account the 
             nesting between the unrandomized factors for the design.}
\usage{fac.layout(unrandomized, nested.factors=NULL, except=NULL, 
                  randomized, seed=NULL, unit.permutation=TRUE)}
\arguments{
 \item{unrandomized}{A \code{\link{data.frame}} or a \code{\link{list}} of 
                     \code{\link{factor}s}, along with their \code{levels}. If 
                     a \code{\link{list}}, the name of each component of the 
                     \code{\link{list}} is a \code{\link{factor}} name and the 
                     component is either a single \code{numeric} value that 
                     is the number of \code{levels}, a numeric \code{vector} 
                     that contains the \code{levels} of the \code{\link{factor}} 
                     or a character code{vector} that contains the \code{labels} 
                     of the \code{levels} of the \code{\link{factor}}.}
 \item{nested.factors}{A \code{\link{list}} of the unrandomized 
                       \code{\link{factor}s} 
                       that are nested in other \code{\link{factor}s} in 
                       \code{unrandomized}. The name of each component is the 
                       name of a \code{\link{factor}} that is nested and the 
                       component is a character \code{vector} containing the 
                       \code{\link{factor}s} within which it is nested. 
                       It is emphasized that the nesting is a property of the 
                       design that is being employed (it is only partly based 
                       on the intrinsic nesting.}
 \item{except}{A \code{\link{character} vector} containing the names 
               unrandomized \code{\link{factor}s} that are to be excepted 
               from the randomization.}
 \item{randomized}{A \code{\link{factor}} or a \code{\link{data.frame}} 
                   containing the values of the \code{\link{factor}(s)} to be 
                   randomized.}
 \item{seed}{A single value, interpreted as an integer, that specifies the 
                       starting value of the random number generator.}
 \item{unit.permutation}{A \code{logical} indicating whether to include the 
                    \code{.Unit} and \code{.Permutation} columns in the 
                    \code{\link{data.frame}}.}
}
\value{A \code{\link{data.frame}} consisting of the values for \code{.Units} 
       and \code{.Permutation} \code{vectors}, provided \code{unit.permutation} 
       is \code{TRUE}, along with the values for the unrandomized and randomized 
       \code{\link{factor}s} that specify the randomized layout for the 
       experiment.}
\details{This \code{functon} uses the method of randomization described by 
         Bailey (1981). That is, a permutation of the units that respects the 
         nesting for the design, but does not permute any of the factors in the 
         \code{except vector}, is obtained. This permutation is 
         applied to the unrandomized \code{\link{factor}s} and then a 
         \code{\link{data.frame}} containing both the permuted unrandomized and 
         unpermuted randomized \code{\link{factor}s} is formed. To produce the 
         randomized layout, the rows of the joint \code{\link{data.frame}} are 
         reordered so that its unrandomized \code{\link{factor}s} are in either 
         standard order or, if a \code{\link{data.frame}} was suppled to 
         \code{unrandomized}, data frame order. 

         The \code{.Units} and \code{.Permutation} \code{vectors} enable one to 
         swap between this permutation and the randomized layout. 
         The \emph{i}th value in \code{.Permutation} gives the unit to which 
         unit \emph{i} was assigned in the randomization. 
}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}} in package \pkg{dae}.}
\examples{
## generate a randomized layout for a 4 x 4 Latin square
## (the nested.factors agument is not needed here as none of the 
## factors are nested)
LS.unit <- data.frame(row = ordered(rep(c("I","II","III","IV"), times=4)), 
                      col = factor(rep(c(0,2,4,6), each=4)))
LS.ran <- data.frame(treat = factor(c(1:4, 2,3,4,1, 3,4,1,2, 4,1,2,3)))
data.frame(LS.unit, LS.ran)
LS.lay <- fac.layout(unrandomized=LS.unit, randomized=LS.ran, seed=7197132) 
LS.lay[LS.lay$.Permutation,]

## generate a randomized layout for a replicated randomized complete 
## block design, with the block factors arranged in standard order for 
## rep then plot and then block
RCBD.unit <- list(rep = 2, plot=1:3, block = c("I","II"))
## specify that plot is nested in block and rep and that block is nested
## in rep
RCBD.nest <- list(plot = c("block","rep"), block="rep")
## generate treatment factor in systematic order so that they correspond
## to plot
tr <- factor(rep(1:3, each=2, times=2))
## obtain randomized layout
RCBD.lay <- fac.layout(unrandomized=RCBD.unit, 
                       nested.factors=RCBD.nest, 
                       randomized=tr, seed=9719532)
#sort into the original standard order
RCBD.perm <- RCBD.lay[RCBD.lay$.Permutation,]
#resort into randomized order
RCBD.lay <- RCBD.perm[order(RCBD.perm$.Units),]

## generate a layout for a split-unit experiment in which: 
## - the main-unit factor is A with 4 levels arranged in 
##   a randomized complete block design with 2 blocks;
## - the split-unit factor is B with 3 levels.
SPL.lay <- fac.layout(unrandomized=list(block = 2, main.unit = 4, split.unit = 3), 
                       nested.factors=list(main.unit = "block", 
                                           split.unit = c("block", "main.unit")), 
                       randomized=fac.gen(list(A = 4, B = 3), times = 2), 
                       seed=155251978, unit.permutation=FALSE)
}
\keyword{design}
\keyword{factor}
\keyword{datagen}
\references{Bailey, R.A. (1981) A unified approach to design of experiments. 
            \emph{Journal of the Royal Statistical Society, Series A}, 
            \bold{144}, 214--223.}