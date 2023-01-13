\name{designRandomize}
\alias{designRandomize}
\title{Randomize allocated to recipient factors to produce a layout for an experiment}
\description{A systematic design is specified by a set of 
             \code{allocated} \code{\link{factor}s} that have been assigned to a set of 
             \code{recipient} \code{\link{factor}s}. In textbook designs the 
             \code{allocated} \code{\link{factor}s} are the treatment factors and the 
             \code{recipient} \code{\link{factor}s} are the \code{\link{factor}s} 
             indexing the units. To obtain a randomized layout for a systematic design 
             it is necessary to provide (i) the systematic arrangement of the 
             \code{allocated} \code{\link{factor}s}, (ii) a \code{\link{list}} of the 
             \code{recipient} \code{\link{factor}s} or a \code{\link{data.frame}} with 
             their values, and (iii) the nesting of the 
             \code{recipient} \code{\link{factor}s} for the design being randomized. 
             Given this information, the \code{allocated} \code{\link{factor}s} 
             will be randomized to the \code{recipient} \code{\link{factor}s}, 
             taking into account the nesting between the \code{recipient} factors 
             for the design.
             However, \code{allocated} \code{\link{factor}s} that 
             have different values associated with those \code{recipient} 
             \code{\link{factor}s} that are in the \code{except} vector will remain 
             unchanged from the systematic design.
             
             Also, if \code{allocated} is \code{NULL} then a random permutation 
             of the \code{recipient} \code{\link{factor}s} is produced 
             that is consistent with their nesting as specified by 
             \code{nested.recipients}.

             For examples of its use also see the vignette accessed via 
             \code{vignette("DesignNotes", package="dae")} and for a discussion of 
             its use see Brien, Sermarini and Demetro (2023).
}
\usage{
designRandomize(allocated = NULL, recipient, nested.recipients = NULL, 
                except = NULL, seed = NULL, unit.permutation = FALSE, ...)}
\arguments{
 \item{allocated}{A \code{\link{factor}} or a \code{\link{data.frame}} 
                   containing the systematically allocated values of the 
                   \code{\link{factor}(s)}. If \code{NULL}, a random permutation 
                   of the \code{recipient} \code{\link{factor}s} is produced 
                   that is consistent with their nesting as specified by 
                   \code{nested.recipients}.}
 \item{recipient}{A \code{\link{data.frame}} or a \code{\link{list}} of 
                     \code{\link{factor}s}, along with their \code{levels} that 
                     specify the set of \code{recipient} \code{\link{factor}s} 
                     that are allocated levels of the 
                     \code{allocated} \code{\link{factor}s}. 

                     If a \code{\link{list}}, the name of each component of the 
                     \code{\link{list}} is a \code{\link{factor}} name and the 
                     component is either (i) a single \code{numeric} value that 
                     is the number of \code{levels}, (ii) a numeric \code{vector} 
                     that contains the \code{levels} of the \code{\link{factor}},  
                     (iii) or a \code{\link{character}} \code{vector} 
                     that contains the \code{labels} of the \code{levels} of the 
                     \code{\link{factor}}. The values of \code{\link{factor}s} 
                     will be generated in standard order using 
                     \code{\link{fac.gen}} and so the values in \code{allocated} 
                     must match this.}
 \item{nested.recipients}{A \code{\link{list}} of the \code{recipient} 
                       \code{\link{factor}s} 
                       that are nested in other \code{\link{factor}s} in 
                       \code{recipient}. The name of each component is the 
                       name of a \code{\link{factor}} that is nested and the 
                       component is a character \code{vector} containing the 
                       \code{\link{factor}s} within which it is nested. 
                       The randomization is controlled by \code{nested.recipients}: 
                       nested recipient factors are permuted within those 
                       factors that nest them. Only the nesting is specified: 
                       it is assumed that if two factors are not nested then 
                       they must be crossed.
                       It is emphasized that the nesting is a property of the 
                       design that is being employed (it is only partly based 
                       on the intrinsic or physical crossing and nesting).}
 \item{except}{A \code{\link{character} vector} containing the names of 
               \code{recipient} \code{\link{factor}s} that are to be excepted 
               from the permutation; any \code{allocated} \code{\link{factor}s} 
               whose values differ between the levels 
               of the \code{\link{factor}s} in this vector will not have those 
               values randomized.}
 \item{seed}{A single \code{\link{numeric}} value, interpreted as an integer, 
             that specifies the starting value of the random number generator.}
 \item{unit.permutation}{A \code{logical} indicating whether to include the 
                    \code{.Unit} and \code{.Permutation} columns in the 
                    \code{\link{data.frame}}.}
 \item{...}{Further arguments passed to or from other methods. Unused at present.}
}
\value{A \code{\link{data.frame}} with the values for the \code{recipient} and 
       \code{allocated} \code{\link{factor}s} that specify the layout for the 
       experiment and, if \code{unit.permutation} is \code{TRUE}, the values 
       for \code{.Units} and \code{.Permutation} \code{vectors}.}
\details{A systematic design is specified by the 
         matching of the supplied \code{allocated} and \code{recipient} 
         \code{\link{factor}s}. If \code{recipient} is a \code{\link{list}} 
         then \code{\link{fac.gen}} is used to generate a \code{\link{data.frame}} 
         with the combinations of the levels of the \code{recipient} 
         \code{\link{factor}s} in standard order. Although, the \code{\link{data.frame}s} 
         are not combined at this stage, the systematic design is 
         the combination, by columns, of the values of the \code{allocated} 
         \code{\link{factor}s} with the values of \code{recipient} 
         \code{\link{factor}s} in the \code{recipient} \code{\link{data.frame}}.

         The method of randomization described by Bailey (1981) is used to 
         randomize the \code{allocated} \code{\link{factor}s} to the 
         \code{recipient} \code{\link{factor}s}. That is, a permutation of the 
         \code{recipient} \code{\link{factor}s} is obtained that respects the 
         nesting for the design, but does not permute any of the factors in 
         the \code{except} vector. This permutation is 
         applied to the \code{recipient} \code{\link{factor}s}. Then the 
         \code{\link{data.frame}} containing the permuted \code{recipient} 
         \code{\link{factor}s} and that containng the unpermuted \code{allocated} 
         \code{\link{factor}s} are combined columnwise, as in \code{cbind}. To produce the 
         randomized layout, the rows of the combined \code{\link{data.frame}} are 
         reordered so that its \code{recipient} \code{\link{factor}s} are in either 
         standard order or, if a \code{\link{data.frame}} was suppled to 
         \code{recipient}, the same order as for the supplied \code{\link{data.frame}}. 

         The \code{.Units} and \code{.Permutation} \code{vectors} enable one to 
         swap between this permutation and the randomized layout. 
         The \emph{i}th value in \code{.Permutation} gives the unit to which 
         unit \emph{i} was assigned in the randomization. 
}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}}, \code{\link{designLatinSqrSys}}, \code{\link{designPlot}}, \code{\link{designAnatomy}} in package \pkg{dae}.}
\examples{
## Generate a randomized layout for a 4 x 4 Latin square
## (the nested.recipients argument is not needed here as none of the 
## factors are nested)
## Firstly, generate a systematic layout
LS.sys <- cbind(fac.gen(list(row = c("I","II","III","IV"), 
                             col = c(0,2,4,6))),
                treat = factor(designLatinSqrSys(4), label = LETTERS[1:4]))
## obtain randomized layout
LS.lay <- designRandomize(allocated = LS.sys["treat"], 
                          recipient = LS.sys[c("row","col")], 
                          seed = 7197132, unit.permutation = TRUE) 
LS.lay[LS.lay$.Permutation,]

## Generate a randomized layout for a replicated randomized complete 
## block design, with the block factors arranged in standard order for 
## rep then plot and then block
## Firstly, generate a systematic order such that levels of the 
## treatment factor coincide with plot
RCBD.sys <- cbind(fac.gen(list(rep = 2, plot=1:3, block = c("I","II"))),
                  tr = factor(rep(1:3, each=2, times=2)))
## obtain randomized layout
RCBD.lay <- designRandomize(allocated = RCBD.sys["tr"], 
                            recipient = RCBD.sys[c("rep", "block", "plot")], 
                            nested.recipients = list(plot = c("block","rep"), 
                                                     block="rep"), 
                            seed = 9719532, 
                            unit.permutation = TRUE)
#sort into the original standard order
RCBD.perm <- RCBD.lay[RCBD.lay$.Permutation,]
#resort into randomized order
RCBD.lay <- RCBD.perm[order(RCBD.perm$.Units),]

## Generate a layout for a split-unit experiment in which: 
## - the main-unit factor is A with 4 levels arranged in 
##   a randomized complete block design with 2 blocks;
## - the split-unit factor is B with 3 levels.
## Firstly, generate a systematic layout
SPL.sys <- cbind(fac.gen(list(block = 2, main.unit = 4, split.unit = 3)),
                 fac.gen(list(A = 4, B = 3), times = 2))
## obtain randomized layout
SPL.lay <- designRandomize(allocated = SPL.sys[c("A","B")], 
                           recipient = SPL.sys[c("block", "main.unit", "split.unit")], 
                           nested.recipients = list(main.unit = "block", 
                                                    split.unit = c("block", "main.unit")), 
                           seed=155251978)

## Generate a permutation of Seedlings within Species
seed.permute <- designRandomize(recipient = list(Species = 3, Seedlings = 4),
                                nested.recipients = list(Seedlings = "Species"),
                                seed = 75724, except = "Species", 
                                unit.permutation = TRUE)
}
\keyword{design}
\keyword{factor}
\keyword{datagen}
\references{Bailey, R.A. (1981) A unified approach to design of experiments. 
            \emph{Journal of the Royal Statistical Society, Series A}, 
            \bold{144}, 214--223.}