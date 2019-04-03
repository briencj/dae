\name{summary.pcanon}
\alias{summary.pcanon}
\alias{summary,pcanon-method}
\title{Summarizes the anatomy of a design, being the decomposition of the sample space based on its 
             canonical analysis, as produced by designAnatomy}
\description{Gives the anatomy of a design in a table; it summarizes the joint decomposition of two or 
             more sets of projectors (Brien and Bailey, 2009) obtained using 
             \code{\link{designAnatomy}}. It includes the efficiency criteria computed 
             from the canonical efficiency factors for the joint decomposition. The labels in 
             the table may be Terms or Sources. The terms are those that would be included in a 
             mixed model for an experiment based on the design. The sources are the orthogonal 
             subspaces, derived from the terms, that make up the decomposition and the degrees 
             of freedom and efficiency factors relate to these subspaces. The table displays 
             how the information for the different sources allowed for in the design are related. 
             For more information about the notation used for sources see the \code{labels} argument of 
             \code{\link{designAnatomy}}.
             
             It is possible to supply an \code{object} that is a \code{\link{pcanon.object}} produced in 
             versions prior to 3.0-0 using \code{projs.canon}.
}
\usage{
\method{summary}{pcanon}(object, labels.swap = FALSE, 
        which.criteria = c("aefficiency", "eefficiency", "order"), ...)}
\arguments{
 \item{object}{A \code{\link{pcanon.object}}.}
  \item{labels.swap}{A \code{\link{logical}} indicating whether to swap between "sources" and 
                `terms' in the output. The default is established by the \code{labels} 
                argument of \code{\link{designAnatomy}} and \code{\link{projs.canon}}.}
 \item{which.criteria}{A \code{\link{character}} \code{vector} nominating the efficiency 
                       criteria to be included in the summary. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency},
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see 
                       \code{\link{efficiency.criteria}}. 
                       If there is only one formula, this 
                       argument is ignored.}
 \item{...}{further arguments affecting the summary produced.}
}
\value{An object of class \code{summary.pcanon} that is a \code{\link{list}} with the two
       components \code{decomp} and \code{aliasing}. 

       The component \code{decomp} is a \code{data.frame} whose rows correspond to subspaces 
       in the decomposition for a design. It has the following \code{attribute}s: 
       (i) \code{title} that is the title for printing with the decomposition table; 
       (ii) \code{ntiers} that is equal to the number of tiers; (iii) \code{orthogonal} that is 
       \code{TRUE} if the design is orthogonal; (iv) \code{labels} that is either "terms" or 
       "sources" depending on the \code{labels} that have resulted from the setting 
       of \code{label.swap}.

      The component \code{aliasing} is a \code{data.frame} that is also of \code{class} 
      \code{aliasing}. It contains information about the aliasing between terms that are 
      derived from the same formula and has the attribute \code{title} that is the title 
      to be printed with the aliasing table.

      However, if the \code{object} supplied is a \code{\link{pcanon.object}} produced with  
      versions prior to 3.0-0 using \code{projs.canon}, the value is a \code{data.frame}, 
      instead of a \code{list}, that has the same \code{attribute}s as the \code{decomp} 
      component of the \code{summary.pcanon} object now produced, except that \code{labels} 
      is always set to "terms".
}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{designAnatomy}}, \code{\link{designAnatomy}}, ,
         \code{\link{pstructure}}, \code{\link{efficiency.criteria}}, 
         \code{\link{proj2.combine}}, \cr 
         \code{\link{proj2.efficiency}}, \code{\link{proj2.eigen}}, 
          \code{\link{print.summary.pcanon}}in package \pkg{dae}, 
         \code{\link{eigen}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain combined decomposition and summarize
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), 
                                data = PBIBD2.lay)
summary(unit.trt.canon, which = c("aeff","eeff","order"))
summary(unit.trt.canon, which = c("aeff","eeff","order"), labels.swap = TRUE)
}
\keyword{array}
\keyword{design}
\keyword{projector}
