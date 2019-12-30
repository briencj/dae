\name{designAnatomy}
\alias{designAnatomy}
\title{Given the layout for a design, obtain its anatomy via the canonical analysis 
       of its projectors to show the confounding and 
       aliasing inherent in the design.}
\description{Computes the canonical efficiency factors for the joint 
       decomposition of two or more structures or sets of mutually orthogonally 
       projectors (Brien and Bailey, 2009; Brien, 2017; Brien, 2019), orthogonalizing 
       projectors in a set to those earlier in the set of projectors with 
       which they are partially aliased. The results can be summarized in the 
       form of a decomposition table that shows the confounding between sources 
       from different sets. For examples of its use also see the vignette 
       \cite{daeDesignNotes.pdf}.
}
\usage{
designAnatomy(formulae, data, keep.order = TRUE, grandMean = FALSE, 
              orthogonalize = "hybrid", labels = "sources", 
              marginality = NULL, check.marginality = TRUE, 
              which.criteria = c("aefficiency","eefficiency","order"), 
              aliasing.print = FALSE, 
              omit.projectors = c("pcanon", "combined"), ...)}
\arguments{
 \item{formulae}{An object of class \code{\link{list}} whose components are of class 
                      \code{\link{formula}}. 
                      Usually, the terms in a single formula have the same status in the 
                      allocation of factors in the design. For example, all involve only 
                      factors that were allocated, or all involve factors that were 
                      recipients of allocated factors. The names of the components are 
                      used to identify the 
                      sources in the \code{\link{summary.pcanon}} object. They will also 
                      be used to name the \code{terms}, \code{sources} and \code{marginality} 
                      lists in the \code{\link{pcanon.object}}.}
 \item{data}{A \code{\link{data.frame}} contains the values of the factors and variables 
             that occur in \code{formulae}.}
 \item{keep.order}{A \code{\link{logical}} indicating whether the terms should keep 
                   their position in the expanded \code{formula} projector, or 
                   reordered so that main effects precede two-factor 
                   interactions, which precede three-factor interactions and 
                   so on.}
 \item{grandMean}{A \code{\link{logical}} indicating whether the projector for the 
                       grand mean is to be included for each structure.}
 \item{orthogonalize}{A \code{\link{character}} vector indicating the method for orthogonalizing a
                      projector to those for terms that occurred previously 
                      in a single formula. Three options are available:  
                      \code{hybrid}; \code{differencing}; \code{eigenmethods}.
                      The \code{hybrid} option is the most general and uses the 
                      relationships between the projection operators for the 
                      terms in the \code{formula} to decide which projectors to 
                      substract and which to orthogonalize using eigenmethods. 
                      The \code{differencing} option subtracts, from the current 
                      projector, those previously orthogonalized projectors for 
                      terms whose factors are a subset of the  current projector's 
                      factors. The \code{eigemethods} option recursively 
                      orthogonalizes the projects using an eigenanalysis of 
                      each projector with previously orthogonalized projectors.
                      If a single value is given, it is used for all formulae.}
  \item{labels}{A \code{\link{character}} nominating the type of labels to be used in 
                      labelling the projectors, and which will be used also 
                      in the output tables, such the tables of the aliasing 
                      in the structure. The two alternatives are \code{terms} and 
                      \code{sources}. Terms have all factors/variables in it 
                      separated by colons (\code{:}). Sources have factors/variables 
                      in them that represent interactions separated by hashes 
                      (\code{#}); if some factors are nested within others, the 
                      nesting factors are surrounded by square brackets 
                      (\code{[} and \code{]}) and separated by colons (\code{:}). 
                      If some generalized, or combined, factors have no marginal 
                      terms, the constituent factors are separated by colons 
                      (\code{:}) and if they interact with other factors in 
                      the source they will be parenthesized.}
 \item{marginality}{A \code{\link{list}} that can be used to supply some or all of the 
                      marginality matrices when it is desired to overwrite calculated 
                      marginality matrices or when they are not calculated. 
                      If the \code{\link{list}} is the same length as the 
                      \code{formulae} \code{\link{list}}, they will be associated in 
                      parallel with the components of \code{formulae}, irespective of 
                      the naming of the two \code{\link{list}}s. If the number of 
                      components in \code{marginlaity} is less than the number of 
                      components in \code{formulae} then both \code{\link{list}}s 
                      must be named so that those in the \code{marginality} 
                      \code{\link{list}} can be matched with those in the 
                      \code{formulae} \code{\link{list}}.

                      Each component of the \code{marginality} \code{\link{list}} 
                      must be either \code{NULL} or a square 
                      \code{\link{matrix}} consisting of zeroes and ones that 
                      gives the marginalites of the terms in the formula. It must have 
                      the row and column names set to the terms from the expanded 
                      \code{formula}, including being in the same order as these terms. 
                      The entry in the ith row and jth column will be one if the 
                      ith term is marginal to the jth term i.e. the column space of the 
                      ith term is a subspace of that for the jth term and so the source 
                      for the jth term is to be made orthogonal to that for the ith term. 
                      Otherwise, the entries are zero. A row and column should not be 
                      included for the grand mean even if \code{grandMean} is \code{TRUE}.}
 \item{check.marginality}{A \code{\link{logical}} indicating whether the marginality matrix, 
                       when it is supplied, is to be checked against that computed 
                       by \code{\link{pstructure.formula}}. It is ignored when 
                       \code{orthogonalize} is set to \code{eigenmethods}.}
 \item{which.criteria}{A \code{\link{character}} vector nominating the efficiency 
                       criteria to be included in the summary of aliasing 
                       between terms within a structure. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency},
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}. 
                       If \code{none}, no summary is printed.}
 \item{aliasing.print}{A \code{\link{logical}} indicating whether the aliasing between sources 
                       is to be printed.}
 \item{omit.projectors}{A \code{\link{character}} vector of the types of  projectors to 
                      omit from the returned \code{pcanon} object. If 
                      \code{pcanon} is included in the vector then the 
                      projectors in these objects will be replaced with a 
                      \code{numeric} containing their degrees of freedom. 
                      If \code{combined} is included in the vector 
                      then the projectors for the combined decomposition will be 
                      replaced with a \code{numeric} containing their degrees of 
                      freedom. If \code{none} is included in the vector then no 
                      projectors will be omitted.}
 \item{...}{further arguments passed to \code{terms}.}
}
\value{A \code{\link{pcanon.object}}.}
\details{For each formula supplied in \code{formulae}, the set of projectors is 
         obtained using \code{\link{pstructure}}; there is one projector 
         for each term in a formula. Then \code{\link{projs.2canon}} is used 
         to perform an analysis of the canonical relationships between two sets 
         of projectors for the first two formulae. If there are further formulae, 
         the relationships between its projectors and the already established 
         decomposition is obtained using \code{\link{projs.2canon}}. The core 
         of the analysis is the determination of eigenvalues of the product of 
         pairs of projectors using the results of James and Wilkinson (1971). 
         However, if the order of balance between two projection matrices is 
         10 or more or the James and Wilkinson (1971) methods fails to produce 
         an idempotent matrix, equation 5.3 of Payne and Tobias (1992) is used 
         to obtain the projection matrices for their joint decompostion.}
\references{Brien, C. J. (2017) Multiphase experiments in practice: A look back. 
            \emph{Australian & New Zealand Journal of Statistics}, \bold{59}, 327-352.

            Brien, C. J. (2019) Multiphase experiments with at least one later 
            laboratory phase . II. Northogonal designs. 
            \emph{Australian & New Zealand Journal of Statistics}, accepted for publication.

            Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.

           James, A. T. and Wilkinson, G. N. (1971) Factorization of the residual 
           operator and canonical decomposition of nonorthogonal factors in the 
           analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.
           
           Payne, R. W. and R. D. Tobias (1992). General balance, combination of 
            information and the analysis of covariance. 
            \emph{Scandinavian Journal of Statistics}, \bold{19}, 3-23.}
\author{Chris Brien}
\seealso{\code{\link{designRandomize}}, \code{\link{designLatinSqrSys}}, \code{\link{designPlot}}, \cr
         \code{\link{pcanon.object}}, \code{\link{p2canon.object}}, 
         \code{\link{summary.pcanon}}, \code{\link{efficiencies.pcanon}}, 
         \code{\link{pstructure}} , \cr
         \code{\link{projs.2canon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, 
         \code{\link{proj2.eigen}},  \code{\link{efficiency.criteria}}, 
         in package \pkg{dae}, \cr 
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
unit.trt.canon <- designAnatomy(formulae = list(unit=~ Block/Unit, trt=~ trt),
                                data = PBIBD2.lay)
summary(unit.trt.canon, which.criteria = c("aeff","eeff","order"))
summary(unit.trt.canon, which.criteria = c("aeff","eeff","order"), labels.swap = TRUE)

## Three-phase sensory example from Brien and Payne (1999)
\dontrun{
data(Sensory3Phase.dat)
Eval.Field.Treat.canon <- designAnatomy(formulae = list(
                              eval= ~ ((Occasions/Intervals/Sittings)*Judges)/Positions, 
                              field= ~ (Rows*(Squares/Columns))/Halfplots,
                              treats= ~ Trellis*Method),
                                        data = Sensory3Phase.dat)
summary(Eval.Field.Treat.canon, which.criteria =c("aefficiency", "order"))
}
}
\keyword{array}
\keyword{design}
\keyword{projector}
