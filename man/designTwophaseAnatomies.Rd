\name{designTwophaseAnatomies}
\alias{designTwophaseAnatomies}
\title{Given the layout for a design and three structure formulae, obtain the anatomies for the two-phase, first-phase, cross-phase and second-phase designs.}
\description{Uses \code{\link{designAnatomy}} to obtain the four species of designs associated with a two-phase design that are described by Brien (2019): the anatomies for the two-phase, first-phase, cross-phase and second-phase designs.
}
\usage{
designTwophaseAnatomies(formulae, data, 
                       which.designs = "all", printAnatomies = TRUE, 
                       orthogonalize = "hybrid", 
                       marginality = NULL, 
                       which.criteria = c("aefficiency", "eefficiency", 
                                          "order"), ...)
}
\arguments{
 \item{formulae}{An object of class \code{\link{list}} with three components of class 
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
 \item{which.designs}{A \code{\link{character}} vector indicating the species of 
               designs that are to be obtained. It should include one or more of
               \code{two-phase}, \code{first-phase}, \code{cross-phase} and 
               \code{second-phase}; \code{all}, the default, results in all 
               four being obtained.}
 \item{printAnatomies}{A \code{\link{logical}} indicating whether or not the anatomies 
               are to be printed.}
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
 \item{which.criteria}{A \code{\link{character}} vector nominating the efficiency 
                       criteria to be included in the summary of aliasing 
                       between terms within a structure. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency},
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}. 
                       If \code{none}, no summary is printed.}
 \item{...}{further arguments passed to \code{designAnatomy}.}
}
\value{A \code{\link{list}} containing the components \code{twophase}, \code{first}, 
         \code{cross} and \code{second}.Each contains the \code{\link{pcanon.object}} 
         for one of the four designs produced by designTwophaseAnatomies, unless it is 
         \code{NULL} because the design was omitted from the \code{which.designs} 
         argument.}
\details{To produce the anatomies, \code{\link{designAnatomy}} is called. The 
         two-phase anatomy is based on the three formulae supplied in \code{formulae}, 
         the first-phase anatomy uses the second and third formulae, the cross-phase 
         anatomy derives from the first and third formulae and the second-phase 
         anatomy is obtained with the first and second formulae.}
\references{Brien, C. J. (2017) Multiphase experiments in practice: A look back. 
            \emph{Australian & New Zealand Journal of Statistics}, \bold{59}, 327-352.

            Brien, C. J. (2019) Multiphase experiments with at least one later 
            laboratory phase . II. Northogonal designs. 
            \emph{Australian & New Zealand Journal of Statistics}, accepted for publication.
}
\author{Chris Brien}
\seealso{\code{\link{designAnatomy}}, 
         \code{\link{pcanon.object}}, \code{\link{p2canon.object}}, 
         \code{\link{summary.pcanon}}, \code{\link{efficiencies.pcanon}}, \cr
         \code{\link{pstructure}} , 
         \code{\link{projs.2canon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, 
         \code{\link{proj2.eigen}},  \cr 
         \code{\link{efficiency.criteria}}, in package \pkg{dae},  
         \code{\link{eigen}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
  #'## Microarray example from Jarrett & Ruggiero (2008) - see Brien (2019)
  jr.lay <- fac.gen(list(Set = 7, Dye = 2, Array = 3))
  jr.lay <- within(jr.lay, 
                   { 
                     Block <- factor(rep(1:7, each=6))
                     Plant <- factor(rep(c(1,2,3,2,3,1), times=7))
                     Sample <- factor(c(rep(c(2,1,2,2,1,1, 1,2,1,1,2,2), times=3), 
                                        2,1,2,2,1,1))
                     Treat <- factor(c(1,2,4,2,4,1, 2,3,5,3,5,2, 3,4,6,4,6,3, 
                                       4,5,7,5,7,4, 5,6,1,6,1,5, 6,7,2,7,2,6, 
                                       7,1,3,1,3,7),
                                     labels=c("A","B","C","D","E","F","G"))
                   })
  
  jr.anat <- designTwophaseAnatomies(formulae = list(array = ~ (Set:Array)*Dye,
                                                     plot = ~ Block/Plant/Sample,
                                                     trt = ~ Treat),
                                     which.designs = c("first","cross"), 
                                     data = jr.lay)  

## Three-phase sensory example from Brien and Payne (1999)
\dontrun{
data(Sensory3Phase.dat)
Sensory.canon <- designTwophaseAnatomies(formulae = list(
                              eval= ~ ((Occasions/Intervals/Sittings)*Judges)/Positions, 
                              field= ~ (Rows*(Squares/Columns))/Halfplots,
                              treats= ~ Trellis*Method),
                                        data = Sensory3Phase.dat)
}
}
\keyword{array}
\keyword{design}
\keyword{projector}
