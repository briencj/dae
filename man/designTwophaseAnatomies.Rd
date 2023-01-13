\name{designTwophaseAnatomies}
\alias{designTwophaseAnatomies}
\title{Given the layout for a design and three structure formulae, obtain the anatomies for the (i) two-phase, (ii) first-phase, (iii) cross-phase, treatments, and (iv) combined-units designs.}
\description{Uses \code{\link{designAnatomy}} to obtain the four species of designs, described by Brien (2019), that are associated with a standard two-phase design: the anatomies for the (i) two-phase, (ii) first-phase, (iii) cross-phase, treatments, and (iv) combined-units designs. (The names of the last two designs in  Brien (2019) were cross-phase and second-phase designs.) For the standard two-phase design, the first-phase design is the design that allocates first-phase treatments to first-phase units. The cross-phase, treatments design allocates the first-phase treatments to the second-phase units and the combined-units design allocates the the first-phase units to the second-phase units. The two-phase design combines the other three species of designs. However, it is not mandatory that the three formula correspond to second-phase units, first-phase units and first-phase treatments, respectively, as is implied above; this is just the correspondence for a standard two-phase design. The essential requirement is that three structure formulae are supplied. For example, if there are both first- and second-phase treatments in a two-phase design, the third formula might involve the treatment factors from both phases. In this case, the default anatomy titles when printing occurs will not be correct, but can be modifed using the \code{titles} argument.
}
\usage{
designTwophaseAnatomies(formulae, data, which.designs = "all", 
                        printAnatomies = TRUE, titles,
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
               \code{combined-units}; \code{all}, the default, results in all 
               four being obtained.}
 \item{printAnatomies}{A \code{\link{logical}} indicating whether or not the anatomies 
               are to be printed.}
 \item{titles}{A \code{\link{character}} vector of length four providing titles for 
                      the printed anatomies. It should have the titles, in the following 
                      order, for the antomies based on : (i) all three \code{formulae}, 
                      (ii) the second and third \code{formulae}, (iii) the first and 
                      third \code{formulae}, and (iv) the first and second \code{formulae}. 
                      If any element of \code{titles} is \code{NA} then that element is 
                      replaced with the corresponding default element of \code{titles}, 
                      these being, in order: 
                      \code{Anatomy for the full two-phase design};  
                      \code{Anatomy for the first-phase design}; 
                      \code{Anatomy for the cross-phase, treatments design}; and  
                      \code{Anatomy for the combined-units design}. The titles generated will 
                      be saved as an attribute of the returned \code{list}.}
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
                      parallel with the components of \code{formulae}, irrespective of 
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
         \code{cross} and \code{combined}.Each contains the \code{\link{pcanon.object}} 
         for one of the four designs produced by \code{designTwophaseAnatomies}, unless it is 
         \code{NULL} because the design was omitted from the \code{which.designs} 
         argument. The returned \code{list} has an attribute \code{titles}, being a 
         \code{character} vector of length four and containing the titles used in 
         printing the anatomies.}
\details{To produce the anatomies, \code{\link{designAnatomy}} is called. The 
         two-phase anatomy is based on the three \code{formulae} supplied in \code{formulae}, 
         the first-phase anatomy uses the second and third \code{formulae}, the cross-phase, 
         treatments anatomy derives from the first and third \code{formulae} and the combined-units 
         anatomy is obtained with the first and second \code{formulae}.}
\references{Brien, C. J. (2017) Multiphase experiments in practice: A look back. 
            \emph{Australian & New Zealand Journal of Statistics}, \bold{59}, 327-352.

            Brien, C. J. (2019) Multiphase experiments with at least one later 
            laboratory phase . II. Northogonal designs. 
            \emph{Australian & New Zealand Journal of Statistics}\bold{61}, 
            234-268.
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
