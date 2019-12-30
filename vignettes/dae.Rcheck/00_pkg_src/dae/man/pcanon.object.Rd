\name{pcanon.object}
\alias{pcanon.object}
\title{Description of a pcanon object}
\description{An object of class \code{pcanon} that contains information derived from several 
\code{formulae} using \code{\link{designAnatomy}}.}
\value{A \code{list} of class \code{pcanon} that has four components: (i) \code{Q}, 
       (ii) \code{terms}, (iii) \code{sources}, (iv) \code{marginality}, and (v) \code{aliasing}. 
       Each component is a \code{list} with as many components as there 
       are  formulae in the \code{formulae} \code{list} supplied to \code{\link{designAnatomy}}.
       
       The \code{Q} \code{list} is made up of the following components:
      \enumerate{
         \item The first component is the joint decomposition of two 
                structures derived from the first two formulae, being the 
                \code{\link{p2canon.object}} produced by \code{\link{projs.2canon}}. 
         \item Then there is a component for each further formulae; it contains 
               the \code{\link{p2canon.object}} obtained by applying 
               \code{\link{projs.2canon}} to the structure for a formula and 
               the already established joint decomposition of the structures 
               for the previous formulae in the \code{formulae}. 
         \item The last component contains the the \code{list} of the 
               projectors that give the combined canonical decomposition derived from 
               all of the \code{formulae}. 
       }
       The \code{terms}, \code{sources}, \code{marginalty} and \code{aliasing} 
       \code{\link{list}}s  have a component for each \code{\link{formula}} in the  
       \code{formulae} argument to \code{\link{designAnatomy}}, 
       Each component of the \code{terms} and \code{sources} \code{\link{list}}s  has 
       a \code{\link{character}} vector containing the terms or sources derived from its 
       \code{\link{formula}}. For the \code{marginality} component, each component is the 
       marginality \code{\link{matrix}} for the \code{terms} derived from its 
       \code{\link{formula}}. For the \code{aliasing} component, each component is the 
       aliasing \code{\link{data.frame}} for the \code{source} derived from its 
       \code{\link{formula}}. The components of these four \code{\link{list}}s  are 
       produced by \code{\link{pstructure.formula}} and are copied from the 
       \code{\link{pstructure.object}} for the \code{\link{formula}}.
       The names of the components of these four lists will be the names of the components 
       in the \code{formulae} list.

       The object has the attribute \code{labels}, which is set to \code{"terms"} or 
       \code{"sources"} according to which of these were used to label the projectors 
       when the object was created.
}

\author{Chris Brien}
\seealso{\code{\link{designAnatomy}}, \code{\link{p2canon.object}}.}
\keyword{array}
\keyword{design}
\keyword{projector}
