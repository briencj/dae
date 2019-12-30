\name{p2canon.object}
\alias{p2canon.object}
\title{Description of a p2canon object}
\description{An object of class \code{p2canon} that contains information derived from two 
\code{formulae} using \code{\link{projs.2canon}}.}
\value{A \code{list} of class \code{p2canon}. It has two components: \code{decomp} and 
       \code{aliasing}. The \code{decomp} component iscomposed as follows:
       \itemize{
          \item It has a component for each component of \code{Q1}. 
          \item Each of the components for \code{Q1} is a \code{list}; 
                 each of these \code{lists} has one component for each of 
                \code{Q2} and a component \code{Pres}. 
          \item Each of the \code{Q2} components is a \code{list} of three 
                components: \code{pairwise}, \code{adjusted} and \code{Qproj}. 
                 These components are  based on an eigenalysis of the 
                 relationship between the projectors for the parent \code{Q1} 
                 and \code{Q2} components. 
              \enumerate{
                 \item Each \code{pairwise} component is based on the nonzero 
                       canonical efficiency factors for the joint decomposition of 
                       the two parent projectors (see \code{\link{proj2.eigen}}). 
                 \item An \code{adjusted} component is based on the nonzero 
                       canonical efficiency factors for the joint decomposition of 
                       the \code{Q1} component and the \code{Q2} component, 
                       the latter adjusted for all \code{Q2} projectors that have 
                       occured previously in the \code{list}. 
                  \item The \code{Qproj} component is the adjusted projector for the 
                       parent \code{Q2} component.
              } 
          \item The \code{pairwise} and \code{adjusted} components have the 
                following components: \code{efficiencies}, \code{aefficiency}, 
                \code{mefficiency}, \code{sefficiency}, \code{eefficiency}, 
                \code{xefficiency}, \code{order} and \code{dforthog} 
                -- for details see \code{\link{efficiency.criteria}}.
        }
        The \code{aliasing} component is a data.frame decribing the aliasing between 
        terms corresponding to two \code{Q2} projectors when estimated in subspaces 
        corresponding to a \code{Q1} projector.
}
\author{Chris Brien}
\seealso{\code{\link{projs.2canon}}, \code{\link{designAnatomy}}, \code{\link{pcanon.object}}.}
\keyword{array}
\keyword{design}
\keyword{projector}
