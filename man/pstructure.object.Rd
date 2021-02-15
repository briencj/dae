\name{pstructure.object}
\alias{pstructure.object}
\title{Description of a pstructure object}
\description{An object of class \code{pstructure} that contains information derived from a 
\code{\link{formula}} using \code{\link{pstructure.formula}}. It also inherits from class \code{list}.}
\value{A \code{list} of class \code{pstructure}  with the following components: 
       \enumerate{
          \item Q: a list with a component of class \code{projector}, being  
                  the orthogonalized projectors for each non-aliased term/source 
                  in the \code{formula}; if \code{grandMean} is \code{TRUE} in 
                  the call to \code{\link{pstructure.formula}} then it also 
                  includes the \code{projector} for it;
          \item terms: a \code{\link{character}} vector with the non-aliased 
                  term names; if \code{grandMean} is \code{TRUE} in 
                  the call to \code{\link{pstructure.formula}} then the first 
                  term will be "\code{Mean}"; 
          \item sources: a \code{\link{character}} vector with the non-aliased
                  source names;
          \item marginality: a \code{\link{matrix}} of zeroes and ones with the same 
                number of rows and columns as number of non-aliased terms, excluding the term for the grand mean 
                even when \code{grandMean} is \code{TRUE}; the row names and column names 
                are the elements \code{terms}, excluding "\code{Mean}"; 

                the entry in the ith row and jth column will be one if the 
                ith term is marginal to the jth term i.e. the column space of the 
                ith term is a subspace of that for the jth term and so the source for 
                the jth term will have been made orthogonal to that for the ith term; 
                otherwise, the entries are zero. 
          \item aliasing: a \code{\link{data.frame}} containing the information about the 
                (partial) aliasing between the sources in the \code{formula}. 
                The columns are:
                \itemize{
                  \item Source: the source names, or associated term name, for those that are 
                        (partially) aliased with previous sources;
                  \item df: the remaining degrees of freedom for the source;
                  \item Alias: the source with which the current entry is (partially) aliased;
                  \item efficiency criteria: a set of columns for the complete set of criteria 
                        calculated by \code{\link{efficiency.criteria}}; the criteria reflect 
                        the amount of information that is aliased with previous sources and a 
                        line is included in the component that reports the informaton remaining 
                        after adjustment for previous sources.
                }
                The information provided depends on the setting of \code{orthogonalize}. 
                All the information is provided for the \code{"hybrid"} option. For the 
                option \code{"differencing"}, no efficiency criteria are included and either 
                the terms/sources of the \code{Alias} are set to \code{"unknown"} and the 
                \code{df} are set to \code{NA} 
                when these are unknown. For the option \code{"eigenmethods"}, the previous 
                terms/sources cannot be identified and so all values of \code{Alias} are set 
                to \code{NA}. If there is no (partial) aliasing then the component is set to 
                \code{NULL}.
       }
       The object has the attribute \code{labels}, which is set to \code{"terms"} or 
       \code{"sources"} according to which of these label the projectors.
}
\author{Chris Brien}
\seealso{\code{\link{pstructure.formula}} and, for further information about the projector classs, 
         \code{\link[=projector-class]{projector}}.}
\keyword{array}
\keyword{design}
\keyword{projector}
