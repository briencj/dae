\name{fac.genfactors}
\alias{fac.genfactors}
\title{Generate all combinations of the levels of the supplied factors, without replication}
\description{Generate all combinations of the levels of the supplied \code{factors}, without replication. 
       This function extracts the \code{levels} from the supplied \code{factors} and uses them to 
       generate the new factors. On the other hand, the levels must supplied in the \code{generate} 
       argument of the function \code{\link{fac.gen}}.}
\usage{fac.genfactors(factors, ...)}
\arguments{
 \item{factors}{A \code{\link{list}} of \code{\link{factor}s}, or an object of \code{\link{factor}s} 
                that is coercible to a \code{\link{list}}.}
 \item{...}{Further arguments passed to the \code{\link{fac.gen}} in creating 
            the \code{\link{data.frame}} of new \code{\link{factor}s}.}
}
\value{A \code{\link{data.frame}} whose columns correspond to \code{\link{factor}s} in the 
       \code{factors} \code{\link{list}}. The values in a column are the generated \code{levels} 
       of the \code{\link{factor}}. The number of rows in the \code{\link{data.frame}} will equal 
       the product of the numbers of levels of the supplied \code{\link{factor}s}.}
\details{The \code{levels} of each \code{\link{factor}} are generated in standard order,
      unless \code{order} is supplied to \code{\link{fac.gen}} via the `...' argument. 
      The \code{levels} of the new \code{\link{factor}s} will be in the same order as 
      in the supplied \code{factors}.}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}} in package \pkg{dae}}
\examples{
## generate a treatments key for the Variety and Nitrogen treatments factors in Oats.dat
data(Oats.dat)
trts.key <- fac.genfactors(factors = Oats.dat[c("Variety", "Nitrogen")])
trts.key$Treatment <- factor(1:nrow(trts.key))
}
\keyword{design}
\keyword{factor}
\keyword{datagen}
