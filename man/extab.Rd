\name{extab}
\alias{extab}
\title{Expands the values in table to a vector}
\description{Expands the values in \code{table} to a \code{vector}
             according to the \code{index.factors} that are considered to index 
             the \code{table}, either in standard or Yates order.  The order 
             of the values in the \code{vector} is determined by the order of 
             the values of the \code{index.factors}.}
\usage{extab(table, index.factors, order="standard")}
\arguments{
 \item{table}{A numeric \code{vector} containing the values to be expanded. Its length
  must equal the product of the number of used levels for the \code{\link{factor}s} in 
  \code{index.factors} and the values in it correspond to all levels combinations 
  of these \code{\link{factor}s}. That is, the values of the \code{index.factors} are 
  irrelevant to \code{table}.}
 \item{index.factors}{A list of \code{\link{factor}s} that index the \code{table}. 
   All the \code{\link{factor}s} must be the same length.}
 \item{order}{The order in which the levels combinations of the \code{index.factors} 
   are to be considered as numbered in indexing \code{table}; \code{standard} numbers them as 
   if they are arranged in standard order, that is with the first factor moving 
   slowest and the last factor moving fastest; \code{yates} numbers them as if they are 
   arranged in Yates order, that is with the first factor moving fastest and 
   last factor moving slowest.}
}
\value{A \code{vector} of length equal to the \code{\link{factor}s} in 
       \code{index.factor} whose values are  taken from \code{table}.}
\author{Chris Brien}
\examples{
## generate a small completely randomized design with the two-level 
## factors A and B 
n <- 12
CRD.unit <- list(Unit = n)
CRD.treat <- fac.gen(list(A = 2, B = 2), each = 3)
CRD.lay <- designRandomize(allocated = CRD.treat, recipient = CRD.unit, 
                           seed = 956)

## set up a 2 x 2 table of A x B effects	
AB.tab <- c(12, -12, -12, 12)

## add a unit-length vector of expanded effects to CRD.lay
attach(CRD.lay)
CRD.lay$AB.effects <- extab(table=AB.tab, index.factors=list(A, B))
}
\keyword{manip}