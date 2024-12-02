\name{fac.multinested}
\alias{fac.multinested}
\title{Creates several factors, one for each level of the nesting factor and each of whose values are either generated within those of a level of the nesting factor or using the values of the nested factor within  the levels of the nesting factor.}
\description{Creates several \code{\link{factor}}s, one for each level of \code{nesting.fac} and each of whose values are either (i) generated within those of the level of \code{nesting.fac} or (ii) using the values of \code{nested.fac} within the levels of the \code{nesting.fac}. For (i), all elements having the same level of \code{nesting.fac} are numbered from 1 to the number of different elements having that level. For (ii), the values of \code{nested.fac} for a level of \code{nesting.fac} are copied. In both cases, for the values of \code{nested.fac} not equal to the level of the values of \code{nested.fac} for which a nested \code{\link{factor}} is being created, the levels are set to \code{outlevel} and labelled using \code{outlabel}. A \code{\link{factor}} is not created for a level of \code{nesting.fac} with label equal to \code{outlabel}. The names of the \code{\link{factor}}s are equal to the levels of \code{nesting.fac}; optionally \code{fac.prefix} is added to the beginning of the names of the \code{\link{factor}}s. The function is used to split up a nested term into separate terms for each level of \code{nesting.fac}.}
\usage{
fac.multinested(nesting.fac, nested.fac = NULL, fac.prefix = NULL, 
                nested.levs = NA, nested.labs = NA, 
                outlevel = 0, outlabel = "rest", ...)}
\arguments{
 \item{nesting.fac}{The \code{\link{factor}} for each of whose \code{levels} a nested facor                      is to be generated, except one is not generated for the \code{outlabel} level.}
 \item{nested.fac}{The \code{\link{factor}} whose values for a level are to be used for the 
                   \code{\link{factor}} being created for that level. If \code{nested.fac} 
                   is \code{NULL}, then the values of the \code{levels} is the the list of 
                   numbers from 1 to the replication of the level of \code{nesting.fac}, 
                   represented as characters.}
 \item{fac.prefix}{The prefix to be added to a level in naming a nested \code{\link{factor}}.}
 \item{nested.levs}{Optional \code{vector} of levels for the nested \code{\link{factor}}s. Any data 
       value that does not match a value in \code{nested.levs} will be \code{NA} in 
       the nested \code{\link{factor}} being generated. The default value of \code{nested.levs} is 
       the list of numbers from 1 to the replication of the \code{levels} of 
       \code{nesting.fac}, represented as characters. When \code{nested.fac} is not 
       \code{NULL}, \code{nested.levs} is not used.}
 \item{nested.labs}{Optional \code{vector} of values to use as \code{labels} for the 
       \code{levels} of the new nested \code{\link{factor}}s; they are combined with \code{outlabel}. 
       The default is \code{as.character(levels)}.} 
 \item{outlevel}{The level to use in the new \code{\link{factor}} for all values of old \code{\link{factor}} that do not correspond to the level of the \code{nesting.fac} to which the new \code{\link{factor}} corresponds.}
 \item{outlabel}{The label to use the \code{outlevel} level.} 
 \item{...}{Further arguments passed to the \code{\link{factor}} call creating a new \code{\link{factor}}.}
}
\value{A \code{\link{data.frame}} containing a \code{\link[=factor-class]{factor}}  
       for each level of \code{nesting.fac}.}
\note{The levels of \code{nesting.fac} do not have to be equally replicated.}
\author{Chris Brien}
\seealso{\code{\link{fac.gen}}, \code{\link{fac.nested}} in package \pkg{dae}, \code{\link{factor}}.}
\examples{
  lay <- data.frame(A = factor(rep(c(1:3), c(3,6,4)), labels = letters[1:3]))
  lay$B <-fac.nested(lay$A)

  #Add factors for B within each level of A
  lay2 <- cbind(lay, fac.multinested(lay$A))
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summary(canon2)

  #Add factors for B within each level of A, but with levels and outlabel given
  lay2 <- cbind(lay, fac.multinested(lay$A, nested.levs = seq(10,60,10), outlabel = "other"))
  \donttest{
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summary(canon2)
  }

  #Replicate the combinations of A and B three times and index them with the factor sample
  lay3 <- rbind(lay,lay,lay)
  lay3$sample <- with(lay3, fac.nested(fac.combine(list(A,B))))
  
  #Add factors for B within each level of A
  lay4 <- cbind(lay3, fac.multinested(nesting.fac = lay$A, nested.fac = lay$B))
  \donttest{
  canon4 <- designAnatomy(list(~(A/(a+b+c))/sample), data = lay4)
  summary(canon4)
  }

  #Add factors for sample within each combination of A and B
  lay5 <- with(lay4, cbind(lay4, 
                           fac.multinested(nesting.fac = a, fac.prefix = "a"),
                           fac.multinested(nesting.fac = b, fac.prefix = "b"),
                           fac.multinested(nesting.fac = c, fac.prefix = "c")))
  \donttest{
  canon5 <- designAnatomy(list(~A/(a/(a1+a2+a3)+b/(b1+b2+b3+b4+b5+b6)+c/(c1+c2+c3))), data = lay5)
  summary(canon5)}

  #Add factors for sample within each level of A
  lay6 <- cbind(lay4, 
                fac.multinested(nesting.fac = lay4$A, nested.fac = lay$sample, fac.prefix = "samp"))
  canon6 <- designAnatomy(list(~A/(a/sampa+b/sampb+c/sampc)), data = lay6)
  summary(canon6)
}
\keyword{factor}
\keyword{manip}