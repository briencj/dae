\name{elements}
\alias{elements}
\title{Extract the elements of an array specified by the subscripts}
\description{Elements of the \code{array} \code{x} corresponding to the rows of the two dimensional 
             object \code{subscripts} are extracted. The number of columns of \code{subscripts} 
             corresponds to the number of dimensions of \code{x}. 
             The effect of supplying less columns in \code{subscripts} than the 
             number of dimensions in \code{x} is the same as for \code{"["}.}
\usage{elements(x, subscripts)}
\arguments{
 \item{x}{An \code{array} with at least two dimensions whose elements are to be extracted.}
 \item{subscripts}{A two dimensional object interpreted as elements by dimensions.}
}
\value{A \code{vector} containing the extracted elements and whose length equals the 
       number of rows in the \code{subscripts} object.}
\author{Chris Brien}
\seealso{\code{Extract}}
\examples{
## Form a table of the means for all combinations of Row and Line.
## Then obtain the values corresponding to the combinations in the data frame x,
## excluding Row 3.
x <- fac.gen(list(Row = 2, Line = 4), each =2)
x$y <- rnorm(16)
RowLine.tab <- tapply(x$y, list(x$Row, x$Line), mean)
xs <- elements(RowLine.tab, subscripts=x[x$"Line" != 3, c("Row", "Line")])
}
\keyword{array}
\keyword{manip}