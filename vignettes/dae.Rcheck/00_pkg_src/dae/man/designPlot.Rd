\name{designPlot}
\alias{designPlot}
\title{A graphical representation of an experimental design using labels stored in a matrix.}
\description{This function uses labels, usually derived from treatment and blocking factors from 
  an experimental design and stored in a matrix, to build a graphical representation of 
  the matrix, highlighting the position of certain labels . 
  It is a modified version of the function supplied with DiGGer. 
  It includes more control over the labelling of the rows and columns 
  of the design and allows for more flexible plotting of designs with 
  unequal block size.}
\usage{
designPlot(designMatrix, labels = NULL, altlabels = NULL, plotlabels = TRUE, 
           rtitle = NULL, ctitle = NULL, 
           rlabelsreverse = FALSE, clabelsreverse = FALSE, 
           font = 1, chardivisor = 2, rchardivisor = 1, cchardivisor = 1, 
           cellfillcolour = NA, plotcellboundary = TRUE, 
           rcellpropn = 1, ccellpropn = 1, 
           blocksequence = FALSE, blockdefinition = NULL, 
           blocklinecolour = 1, blocklinewidth = 2, 
           rotate = FALSE, new = TRUE, ...)
}
\arguments{
  \item{designMatrix}{A \code{\link{matrix}} containing a set of numerics or characters 
                   being the labels as they have been assigned to the cells of the grid 
                   represented by the \code{\link{matrix}}.}
  \item{labels}{A \code{\link{numeric}} or \code{\link{character}} vector giving the 
                cells in \code{designMatrix} that are to be plotted in this call to 
                \code{designPlot}. If \code{NULL} then all the cells are plotted.

                What is actually plotted for a cell is controlled jointly by \code{labels}, 
                \code{plotlabels}, \code{altlabels}, \code{plotcellboundary} and 
                \code{cellfillcolour}. If \code{plotlabels} is TRUE and \code{altlabels} 
                is \code{NULL} then \code{labels} are plotted in the cells, unless 
                \code{labels} is \code{NULL} when the labels in \code{designMatrix} 
                are plotted. 

                Whatever is being plotted, \code{altlabels} and 
                \code{cellfillcolour} must have an appropriate number of values.
                See \code{\link{text}} for more information on specifying the labels.}
  \item{altlabels}{Either a \code{\link{character}} vector containing an alternative set of 
                   labels for the \code{labels} currently being plotted or a single \code{\link{integer}} 
                   specifying an alternative symbol to be used in plotting cells when 
                   \code{plotlabels} is \code{TRUE}. The length of \code{altlabels} must be one or 
                   the same length as \code{labels}, unless \code{labels} is \code{NULL} in which case 
                   it must equal the number of unique labels in \code{designMatrix}.  
                   
                   If \code{altlabels} is \code{NULL}, the labels specified in \code{labels} are 
                   plotted when \code{plotlabels} is \code{TRUE}. If \code{labels} is also \code{NULL}, 
                   the labels in \code{designMatrix} are plotted. 
                   See \code{\link{text}} for more information on specifying the labels.}
  \item{plotlabels}{A \code{\link{logical}} to indicate whether labels are to be 
                    plotted in the cells. If TRUE, print all labels or 
                    the specific labels listed in \code{labels}. If FALSE, no labels are 
                    printed in the cells.}
  \item{rtitle}{A \code{\link{character}} string to use as a title for rows of the plot.
    If \code{rtitle} is \code{NULL} then no title is plotted.}
  \item{ctitle}{A \code{\link{character}} string to use as a title for columns of the plot. 
    If \code{ctitle} is \code{NULL} then no title is plotted.}
  \item{rlabelsreverse}{A \code{\link{logical}} indicating whether to reverse the row labels.}
  \item{clabelsreverse}{A \code{\link{logical}} indicating whether to reverse the column labels.}
  \item{font}{An \code{\link{integer}} specifying the font to be used for row and column labelling. 
    See \code{\link{par}} for further details.}
  \item{chardivisor}{A \code{\link{numeric}} that changes the size of text and symbols in the cells 
        by dividing the default size by it.}
  \item{rchardivisor}{A \code{\link{numeric}} that changes the size of the labels of the rows of the
        design by dividing the default size by it.}
  \item{cchardivisor}{A \code{\link{numeric}} that changes the size of the labels of the columns of the
        design by dividing the default size by it.}
  \item{cellfillcolour}{A \code{\link{character}} string specifying the colour of the fill for 
                        the cells to be plotted in this call. If there is only one colour 
                        then all cells being plotted with that colour. If there is more than one
                        colour then, unless \code{labels} is \code{NULL}, the number of colours must 
                        at least equal the number of labels 
                        and then the fill colours will be matched, one for one from the first colour,
                        with the labels. If \code{labels} is \code{NULL} then the number of colours 
                        must  at least equal the number of unique labels in \code{designMatrix}. 
                        The default, \code{NA}, is to leave ther cells unfilled.
    
    See also \code{Colour specification} under the \code{\link{par}} function.}
  \item{plotcellboundary}{A \code{\link{logical}} indicting whether a boundary is to plotted 
        around a cell.}
  \item{rcellpropn}{a value between 0 and 1 giving the proportion of the standard row size of 
        a cell size to be plotted as a cell.}
  \item{ccellpropn}{a value between 0 and 1 giving the proportion of the standard column size of 
        a cell size to be plotted as a cell.}
  \item{blocksequence}{A \code{\link{logical}} that determines whether block numbers are repetitions 
    or sequences of block numbers.}
  \item{blockdefinition}{A \code{\link{matrix}} of block sizes:
      \itemize{
        \item if there is only one row, then the first element is interpreted as the no. rows in 
        each block and blocks with this number of rows are to be repeated  across the rows of the design.
        \item if there is more than one row, then each row of the matrix specifies a block, 
        with the sequence of rows in the matrix specifying a corresponding 
        sequence of blocks down the rows of the design.}
    Similarly, a single value for a column specifies a repetition of blocks of that size 
    across the columns of the design, while several column values specifies a 
    sequence of blocks across the columns of the size specified.}
  \item{blocklinecolour}{A \code{\link{character}} string specifying the colour of the block boundary.
    
    See also \code{Colour specification} under the \code{\link{par}} function.}
  \item{blocklinewidth}{A \code{\link{numeric}} giving the width of the block boundary to be plotted.}
  \item{rotate}{A \code{\link{logical}} which, if \code{TRUE}, results in the matrix being rotated 
    90 degrees for plotting.}
  \item{new}{A \code{\link{logical}} indicating if a new plot is to be produced or the current 
    plot is added to.}
  \item{...}{further arguments passed to \code{\link{polygon}} in plotting the cell.}
}
\value{no values are returned, but a plot is produced.}

\references{Coombes, N. E. (2009). \emph{DiGGer design search tool in R}. 
  \url{http://nswdpibiom.org/austatgen/software/}}


\author{Chris Brien}
\seealso{\code{\link{blockboundaryPlot}}, \code{\link{designPlotlabels}}, \code{\link{designLatinSqrSys}}, \code{\link{designRandomize}}, \code{\link{designAnatomy}}  
          in package \pkg{dae}.  \cr
          Also, \code{\link{par}}, \code{\link{polygon}}, 
  \code{DiGGer}}
\examples{\dontrun{
  designPlot(des.mat, labels=1:4, cellfillcolour="lightblue", new=TRUE, 
             plotcellboundary = TRUE, chardivisor=3, 
             rtitle="Lanes", ctitle="Positions", 
             rcellpropn = 1, ccellpropn=1)
  designPlot(des.mat, labels=5:87, plotlabels=TRUE, cellfillcolour="grey", new=FALSE,
             plotcellboundary = TRUE, chardivisor=3)
  designPlot(des.mat, labels=88:434, plotlabels=TRUE, cellfillcolour="lightgreen", 
             new=FALSE, plotcellboundary = TRUE, chardivisor=3, 
             blocksequence=TRUE, blockdefinition=cbind(4,10,12), 
             blocklinewidth=3, blockcolour="blue")}}
\keyword{design}
\keyword{plot}