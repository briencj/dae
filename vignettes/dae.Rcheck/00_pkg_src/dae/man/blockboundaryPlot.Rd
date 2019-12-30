\name{blockboundaryPlot}
\alias{blockboundaryPlot}
\title{This function plots a block boundary on a plot produced by \code{\link{designPlot}}.}
\description{This function plots a block boundary on a plot produced by \code{\link{designPlot}}.
             It allows control of the starting unit, through rstart and cstart, 
             and the number of rows (nrows) and columns (ncolumns) from the starting unit 
             that the blocks to be plotted are to cover.}
\usage{
blockboundaryPlot(blockdefinition = NULL, blocksequence = FALSE, 
                  rstart= 0, cstart = 0, nrows, ncolumns, 
                  blocklinecolour = 1, blocklinewidth = 2)}
\arguments{
 \item{blockdefinition}{A \code{\link{matrix}} of block sizes:
    \itemize{
       \item if there is only one row, then the first element is interpreted as the no. rows in 
         each block and blocks with this number of rows are to be repeated  across the rows of 
         the design.
       \item if there is more than one row, then each row of the matrix specifies a block, 
         with the sequence of rows in the matrix specifying a corresponding 
         sequence of blocks down the rows of the design.}
    Similarly, a single value for a column specifies a repetition of blocks of that size 
         across the columns of the design, while several column values specifies a 
         sequence of blocks across the columns of the size specified.}
 \item{blocksequence}{A \code{\link{logical}} that determines whether block numbers are repetitions 
              or sequences of block numbers.}
 \item{rstart}{A \code{\link{numeric}} speccifying the row after which the plotting 
               of block boundaries is to start.}
 \item{cstart}{A \code{\link{numeric}} speccifying the column after which the plotting 
               of block boundaries is to start.}
 \item{nrows}{A \code{\link{numeric}} the number of rows (nrows), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{ncolumns}{A \code{\link{numeric}} the number of columns (ncolumns), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{blocklinecolour}{A \code{\link{character}} string specifying the colour of the block boundary.
  
             See \code{Colour specification} under the \code{\link{par}} function.}
 \item{blocklinewidth}{A \code{\link{numeric}} giving the width of the block boundary to be plotted.}
}
\value{no values are returned, but modifications are made to the currently active plot.}


\author{Chris Brien}
\seealso{\code{\link{designPlot}}, \code{\link{par}}, \code{DiGGer}}
\examples{\dontrun{
    SPL.Lines.mat <- matrix(as.numfac(Lines), ncol=16, byrow=T)
    colnames(SPL.Lines.mat) <- 1:16
    rownames(SPL.Lines.mat) <- 1:10
    SPL.Lines.mat <- SPL.Lines.mat[10:1, 1:16]
    designPlot(SPL.Lines.mat, labels=1:10, new=TRUE,
               rtitle="Rows",ctitle="Columns", 
               chardivisor=3, rcellpropn = 1, ccellpropn=1,
               plotcellboundary = TRUE)
    #Plot Mainplot boundaries
    blockboundaryPlot(blockdefinition = cbind(4,16), rstart = 1, 
                      blocklinewidth = 3, blockcolour = "green", 
                      nrows = 9, ncolumns = 16)
    blockboundaryPlot(blockdefinition = cbind(1,4), 
                      blocklinewidth = 3, blockcolour = "green", 
                      nrows = 1, ncolumns = 16)
    blockboundaryPlot(blockdefinition = cbind(1,4), rstart= 9, nrows = 10, ncolumns = 16, 
                      blocklinewidth = 3, blockcolour = "green")
    #Plot all 4 block boundaries            
    blockboundaryPlot(blockdefinition = cbind(8,5,5,4), blocksequence=T, 
                      cstart = 1, rstart= 1, nrows = 9, ncolumns = 15, 
                      blocklinewidth = 3,blockcolour = "blue")
    blockboundaryPlot(blockdefinition = cbind(10,16), blocklinewidth=3, blockcolour="blue", 
                      nrows=10, ncolumns=16)
    #Plot border and internal block boundaries only
    blockboundaryPlot(blockdefinition = cbind(8,14), cstart = 1, rstart= 1, 
                      nrows = 9, ncolumns =  15,
                      blocklinewidth = 3, blockcolour = "blue")
    blockboundaryPlot(blockdefinition = cbind(10,16), 
                      blocklinewidth = 3, blockcolour = "blue", 
                      nrows = 10, ncolumns = 16)}}

\keyword{design}
\keyword{plot}