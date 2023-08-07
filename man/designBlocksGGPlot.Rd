\name{designBlocksGGPlot}
\alias{designBlocksGGPlot}
\title{Adds block boundaries to a plot produced by \code{\link{designGGPlot}}.}
\description{This function adds block boundaries to a plot produced by \code{\link{designGGPlot}}.
             It allows control of the starting unit, through originrow and origincolumn, 
             and the number of rows (nrows) and columns (ncolumns) from the starting unit 
             that the blocks to be plotted are to cover.}
\usage{
designBlocksGGPlot(ggplot.obj, blockdefinition = NULL, blocksequence = FALSE, 
                   originrow= 0, origincolumn = 0, nrows, ncolumns, 
                   blocklinecolour = "blue", blocklinesize = 2, 
                   facetstrips.placement = "inside", 
                   printPlot = TRUE)}
\arguments{
 \item{ggplot.obj}{An object produced by \code{ggplot}.}
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
 \item{blocksequence}{A \code{\link{logical}} that determines whether block numbers are 
         repetitions or sequences of block numbers.}
 \item{originrow}{A \code{\link{numeric}} speccifying the row after which the plotting 
               of block boundaries is to start.}
 \item{origincolumn}{A \code{\link{numeric}} speccifying the column after which the plotting 
               of block boundaries is to start.}
 \item{nrows}{A \code{\link{numeric}} the number of rows (nrows), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{ncolumns}{A \code{\link{numeric}} the number of columns (ncolumns), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{blocklinecolour}{A \code{\link{character}} string specifying the colour of the block boundary.
  
             See \code{Colour specification} under the \code{\link{par}} function.}
 \item{blocklinesize}{A \code{\link{numeric}} giving the width of the block boundary to be plotted.}
 \item{facetstrips.placement}{A \code{character} nominating where the strip is to be placed 
              with respect to axes text and titles, either "inside" both text and titles, 
              "outside.text" or "outside.title". This argument is important only when axes 
              and strips are on the same side of the plot. When this occurs, the default is 
              to place them inside the axis text. \bold{Note:} This argument must be specified 
              only once in the constructon of the plot and after every other aspect of the plot 
              has been finalized. In particular, if \code{designBlocksGGPlot} is to called  
              after \code{designGGPlot} to add block boundaries, then \code{facetstrips.placement} 
              should be specified in the last call to \code{designBlocksGGPlot}, not in the  
              call to \code{designGGPlot}.}
 \item{printPlot}{A \code{\link{logical}} indicating whether to print the plot after adding the block 
             boundaries.}
}
\value{An object of class "\code{ggplot}", formed by adding to the input \code{ggplot.obj} and 
      which can be plotted using \code{print}.}


\author{Chris Brien}
\seealso{\code{\link{designGGPlot}}, \code{\link{par}}, \code{DiGGer}}
\examples{
## Construct a randomized layout for the split-unit design described by 
## Brien et al. (2011, Section 5)
split.sys <- cbind(fac.gen(list(Months = 4, Athletes = 3, Tests = 3)),
                   fac.gen(list(Intensities = LETTERS[1:3], Surfaces = 3), 
                           times = 4))
split.lay <- designRandomize(allocated = split.sys[c("Intensities", "Surfaces")],
                             recipient = split.sys[c("Months", "Athletes", "Tests")], 
                             nested.recipients = list(Athletes = "Months", 
                                                      Tests = c("Months", "Athletes")),
                             seed = 2598)
## Plot the design
cell.colours <- c("lightblue","lightcoral","lightgoldenrod","lightgreen","lightgrey",
                  "lightpink","lightsalmon","lightcyan","lightyellow","lightseagreen")

split.lay <- within(split.lay, 
                    Treatments <- fac.combine(list(Intensities, Surfaces), 
                                              combine.levels = TRUE))
plt <- designGGPlot(split.lay, labels = "Treatments", 
                    row.factors = "Tests", column.factors = c("Months", "Athletes"),
                    colour.values = cell.colours[1:9], size = 6, 
                    blockdefinition = rbind(c(3,1)), blocklinecolour = "darkgreen",
                    printPlot = FALSE)
#Add Month boundaries
designBlocksGGPlot(plt, nrows = 3, ncolumns = 3, blockdefinition = rbind(c(3,3)))



#### A layout for a growth cabinet experiment that allows for edge effects
data(Cabinet1.des)
plt <- designGGPlot(Cabinet1.des, labels = "Combinations", cellalpha = 0.75,
                    title = "Lines and Harvests allocation for Cabinet 1", 
                    printPlot = FALSE)

## Plot Mainplot boundaries
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(4,16), originrow= 1 , 
                          blocklinecolour = "green", nrows = 9, ncolumns = 16, 
                          printPlot = FALSE)
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(1,4), 
                          blocklinecolour = "green", nrows = 1, ncolumns = 16, 
                          printPlot = FALSE)
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(1,4), originrow= 9, 
                          blocklinecolour = "green", nrows = 10, ncolumns = 16, 
                          printPlot = FALSE)
## Plot all 4 block boundaries            
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(8,5,5,4), blocksequence = TRUE, 
                          origincolumn = 1, originrow= 1, 
                          blocklinecolour = "blue", nrows = 9, ncolumns = 15, 
                          printPlot = FALSE)
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(10,16), 
                          blocklinecolour = "blue", nrows = 10, ncolumns = 16, 
                          printPlot = FALSE)
## Plot border and internal block boundaries only
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(8,14), origincolumn = 1, originrow= 1, 
                          blocklinecolour = "blue", nrows = 9, ncolumns = 15, 
                          printPlot = FALSE)
plt <- designBlocksGGPlot(plt, blockdefinition = cbind(10,16), 
                          blocklinecolour = "blue", nrows = 10, ncolumns = 16)
}
\source{
Brien, C.J., Harch, B.D., Correll, R.L., and Bailey, R.A. (2011) Multiphase experiments with at least one later laboratory phase. I. Orthogonal designs. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 16:422-450.}

\keyword{design}
\keyword{plot}