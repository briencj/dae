\name{designGGPlot}
\alias{designGGPlot}
\title{Plots labels on two-way grids of coloured cells using \code{ggplot2} to represent an experimental design}
\description{Plots the \code{labels} in a grid of cells specified by  
       \code{row.factors} and \code{column.factors}. The cells can be coloured by the values of 
       the column specified by \code{column.name} and can be divided into facets by 
       specifying multiple row and or column factors.}
\usage{
designGGPlot(design, labels = NULL, label.size = NULL, 
             row.factors = "Rows", column.factors = "Columns", 
             scales.free = "free", facetstrips.switch = NULL, 
             facetstrips.placement = "inside", 
             cellfillcolour.column = NULL, colour.values = NULL, 
             cellalpha = 1, celllinetype = "solid", celllinesize = 0.5, 
             celllinecolour = "black", cellheight = 1, cellwidth = 1,
             reverse.x = FALSE, reverse.y = TRUE, x.axis.position = "top", 
             xlab, ylab, title, labeller = label_both, 
             title.size = 15, axis.text.size = 15, 
             blocksequence = FALSE, blockdefinition = NULL, 
             blocklinecolour = "blue", blocklinesize = 2, 
             printPlot = TRUE, ggplotFuncs = NULL, ...)}
\arguments{
 \item{design}{A \code{\link{data.frame}} containing \code{labels}, \code{column.factors}, 
             \code{row.factors} and, if specified, \code{colour.column}.}
 \item{labels}{A \code{character} giving the name of the column in 
               \code{data} containing the labels to be plotted on the grid. 
               If \code{labels} is \code{NULL}, no labels are added.}
 \item{label.size}{A \code{numeric} giving the size of the labels.}
 \item{row.factors}{A \code{character} giving the names of the \code{\link{factor}}s 
               (or \code{\link{numeric}}s) in \code{data} that index the rows of the 
               plot grid used to represent the design. If there is more than one name, 
               then \code{facet_grid} is used to \code{facet} the plot in the \code{y} 
               direction, based on all but the last name. The \code{\link{factor}} 
               corresponding to the last name will index the rows in each \code{facet}.}
 \item{column.factors}{A \code{character} giving the names of the \code{\link{factor}}s 
               (or \code{\link{numeric}}s) in \code{data} that index the columns of the 
               plot grid used to represent the design. If there is more than one name, 
               then \code{facet_grid} is used to \code{facet} the plot in the \code{x} 
               direction, based on all but the last name. The \code{\link{factor}} 
               corresponding to the last name will index the columns in each \code{facet}.}
 \item{scales.free}{When plots are facetted, a \code{character} specifying whether scales 
              are shared across all facets (\code{fixed}), or vary across rows 
              (\code{free_x}), columns (\code{free_y}), or both rows and columns 
              (the default, \code{free}).  
              The \code{free_x}, \code{free_y} and \code{free} options may not work when 
              the plot grid is indexed using \code{\link{numeric}}s.}
 \item{facetstrips.switch}{When plots are facetted, the strip text are displayed on the 
              top and right of the plot by default.  If \code{facetstrips.switch} is "x", 
              the top strip text will be switched to the bottom. If "y", the right-hand 
              side labels will be switched to the left. The argument can also be set to 
              "both". The argument \code{facetstrips.placement} can be 
              used to change the relationship between the \code{strip text} and the 
              \code{axis.text} and the \code{axis.title}.}
 \item{facetstrips.placement}{A \code{character} nominating where the strip is to be placed 
              with respect to axes text and titles, either "inside" both text and titles, 
              "outside.text" or "outside.title". This argument is important only when axes 
              and strips are on the same side of the plot. When this occurs, the default is 
              to place them inside the axis text. \bold{Note:} This argument must be specified 
              only once in the constructon of the plot and after every other aspect of the plot 
              has been finalized. In particular, if \code{\link{designBlocksGGPlot}} is to called  
              after \code{designGGPlot} to add block boundaries, then \code{facetstrips.placement} 
              should be specified in the call to \code{\link{designBlocksGGPlot}}, not in the 
              call to \code{designGGPlot}.}
 \item{reverse.x}{A \code{logical} which, if true, causes the order of values on the 
              x-axis to be reversed, the natural order being to increase from left to right.}
 \item{reverse.y}{A \code{logical} which, if true, causes the order of values on the 
              y-axis to be reversed, the natural order being to increase from bottom to top.}
 \item{x.axis.position}{A \code{character} giving the position of the x-axis; can be \code{top} 
               or \code{bottom}.}
 \item{cellfillcolour.column}{A \code{character} giving the name of the column in 
               \code{data} that is to be used to vary the colour the used to fill a cell.}
 \item{colour.values}{A \code{character} giving the name or names of the colours to be used 
               in filling the cell. If cellfillcolour.column is not \code{NULL} then the 
               number of colours sepcified needs to match the number of unique values in 
               the cellfillcolour.column.}
  \item{cellalpha}{A \code{numeric} specifying the degree of transparency to be used in cell fill. 
               It is a ratio in which the denominator specifies the number of points (or 
               lines) that must be overplotted to give a solid cover.}
 \item{celllinetype}{A \code{numeric} or \code{character} giving the type of line for the 
               cell border. An integer or name: 0 = blank, 1 = solid, 2 = dashed, 
               3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash. For more information 
               see \code{vignette("ggplot2-specs")}.}
 \item{celllinesize}{A \code{numeric} specifying the size of the line in mm.}
 \item{celllinecolour}{A \code{character} giving the name of the colour to use for the cell 
               outline.}
 \item{cellheight}{A \code{numeric} specifying the height of a cell.}
 \item{cellwidth}{A \code{numeric} specifying the width of of a cell.}
 \item{xlab}{\code{Label} for the x-axis. By default it is the last name in the
             \code{column.factors}.}
 \item{ylab}{\code{Label} for the y-axis. By default it is the last name in the
             \code{row.factors}.}
 \item{title}{Title for plot window. By default it is "Plot of labels".}
 \item{labeller}{A \code{function} for specifying the formatting of the strip labels of 
                 the facet grids used when there is more than one \code{row.factors} or 
                 \code{column.factors}. (See \code{\link{labellers}} for \code{ggplot2}.)}
 \item{title.size}{A \code{numeric} giving the size for all titles: \code{xlab}, 
              \code{ylab} and \code{title}.}
 \item{axis.text.size}{A \code{numeric} giving the size for tick labels.}
 \item{blocksequence}{A \code{\link{logical}} that determines whether block numbers 
              are repetitions or sequences of block numbers.}
 \item{blockdefinition}{A \code{\link{matrix}} of block sizes:
      \itemize{
        \item if there is only one row, then the first element is interpreted as 
               the no. rows in each block and blocks with this number of rows are to be 
               repeated  across the rows of the design.
        \item if there is more than one row, then each row of the matrix specifies a block, 
               with the sequence of rows in the matrix specifying a corresponding 
                sequence of blocks down the rows of the design.}
       Similarly, a single value for a column specifies a repetition of blocks of that size 
       across the columns of the design, while several column values specifies a 
       sequence of blocks across the columns of the size specified.}
 \item{blocklinecolour}{A \code{\link{character}} string specifying the colour of the 
       block boundary.
    
       See also the \code{scale_colour_*} functions or \code{Colour specification} under the
       \code{\link{par}} function.}
 \item{blocklinesize}{A \code{\link{numeric}} giving the width of the block boundary to be
        plotted.}
 \item{printPlot}{A \code{\link{logical}} indicating whether to print the plot produced.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{\dots}{Other arguments that are passed down to the \code{\link{geom_text}} call that plots the \code{labels}.}
}
\value{An object of class "\code{ggplot}", which can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{designBlocksGGPlot}}, \code{\link{fac.combine}} in package \pkg{dae}, 
         \code{\link{designPlot}}.}
\examples{
#### Plot a randomized complete block design
Treatments <- factor(rep(1:6, times = 5))
RCBD.lay <- designRandomize(allocated = Treatments,
                            recipient = list(Blocks = 5, Units = 6),
                            nested.recipients = list(Units = "Blocks"),
                            seed = 74111)
designGGPlot(RCBD.lay, labels = "Treatments", size = 5, 
             row.factors = "Blocks", column.factors = "Units", 
             blockdefinition = cbind(1,5))
             
## Plot without labels
designGGPlot(RCBD.lay, cellfillcolour.column = "Treatments", 
             row.factors = "Blocks", column.factors = "Units", 
             colour.values = c("lightblue","lightcoral","lightgoldenrod",
                               "lightgreen","lightgrey", "lightpink"), 
             blockdefinition = cbind(1,6))

             
#### Plot a lattice square design
data(LatticeSquare_t49.des)
designGGPlot(LatticeSquare_t49.des, labels = "Lines", size = 5, 
             row.factors = c("Intervals", "Runs"), column.factors = "Times", 
             blockdefinition = cbind(7,7))
}
\keyword{aplot}
\keyword{hplot}
\keyword{design}
