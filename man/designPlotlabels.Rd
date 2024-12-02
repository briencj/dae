\name{designPlotlabels}
\alias{designPlotlabels}
\title{Plots labels on a two-way grid using \pkg{ggplot2}}
\description{Plots the \code{labels} in a grid specified by  
       \code{grid.x}and \code{grid.y}. The labels can be coloured by the values of 
       the column specified by \code{column.name}.}
\usage{
designPlotlabels(data, labels, grid.x = "Columns", grid.y = "Rows", 
                 colour.column=NULL, colour.values=NULL, 
                 reverse.x = FALSE, reverse.y = TRUE, 
                 xlab, ylab, title, printPlot = TRUE, ggplotFuncs = NULL, ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing \code{labels}, \code{grid.x}, 
             \code{grid.y} and, if specified, \code{colour.column}.}
 \item{labels}{A \code{character} giving the name of the column in 
               \code{data} containing the labels to be plotted on the grid.}
 \item{grid.x}{A \code{character} giving the name of the column in \code{data} that 
               specifies the x-coordinates of the plot grid.}
 \item{grid.y}{A \code{character} giving the name of the column in \code{data} that 
               specifies the y-coordinates of the plot grid.}
 \item{reverse.x}{A \code{logical} which, if true, causes the order of values on the 
              x-axis to be reversed.}
 \item{reverse.y}{A \code{logical} which, if true, causes the order of values on the 
              y-axis to be reversed.}
 \item{colour.column}{A \code{character} giving the name of the column in \code{data} that 
               is to be used to colour the values plotted on the grid.}
 \item{colour.values}{A \code{character} giving the name of the column in \code{data} that 
               is to be used to colour the values plotted on the grid.}
 \item{xlab}{\code{Label} for the x-axis. By default it is the name of the
             \code{grid.x}.}
 \item{ylab}{\code{Label} for the y-axis. By default it is the name of the
             \code{grid.y}.}
 \item{title}{Title for plot window. By default it is "Plot of labels".}
 \item{printPlot}{A \code{\link{logical}} indicating whether to print the plot.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \pkg{ggplot2} function. 
                It is created by calling the \code{\link{list}} function with 
                a \pkg{ggplot2} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{\dots}{Other arguments that are passed down to the \code{\link[ggplot2]{geom_text}} call that plots the \code{labels}.}
}
\value{An object of class "\code{ggplot}", which can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}} in package \pkg{dae}, \code{\link{designPlot}}.}
\examples{
Treatments <- factor(rep(1:6, times = 5))
RCBD.lay <- designRandomize(allocated = Treatments,
                            recipient = list(Blocks = 5, Units = 6),
                            nested.recipients = list(Units = "Blocks"),
                            seed = 74111)
designPlotlabels(RCBD.lay, labels = "Treatments", 
                 grid.x = "Units", grid.y = "Blocks",
                 colour.column = "Treatments", size = 5)
}
\keyword{aplot}
\keyword{hplot}
\keyword{design}
