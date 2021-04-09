getLinePosns <- function(axis.posns, endspace = 0.5)
{
  endsonly <- FALSE
  if (length(axis.posns) <= 2) endsonly <- TRUE
  axis.posns <- sort(unique(axis.posns))
  half.diffs <- diff(axis.posns)/2
  if (endsonly)
    line.posns <- c(axis.posns[1]-endspace, 
                    axis.posns[2]+endspace)
  else
    line.posns <- c(axis.posns[1]-endspace, 
                    axis.posns[1:(length(axis.posns)-1)] + half.diffs, 
                    axis.posns[length(axis.posns)]+endspace)
  return(line.posns)
}

"designGGPlot" <- function(design, labels = NULL, label.size = NULL,
                           row.factors = "Rows", column.factors = "Columns", scales.free = "free", 
                           cellfillcolour.column=NULL, colour.values=NULL, cellalpha = 1, 
                           celllinetype = "solid", celllinesize = 0.5, celllinecolour = "black",
                           cellheight = 1, cellwidth = 1, 
                           reverse.x = FALSE, reverse.y = TRUE, x.axis.position = "top", 
                           xlab, ylab, title, labeller = label_both, 
                           title.size = 15, axis.text.size = 15, 
                           blocksequence = FALSE, blockdefinition = NULL, 
                           blocklinecolour = "blue", blocklinesize = 2, 
                           printPlot = TRUE, ggplotFuncs = NULL, ...)
  ## Function that uses ggplot to plot labels on a grid
{
  
  opts <- c("top", "bottom")
  x.axis.position <- opts[check.arg.values(x.axis.position, opts)]
  
  #Check for multiple factors in either direction
  #Rows
  if (length(row.factors) == 1) {
    grid.y <- row.factors
    facet.y <- NULL
    
  }
  else 
  {
    grid.y <- row.factors[length(row.factors)]
    facet.y <- row.factors[-length(row.factors)]
    if (reverse.y)
      for (fac in facet.y)
        design[fac] <- factor(design[[fac]], levels = rev(levels(design[[fac]])))
    facet.y <- paste0("vars(", paste(facet.y, collapse = ","), ")")
  }
  #Columns
  if (length(column.factors) == 1) {
    grid.x <- column.factors
    facet.x <- NULL
    
  }
  else 
  {
    grid.x <- column.factors[length(column.factors)]
    facet.x <- column.factors[-length(column.factors)]
    if (reverse.x)
      for (fac in facet.x)
        design[fac] <- factor(design[[fac]], levels = rev(levels(design[[fac]])))
    facet.x <- paste0("vars(", paste(facet.x, collapse = ","), ")")  
  }
  
  if (missing(xlab)) xlab <- grid.x
  if (missing(ylab)) ylab <- grid.y
  if (missing(title)) title <- paste("Plot of",labels,sep = " ")

  #Set up the plot
  plt <- ggplot(data = design, aes_string(x = grid.x, y = grid.y)) +
    labs(x = xlab, y = ylab, title = title) + 
    theme(panel.background = element_blank(),
          legend.position = "none",
          title = element_text(size = title.size, face = "bold"),
          axis.text = element_text(size = axis.text.size, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size=title.size, face="bold"))
  
  if (!(is.null(colour.values)))
    plt <- plt + scale_fill_manual(values = colour.values)
  
  #Check that one of labels and cellfillcolour.column
  if (is.null(labels) && is.null(cellfillcolour.column))
    stop("At least one of labels and cellfillcolour.column must be set")
  #Create tiles
  if (is.null((cellfillcolour.column)))
    plt <- plt +  geom_tile(aes_string(fill = labels), 
                            colour = celllinecolour, alpha = cellalpha, 
                            linetype = celllinetype, size = celllinesize, 
                            height = cellheight, width = cellwidth)
  else
    plt <- plt +  geom_tile(aes_string(fill = cellfillcolour.column), 
                            colour = celllinecolour, alpha = cellalpha, 
                            linetype = celllinetype, size = celllinesize, 
                            height = cellheight, width = cellwidth)
  
  #Add labels, if specified
  if (!is.null(labels))
  {
    if (!is.null(label.size))
      plt <- plt + geom_text(aes_string(label = labels), size = label.size, 
                             fontface = "bold", ...)
    else
      plt <- plt + geom_text(aes_string(label = labels), fontface = "bold", ...)
  }
  
  #Set up y scale
  if (inherits(design[[grid.y]], what = "factor"))
  {
    nrows <- length(levels(design[[grid.y]]))
    if (reverse.y)
      plt <- plt + scale_y_discrete(limits = rev, expand = c(0,0))
#    plt <- plt + scale_y_discrete(limits = rev(levels(design[[grid.y]])), expand = c(0,0))
    else
      plt <- plt + scale_y_discrete(expand = c(0,0))
  }
  else
  {
    rows <- sort(unique(design[[grid.y]]))
    nrows <- length(rows)
    row.posns <- getLinePosns(rows)
    if (reverse.y)
      plt <- plt + scale_y_reverse(limits = c(row.posns[c(1,length(row.posns))]), 
                                      expand = c(0,0)) 
    else
      plt <- plt + scale_y_continuous(limits = c(row.posns[c(1,length(row.posns))]), 
                                      expand = c(0,0)) 
  }
  
  #Set up x scale
  if (inherits(design[[grid.x]], what = "factor"))
  {
    ncolumns <- length(levels(design[[grid.x]]))
    if (reverse.x)
      plt <- plt + scale_x_discrete(limits = rev, expand = c(0,0), position = x.axis.position)
#    plt <- plt + scale_x_discrete(limits = rev(levels(design[[grid.x]])), expand = c(0,0), 
#                                    position = x.axis.position)
    else
      plt <- plt + scale_x_discrete(expand = c(0,0), position = x.axis.position)
  }
  else
  {
    columns <- sort(unique(design[[grid.x]]))
    ncolumns <- length(columns)
    col.posns <- getLinePosns(columns)
    if (reverse.x)
      plt <- plt + scale_x_reverse(limits = c(col.posns[c(1,length(col.posns))]), 
                                      expand = c(0,0), position = x.axis.position)
    else
      plt <- plt + scale_x_continuous(limits = c(col.posns[c(1,length(col.posns))]), 
                                      expand = c(0,0), position = x.axis.position)
  }

  #Set up facetting if needed
  # if (!is.null(facet.x))
  #   facet.x <- paste0("vars(", paste(facet.x, collapse = ","), ")")
  # if (!is.null(facet.y))
  #   facet.y <- paste0("vars(", paste(facet.y, collapse = ","), ")")
  if (!is.null(facet.x)) 
  {
    if (!is.null(facet.y))
      plt <- plt + facet_grid(rows = eval(parse(text=facet.y)), 
                              cols = eval(parse(text=facet.x)), 
                              labeller = labeller, scales = scales.free, as.table = FALSE)
    else
      plt <- plt + facet_grid(cols = eval(parse(text=facet.x)), 
                              labeller = labeller, scales = scales.free, as.table = FALSE)
  } else
  {
    if (!is.null(facet.y))      
      plt <- plt + facet_grid(rows = eval(parse(text=facet.y)), 
                              labeller = labeller, scales = scales.free, as.table = FALSE)
  }
  
  if (!is.null(ggplotFuncs))
  {
    for (f in ggplotFuncs)
      plt <- plt + f
  }
  
  if (!is.null(blockdefinition))
    plt <- designBlocksGGPlot(plt, nrows = nrows, ncolumns = ncolumns, 
                              blocksequence = blocksequence, blockdefinition = blockdefinition, 
                              blocklinecolour = blocklinecolour, blocklinesize = blocklinesize,
                              printPlot = printPlot)
  else
  { 
    if (printPlot)
      print(plt)
    invisible(plt)
  }
}


"plotarectangle" <- function(plt, xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinesize,blocklinecolour)
{ 
  ncimod <- nci
  nrimod <- nri
  if (xoff + nci > ncolumns) 
  { 
    ncimod <- ncolumns - xoff
  }
  if (yoff + nri > nrows) 
  { 
    nrimod <- nrows - yoff
  }
  lines.dat <- data.frame(x = xi + xoff + c(1, 1, ncimod, ncimod, 1),
                          y = yi + yoff + c(nrimod, 1, 1, nrimod, nrimod))
  plt <- plt + geom_path(data = lines.dat, mapping = aes_string(x="x",y="y"), 
                         colour = blocklinecolour, size = blocklinesize)
  invisible(plt)
}


"designBlocksGGPlot" <- function(ggplot.obj, blockdefinition = NULL, blocksequence = FALSE, 
                                 originrow = 0, origincolumn = 0, nrows, ncolumns, 
                                 blocklinecolour = "blue", blocklinesize = 2, printPlot = TRUE)
{
  #This function is a modified version of code extracted from moddes.plot
  #It allows one to set the rectangle for plotting using 
  #blockdefinition is a matrix of block sizes:
  #if there is only one row, then it is interpreted as the no. rows to be repeated
  #     for a sequence of blocks whose size is specified all but the first element in the row.
  #if there is more than one row, then each row of the matrix specifies a block,
  #     with the sequence of rows in the matrix specifying a corresponding
  #     sequence of blocks down the rows of the design.
  #Similarly, a single value for a column specifies a repetition of blocks of that size
  #     across the columns of the design, while several column values specifies a
  #     sequence of blocks across the columns of the size specified.
  if (!is.null(blockdefinition))
  { 
    rstart <- originrow
    cstart <- origincolumn
    dims <- dim(blockdefinition)
    xi <- c(-0.5, -0.5, 0.5, 0.5, -0.5)
    yi <- c(0.5, -0.5, -0.5, 0.5, 0.5)
    if (!blocksequence) #blockdefinition interpreted as repetitions of blocks of specified size
    { 
      for (i in seq(dims[1]))
      { 
        nri <- blockdefinition[i, 1]
        nci <- blockdefinition[i, 2]
        for (j in seq(ceiling((nrows - rstart)/nri)))
        { 
          for (k in seq(ceiling((ncolumns - cstart)/nci)))
          { 
            xoff <- nci * (k - 1) + cstart
            yoff <- nri * (j - 1) + rstart
            ggplot.obj <- plotarectangle(ggplot.obj, xi, yi, xoff, yoff, 
                                         nrows, ncolumns, nri, nci,
                                         blocklinesize, blocklinecolour)
          }
        }
      }
      # rows <- seq(rstart, nrows, by = nri)
      # cols <- seq(cstart, ncolumns, by = nci)
      # if (rows[length(rows)] != nrows)
      #   rows <- c(rows, nrow)
      # if (cols[length(cols)] != ncolumns)
      #   cols <- c(cols, ncolumns)
      # ggplot.obj <- ggplot.obj +  
      #   geom_hline(yintercept = rows+0.5, colour = blocklinecolour, size = blocklinesize) +
      #   geom_vline(xintercept = cols+0.5, colour = blocklinecolour, size = blocklinesize)
    }
    else #blockdefinition interpreted as a sequence of block specification
    { 
      if (dims[1] > 1) #multiple rows
      { 
        yoff <- rstart
        for (k in seq(dims[1])) { 
          if (dims[2] > 2) #multiple columns
          { 
            xoff <- cstart
            nri <- blockdefinition[k, 1]
            for (i in seq(2,dims[2])) {
              nci <- blockdefinition[k, i]
              ggplot.obj <- plotarectangle(ggplot.obj, xi, yi, xoff, yoff, 
                                           nrows, ncolumns, nri, nci,
                                           blocklinesize, blocklinecolour)
              xoff <- xoff + nci
            }
          }
          else  #single column specifier
          { 
            nri <- blockdefinition[k, 1]
            nci <- blockdefinition[k, 2]
            for (j in seq(ceiling((ncolumns - cstart)/nci)))
            { 
              xoff <- nci * (j - 1) + cstart
              ggplot.obj <- plotarectangle(ggplot.obj, xi, yi, xoff, yoff, 
                                           nrows, ncolumns, nri, nci,
                                           blocklinesize, blocklinecolour)
            }
          }
          yoff <- yoff + nri
        }
      }
      else  #only one row in matrix
      { 
        if (dims[2] > 2) #multiple columns
        { 
          xoff <- cstart
          nri <- blockdefinition[1, 1]
          for (i in seq(2,dims[2]))
          { 
            nci <- blockdefinition[1, i]
            for (j in seq(ceiling(nrows/nri - rstart)))
            { 
              yoff <- nri * (j - 1) + rstart
              ggplot.obj <- plotarectangle(ggplot.obj, xi, yi, xoff, yoff, 
                                           nrows, ncolumns, nri, nci,
                                           blocklinesize, blocklinecolour)
            }
            xoff <- xoff + nci
          }
        }
        else #only one row and one column specified in the matrix
        { 
          nri <- blockdefinition[1, 1]
          nci <- blockdefinition[1, 2]
          for (j in seq(ceiling((nrows - rstart)/nri)))
          { 
            for (k in seq(ceiling((ncolumns - cstart)/nci)))
            { 
              xoff <- nci * (k - 1) + cstart
              yoff <- nri * (j - 1) + rstart
              ggplot.obj <- plotarectangle(ggplot.obj, xi, yi, xoff, yoff, 
                                           nrows, ncolumns, nri, nci,
                                           blocklinesize, blocklinecolour)
            }
          }
        }
      }
    }
  }
  
  if (printPlot)
    print(ggplot.obj)
  invisible(ggplot.obj)
  
  invisible(ggplot.obj)
}
