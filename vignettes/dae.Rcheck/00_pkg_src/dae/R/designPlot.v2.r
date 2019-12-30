"designPlotlabels" <- function(data, labels, grid.x = "Columns", grid.y = "Rows", 
                               colour.column=NULL, colour.values=NULL, 
                               reverse.x = FALSE, reverse.y = TRUE, 
                               xlab, ylab, title, printPlot = TRUE, ggplotFuncs = NULL, ...)
  ## Function that uses ggplot to plot labels on a grid
{
  warning("designPlotlabels has been superseded by designGGPlot and will be deprecated in a future version")
  if (missing(xlab)) xlab <- grid.x
  if (missing(ylab)) ylab <- grid.y
  if (missing(title)) title <- paste("Plot of",labels,sep = " ")
  
  plt <- ggplot(data = data, aes_string(x = grid.x, y = grid.y, label = labels)) +
    labs(x = xlab, y = ylab, title = title)

  if (reverse.x)
  {
    if (inherits(data[[grid.x]], what = "factor"))
      plt <- plt + scale_x_discrete(limits = rev(levels(data[[grid.x]])))
    else
      plt <- plt + scale_x_reverse()
  }
  
  if (reverse.y)
  {
    if (inherits(data[[grid.y]], what = "factor"))
      plt <- plt + scale_y_discrete(limits = rev(levels(data[[grid.y]])))
    else
      plt <- plt + scale_y_reverse()
  }
  
  if (!is.null(colour.column)) 
    if (labels == colour.column)
      plt <- plt + theme(legend.position = "none")
  
  if (is.null((colour.column)))
    plt <- plt + geom_text(aes_string(), ...)
  else
    plt <- plt + geom_text(aes_string(colour = colour.column), ...)
  
  if (!(is.null(colour.values)))
    plt <- plt + scale_colour_manual(values = colour.values)
  
  if (!is.null(ggplotFuncs))
  {
    for (f in ggplotFuncs)
      plt <- plt + f
  }
  
  if (printPlot)
    print(plt)
  
  invisible(plt)
}


"plotablock" <- function(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
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
  lines(xi + xoff + c(1, 1, ncimod, ncimod, 1),
        yi - yoff - c(1, nrimod, nrimod, 1, 1), lwd = blocklinewidth,
        col = blocklinecolour)
  invisible()
}


"blockboundaryPlot" <- function(blockdefinition = NULL, blocksequence = FALSE, rstart= 0, cstart = 0, 
                                nrows, ncolumns, blocklinecolour = 1, blocklinewidth = 2)
  #This function is a modified version of code extracted from moddes.plot
  #It allows one to set the rectangle for plotting using 
{
  #blockdefinition is a matrix of block sizes:
  #if there is only one row, then it is interpreted as the no. rows to be repeated
  #     for a sequence of blocks whose size is specified all but the first element in the row.
  #if there is more than one row, then each row of the matrix specifies a block,
  #     with the sequence of rows in the matrix specifying a corresponding
  #     sequence of blocks down the rows of the design.
  #Similarly, a single value for a column specifies a repetition of blocks of that size
  #     across the columns of the design, while several column values specifies a
  #     sequence of blocks across the columns of the size specified.
  warning("blockboundaryPlot has been superseded by designBlocksGGPlot and will be deprecated in a future version")
  
  if (!is.null(blockdefinition))
  { dims <- dim(blockdefinition)
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
          plotablock(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
        }
      }
    }
  }
  else #blockdefinition interpreted as a sequence of block specification
  { 
    if (dims[1] > 1) #multiple rows
    { 
      yoff <- rstart
      for (k in seq(dims[1]))
      { 
        if (dims[2] > 2) #multiple columns
        { 
          xoff <- cstart
          nri <- blockdefinition[k, 1]
          for (i in seq(2,dims[2]))
          { nci <- blockdefinition[k, i]
          plotablock(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
          xoff <- xoff + nci
          }
        }
        else  #single column specifier
        { nri <- blockdefinition[k, 1]
        nci <- blockdefinition[k, 2]
        for (j in seq(ceiling((ncolumns - cstart)/nci)))
        { 
          xoff <- nci * (j - 1) + cstart
          plotablock(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
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
            plotablock(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
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
            plotablock(xi,yi,xoff,yoff,nrows,ncolumns,nri,nci,blocklinewidth,blocklinecolour)
          }
        }
      }
    }
  }
  }
  invisible()
}

"designPlot" <- function (designMatrix, labels = NULL, altlabels = NULL, plotlabels = TRUE, 
                          rtitle = NULL, ctitle = NULL, 
                          rlabelsreverse = FALSE, clabelsreverse = FALSE, 
                          font = 1, chardivisor = 2, rchardivisor = 1, cchardivisor = 1, 
                          cellfillcolour = NA, plotcellboundary = TRUE, 
                          rcellpropn = 1, ccellpropn = 1, 
                          blocksequence = FALSE, blockdefinition = NULL, 
                          blocklinecolour = 1, blocklinewidth = 2, 
                          rotate = FALSE, new = TRUE, ...)
  #Added blocksequence on 9/5/2013
  #It determines whether block numbers are repetitions or sequences of block numbers 
{
  warning("designPlot has been superseded by designGGPlot and may be deprecated in a future version")
  if (is.null(labels)) 
    labels <- unique(as.vector(designMatrix))
  if (rcellpropn > 1 | rcellpropn <= 0 | ccellpropn > 1 | ccellpropn <=0 )
    stop("rcellpropn and ccellpropn must be positive and less than one")
  if (length(cellfillcolour) > 1 & length(cellfillcolour) < length(labels))
    stop("The number of colours must either be one or more than the number of labels")
  drow <- -1 * as.vector(row(designMatrix))
  drange <- as.vector(col(designMatrix))
  dtrt <- as.vector(designMatrix)
  nrows <- -min(drow)
  ncolumns <- max(drange)
  rowlabs <- rownames(designMatrix)
  collabs <- colnames(designMatrix)
  if (is.null(rowlabs))
    rowlabs <- paste(seq(nrows))
  if (is.null(collabs))
    collabs <- paste(seq(ncolumns))
  charot <- 0
  if (rotate) 
  {   
    dc <- dim(designMatrix)[2]
    designMatrix <- designMatrix[, rev(seq(dc))]
    designMatrix <- t(designMatrix)
    if (!is.null(blockdefinition)) 
    { 
      if (length(blockdefinition == 2)) 
        blockdefinition <- cbind(blockdefinition)
    }
    tmptitle <- ctitle
    ctitle <- rtitle
    rtitle <- tmptitle
    charot <- 90
    tmplabs <- collabs
    collabs <- rowlabs
    rowlabs <- tmplabs
  }
  csival <- min(par()$fin/c(ncolumns, nrows))/chardivisor
  if (rotate) 
  {  
    csival <- min(par()$fin/c(nrows, ncolumns))/chardivisor
  }
  cexval <- csival/par()$csi/0.7
  rcexval <- cexval*chardivisor/rchardivisor
  ccexval <- cexval*chardivisor/cchardivisor
  lineval = (max(nchar(rowlabs))+1)*rcexval*0.5
  if (new) 
  { 
    plot(range(drange) + c(-1, 1), range(drow) + c(-1, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    if (rotate) 
    { 
      if (!is.null(rtitle)) 
      { #Modification to implement rlabelsreverse - 15/12/2012
        #else option is original code
        if (rlabelsreverse)
          mtext(rowlabs, side = 2, line = 0,
                at = -seq(nrows), cex = rcexval, adj = 1, las = 1)
        else
          mtext(rev(rowlabs), side = 2, line = 0,
                at = -seq(nrows), cex = rcexval, adj = 1, las = 1)
      }
      mtext(rtitle, side = 2, line = lineval, at = -nrows/2 - 1/2,
            adj = 0.5, cex = rcexval*1.25, font = font)
      mtext(ctitle, side = 3, line = 2, at = ncolumns/2 + 1/2,
            adj = 0.5, cex = ccexval*1.25, font = font)
      if (!is.null(ctitle)) 
      { # Modification to implement clabelsreverse - 15/12/2012
        #else option is original code
        if (clabelsreverse)
          mtext(rev(collabs), side = 3, line = 0, at = seq(ncolumns),
                cex = ccexval)
        else
          mtext(collabs, side = 3, line = 0, at = seq(ncolumns),
                cex = ccexval)
      }
    }
    else 
    { 
      if (!is.null(rtitle)) 
      { 
        if (rlabelsreverse)
          mtext(rev(rowlabs), side = 2, line = 0, at = -seq(nrows),
                cex = rcexval, adj = 1, las = 1)
        else
          mtext(rowlabs, side = 2, line = 0, at = -seq(nrows),
                cex = rcexval, adj = 1, las = 1)
      }
      mtext(rtitle, side = 2, line = lineval, at = -nrows/2 - 1/2,
            adj = 0.5, cex = rcexval*1.25, font = font)
      mtext(ctitle, side = 3, line = 2, at = ncolumns/2 + 1/2,
            adj = 0.5, cex = ccexval*1.25, font = font)
      if (!is.null(ctitle)) 
      { 
        if (clabelsreverse)
          mtext(rev(collabs), side = 3, line = 0, at = seq(ncolumns),
                cex = ccexval)
        else
          mtext(collabs, side = 3, line = 0, at = seq(ncolumns),
                cex = ccexval)
      }
    }
  }
  for (i in labels) 
  {   
    x <- drange[dtrt == i]
    y <- drow[dtrt == i]
    for (j in seq(x)) 
    {   
      xo <- x[j] + c(0.5, 0.5, -0.5, -0.5, 0.5) * ccellpropn
      yo <- y[j] + c(-0.5, 0.5, 0.5, -0.5, -0.5) * rcellpropn
      if (plotcellboundary) 
      {  
        if (length(cellfillcolour) > 1)
          polygon(xo, yo, col=cellfillcolour[match(i, labels)], ...)
        else
          polygon(xo, yo, col=cellfillcolour, ...)
      }
      if (plotlabels) 
      {   
        if (!is.null(altlabels)) 
        { text(x, y, labels = altlabels[match(i, labels)], cex = cexval)#/0.7)
        }
        else 
        { text(x, y, labels = i, srt = charot, adj = 0.5,
               cex = cexval)#/0.7)
        }
      }
    }
  }

  blockboundaryPlot(blockdefinition = blockdefinition, blocksequence = blocksequence, 
                    rstart= 0, cstart = 0, nrows = nrows, ncolumns = ncolumns, 
                    blocklinecolour = blocklinecolour, blocklinewidth = blocklinewidth)
  invisible()
}
