pkgname <- "dae"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "dae-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('dae')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Zncsspline")
### * Zncsspline

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Zncsspline
### Title: Calculates the design matrix for fitting the random component of
###   a natural cubic smoothing spline
### Aliases: Zncsspline
### Keywords: array design

### ** Examples

Z <- Zncsspline(knot.points = 1:10, Gpower = 0.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Zncsspline", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as.data.frame.pstructure")
### * as.data.frame.pstructure

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as.data.frame.pstructure
### Title: Coerces a pstructure.object to a data.frame.
### Aliases: as.data.frame.pstructure
### Keywords: array projector

### ** Examples

## Generate a data.frame with 4 factors, each with three levels, in standard order
ABCD.lay <- fac.gen(list(A = 3, B = 3, C = 3, D = 3))

## create a pstructure object based on the formula ((A*B)/C)*D
ABCD.struct <- pstructure.formula(~ ((A*B)/C)*D, data =ABCD.lay)

## print the object either using the Method function or the generic function 
ABCS.dat <- as.data.frame.pstructure(ABCD.struct)
as.data.frame(ABCD.struct)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as.data.frame.pstructure", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as.numfac")
### * as.numfac

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as.numfac
### Title: Convert a factor to a numeric vector
### Aliases: as.numfac
### Keywords: factor manip

### ** Examples

## set up a factor and convert it to a numeric vector
a <- factor(rep(1:3, 4))
x <- as.numfac(a)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as.numfac", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blockboundaryPlot")
### * blockboundaryPlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blockboundaryPlot
### Title: This function plots a block boundary on a plot produced by
###   'designPlot'.
### Aliases: blockboundaryPlot
### Keywords: design plot

### ** Examples
## Not run: 
##D     SPL.Lines.mat <- matrix(as.numfac(Lines), ncol=16, byrow=T)
##D     colnames(SPL.Lines.mat) <- 1:16
##D     rownames(SPL.Lines.mat) <- 1:10
##D     SPL.Lines.mat <- SPL.Lines.mat[10:1, 1:16]
##D     designPlot(SPL.Lines.mat, labels=1:10, new=TRUE,
##D                rtitle="Rows",ctitle="Columns", 
##D                chardivisor=3, rcellpropn = 1, ccellpropn=1,
##D                plotcellboundary = TRUE)
##D     #Plot Mainplot boundaries
##D     blockboundaryPlot(blockdefinition = cbind(4,16), rstart = 1, 
##D                       blocklinewidth = 3, blockcolour = "green", 
##D                       nrows = 9, ncolumns = 16)
##D     blockboundaryPlot(blockdefinition = cbind(1,4), 
##D                       blocklinewidth = 3, blockcolour = "green", 
##D                       nrows = 1, ncolumns = 16)
##D     blockboundaryPlot(blockdefinition = cbind(1,4), rstart= 9, nrows = 10, ncolumns = 16, 
##D                       blocklinewidth = 3, blockcolour = "green")
##D     #Plot all 4 block boundaries            
##D     blockboundaryPlot(blockdefinition = cbind(8,5,5,4), blocksequence=T, 
##D                       cstart = 1, rstart= 1, nrows = 9, ncolumns = 15, 
##D                       blocklinewidth = 3,blockcolour = "blue")
##D     blockboundaryPlot(blockdefinition = cbind(10,16), blocklinewidth=3, blockcolour="blue", 
##D                       nrows=10, ncolumns=16)
##D     #Plot border and internal block boundaries only
##D     blockboundaryPlot(blockdefinition = cbind(8,14), cstart = 1, rstart= 1, 
##D                       nrows = 9, ncolumns =  15,
##D                       blocklinewidth = 3, blockcolour = "blue")
##D     blockboundaryPlot(blockdefinition = cbind(10,16), 
##D                       blocklinewidth = 3, blockcolour = "blue", 
##D                       nrows = 10, ncolumns = 16)
## End(Not run)


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blockboundaryPlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("correct.degfree")
### * correct.degfree

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: correct.degfree
### Title: Check the degrees of freedom in an object of class projector
### Aliases: correct.degfree
### Keywords: array projector

### ** Examples

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom
degfree(proj.m) <- 1
    
## check degrees of freedom are correct
correct.degfree(proj.m)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("correct.degfree", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("decomp.relate")
### * decomp.relate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: decomp.relate
### Title: Examines the relationship between the eigenvectors for two
###   decompositions
### Aliases: decomp.relate
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## obtain intra- and inter-block decompositions
decomp.inter <- proj2.eigen(unit.struct$Q[["Block"]], trt.struct$Q[["trt"]])
decomp.intra <- proj2.eigen(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])

## check that intra- and inter-block decompositions are orthogonal
decomp.relate(decomp.intra, decomp.inter) 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("decomp.relate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("degfree")
### * degfree

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: degfree
### Title: Degrees of freedom extraction and replacement
### Aliases: degfree degfree<-
### Keywords: array projector

### ** Examples

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## coerce to a projector
proj.m <- projector(m)

## extract its degrees of freedom
degfree(proj.m)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom and print the projector
degfree(proj.m) <- proj.m
print(proj.m)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("degfree", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designAmeasures")
### * designAmeasures

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designAmeasures
### Title: Calculates the average variance of pairwise differences from the
###   variance matrix for predictions that can be obtained using
###   mat.Vpredicts
### Aliases: designAmeasures
### Keywords: design

### ** Examples

## Reduced example from Smith et al. (2015) 
## Generate two-phase design
mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                            levels = c("Y","W","G","M","D","E"))
start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
rownames(start.design) <- NULL

## Set up matrices
n <- nrow(start.design)
W <- model.matrix(~ -1+ Variety, start.design)
ng <- ncol(W)
Gg<- diag(1, ng)
Vu <- with(start.design, fac.vcmat(Mrep, 0.3) + 
                         fac.vcmat(fac.combine(list(Mrep, Mday)), 0.2) + 
                         fac.vcmat(Frep, 0.1) + 
                         fac.vcmat(fac.combine(list(Frep, Fplot)), 0.2))
R <- diag(1, n)
  
## Calculate the varaince matrix of the predicted random Variety effects
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
  
## Calculate A-optimality measure
designAmeasures(Vp)
designAmeasures(Vp, groups=list(fldUndup = c(1:4), fldDup = c(5,6)))
grpsizes <- c(4,2)
names(grpsizes) <- c("fldUndup", "fldDup")
designAmeasures(Vp, groupsizes = grpsizes)
designAmeasures(Vp, groupsizes = c(4))
designAmeasures(Vp, groups=list(c(1,4),c(5,6)))

## Calculate the variance matrix of the predicted fixed Variety effects, elminating the grand mean
Vp.reduc <- mat.Vpred(W = W, Gg = 0, Vu = Vu, R = R, 
                      eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
## Calculate A-optimality measure
designAmeasures(Vp.reduc)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designAmeasures", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designAnatomy")
### * designAnatomy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designAnatomy
### Title: Given the layout for a design, obtain its anatomy via the
###   canonical analysis of its projectors to show the confounding and
###   aliasing inherent in the design.
### Aliases: designAnatomy
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain combined decomposition and summarize
unit.trt.canon <- designAnatomy(formulae = list(unit=~ Block/Unit, trt=~ trt),
                                data = PBIBD2.lay)
summary(unit.trt.canon, which.criteria = c("aeff","eeff","order"))
summary(unit.trt.canon, which.criteria = c("aeff","eeff","order"), labels.swap = TRUE)

## Three-phase sensory example from Brien and Payne (1999)
## Not run: 
##D data(Sensory3Phase.dat)
##D Eval.Field.Treat.canon <- designAnatomy(formulae = list(
##D                               eval= ~ ((Occasions/Intervals/Sittings)*Judges)/Positions, 
##D                               field= ~ (Rows*(Squares/Columns))/Halfplots,
##D                               treats= ~ Trellis*Method),
##D                                         data = Sensory3Phase.dat)
##D summary(Eval.Field.Treat.canon, which.criteria =c("aefficiency", "order"))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designAnatomy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designBlocksGGPlot")
### * designBlocksGGPlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designBlocksGGPlot
### Title: Adds block boundaries to a plot produced by 'designGGPlot'.
### Aliases: designBlocksGGPlot
### Keywords: design plot

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designBlocksGGPlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designGGPlot")
### * designGGPlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designGGPlot
### Title: Plots labels on two-way grids of coloured cells using 'ggplot2'
### Aliases: designGGPlot
### Keywords: aplot hplot design

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designGGPlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designLatinSqrSys")
### * designLatinSqrSys

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designLatinSqrSys
### Title: Generate a systematic plan for a Latin Square design
### Aliases: designLatinSqrSys
### Keywords: array design

### ** Examples

   matrix(designLatinSqrSys(5, start = c(seq(1, 5, 2), seq(2, 5, 2))), nrow=5)
   designLatinSqrSys(3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designLatinSqrSys", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designPlot")
### * designPlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designPlot
### Title: A graphical representation of an experimental design using
###   labels stored in a matrix.
### Aliases: designPlot
### Keywords: design plot

### ** Examples
## Not run: 
##D   designPlot(des.mat, labels=1:4, cellfillcolour="lightblue", new=TRUE, 
##D              plotcellboundary = TRUE, chardivisor=3, 
##D              rtitle="Lanes", ctitle="Positions", 
##D              rcellpropn = 1, ccellpropn=1)
##D   designPlot(des.mat, labels=5:87, plotlabels=TRUE, cellfillcolour="grey", new=FALSE,
##D              plotcellboundary = TRUE, chardivisor=3)
##D   designPlot(des.mat, labels=88:434, plotlabels=TRUE, cellfillcolour="lightgreen", 
##D              new=FALSE, plotcellboundary = TRUE, chardivisor=3, 
##D              blocksequence=TRUE, blockdefinition=cbind(4,10,12), 
##D              blocklinewidth=3, blockcolour="blue")
## End(Not run)


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designPlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designPlotlabels")
### * designPlotlabels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designPlotlabels
### Title: Plots labels on a two-way grid using 'ggplot2'
### Aliases: designPlotlabels
### Keywords: aplot hplot design

### ** Examples

Treatments <- factor(rep(1:6, times = 5))
RCBD.lay <- designRandomize(allocated = Treatments,
                            recipient = list(Blocks = 5, Units = 6),
                            nested.recipients = list(Units = "Blocks"),
                            seed = 74111)
designPlotlabels(RCBD.lay, labels = "Treatments", 
                 grid.x = "Units", grid.y = "Blocks",
                 colour.column = "Treatments", size = 5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designPlotlabels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designRandomize")
### * designRandomize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designRandomize
### Title: Randomize allocated to recipient factors to produce a layout for
###   an experiment
### Aliases: designRandomize
### Keywords: design factor datagen

### ** Examples

## Generate a randomized layout for a 4 x 4 Latin square
## (the nested.recipients argument is not needed here as none of the 
## factors are nested)
## Firstly, generate a systematic layout
LS.sys <- cbind(fac.gen(list(row = c("I","II","III","IV"), 
                             col = c(0,2,4,6))),
                treat = factor(designLatinSqrSys(4), label = LETTERS[1:4]))
## obtain randomized layout
LS.lay <- designRandomize(allocated = LS.sys["treat"], 
                          recipient = LS.sys[c("row","col")], 
                          seed = 7197132, unit.permutation = TRUE) 
LS.lay[LS.lay$.Permutation,]

## Generate a randomized layout for a replicated randomized complete 
## block design, with the block factors arranged in standard order for 
## rep then plot and then block
## Firstly, generate a systematic order such that levels of the 
## treatment factor coincide with plot
RCBD.sys <- cbind(fac.gen(list(rep = 2, plot=1:3, block = c("I","II"))),
                  tr = factor(rep(1:3, each=2, times=2)))
## obtain randomized layout
RCBD.lay <- designRandomize(allocated = RCBD.sys["tr"], 
                            recipient = RCBD.sys[c("rep", "block", "plot")], 
                            nested.recipients = list(plot = c("block","rep"), 
                                                     block="rep"), 
                            seed = 9719532, 
                            unit.permutation = TRUE)
#sort into the original standard order
RCBD.perm <- RCBD.lay[RCBD.lay$.Permutation,]
#resort into randomized order
RCBD.lay <- RCBD.perm[order(RCBD.perm$.Units),]

## Generate a layout for a split-unit experiment in which: 
## - the main-unit factor is A with 4 levels arranged in 
##   a randomized complete block design with 2 blocks;
## - the split-unit factor is B with 3 levels.
## Firstly, generate a systematic layout
SPL.sys <- cbind(fac.gen(list(block = 2, main.unit = 4, split.unit = 3)),
                 fac.gen(list(A = 4, B = 3), times = 2))
## obtain randomized layout
SPL.lay <- designRandomize(allocated = SPL.sys[c("A","B")], 
                           recipient = SPL.sys[c("block", "main.unit", "split.unit")], 
                           nested.recipients = list(main.unit = "block", 
                                                    split.unit = c("block", "main.unit")), 
                           seed=155251978)

## Generate a permutation of Seedlings within Species
seed.permute <- designRandomize(recipient = list(Species = 3, Seedlings = 4),
                                nested.recipients = list(Seedlings = "Species"),
                                seed = 75724, except = "Species", 
                                unit.permutation = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designRandomize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("designTwophaseAnatomies")
### * designTwophaseAnatomies

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: designTwophaseAnatomies
### Title: Given the layout for a design and three structure formulae,
###   obtain the anatomies for the two-phase, first-phase, cross-phase and
###   second-phase designs.
### Aliases: designTwophaseAnatomies
### Keywords: array design projector

### ** Examples

  #'## Microarray example from Jarrett & Ruggiero (2008) - see Brien (2019)
  jr.lay <- fac.gen(list(Set = 7, Dye = 2, Array = 3))
  jr.lay <- within(jr.lay, 
                   { 
                     Block <- factor(rep(1:7, each=6))
                     Plant <- factor(rep(c(1,2,3,2,3,1), times=7))
                     Sample <- factor(c(rep(c(2,1,2,2,1,1, 1,2,1,1,2,2), times=3), 
                                        2,1,2,2,1,1))
                     Treat <- factor(c(1,2,4,2,4,1, 2,3,5,3,5,2, 3,4,6,4,6,3, 
                                       4,5,7,5,7,4, 5,6,1,6,1,5, 6,7,2,7,2,6, 
                                       7,1,3,1,3,7),
                                     labels=c("A","B","C","D","E","F","G"))
                   })
  
  jr.anat <- designTwophaseAnatomies(formulae = list(array = ~ (Set:Array)*Dye,
                                                     plot = ~ Block/Plant/Sample,
                                                     trt = ~ Treat),
                                     which.designs = c("first","cross"), 
                                     data = jr.lay)  

## Three-phase sensory example from Brien and Payne (1999)
## Not run: 
##D data(Sensory3Phase.dat)
##D Sensory.canon <- designTwophaseAnatomies(formulae = list(
##D                               eval= ~ ((Occasions/Intervals/Sittings)*Judges)/Positions, 
##D                               field= ~ (Rows*(Squares/Columns))/Halfplots,
##D                               treats= ~ Trellis*Method),
##D                                         data = Sensory3Phase.dat)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("designTwophaseAnatomies", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("detect.diff")
### * detect.diff

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: detect.diff
### Title: Computes the detectable difference for an experiment
### Aliases: detect.diff
### Keywords: design

### ** Examples

## Compute the detectable difference for a randomized complete block design 
## with four treatments given power is 0.8 and alpha is 0.05. 
rm <- 5
detect.diff(rm = rm, df.num = 3, df.denom = 3 * (rm - 1),sigma = sqrt(20))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("detect.diff", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("efficiencies")
### * efficiencies

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: efficiencies
### Title: Extracts the canonical efficiency factors from a 'pcanon.object'
###   or a 'p2canon.object'.
### Aliases: efficiencies efficiencies.pcanon efficiencies.p2canon
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain combined decomposition using designAnatomy and get the efficiencies
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), data = PBIBD2.lay)
efficiencies.pcanon(unit.trt.canon)

##obtain the projectors for each formula using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition projs.2canon and get the efficiencies
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
efficiencies.p2canon(unit.trt.p2canon)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("efficiencies", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("efficiency.criteria")
### * efficiency.criteria

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: efficiency.criteria
### Title: Computes efficiency criteria from a set of efficiency factors
### Aliases: efficiency.criteria
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## save intrablock efficiencies
eff.inter <- proj2.efficiency(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])

## calculate efficiency criteria
efficiency.criteria(eff.inter)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("efficiency.criteria", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("elements")
### * elements

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: elements
### Title: Extract the elements of an array specified by the subscripts
### Aliases: elements
### Keywords: array manip

### ** Examples

## Form a table of the means for all combinations of Row and Line.
## Then obtain the values corresponding to the combinations in the data frame x,
## excluding Row 3.
x <- fac.gen(list(Row = 2, Line = 4), each =2)
x$y <- rnorm(16)
RowLine.tab <- tapply(x$y, list(x$Row, x$Line), mean)
xs <- elements(RowLine.tab, subscripts=x[x$"Line" != 3, c("Row", "Line")])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("elements", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extab")
### * extab

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extab
### Title: Expands the values in table to a vector
### Aliases: extab
### Keywords: manip

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extab", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.ar1mat")
### * fac.ar1mat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.ar1mat
### Title: forms the ar1 correlation matrix for a (generalized) factor
### Aliases: fac.ar1mat
### Keywords: array

### ** Examples

## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a 12 x 12 ar1 matrix corrresponding to B
ar1.B <- fac.ar1mat(B, 0.6)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.ar1mat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.combine")
### * fac.combine

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.combine
### Title: Combines several factors into one
### Aliases: fac.combine
### Keywords: factor manip

### ** Examples

## set up two factors
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## obtain six-level factor corresponding to the combinations of A and B
AB <- fac.combine(list(A,B))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.combine", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.divide")
### * fac.divide

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.divide
### Title: Divides a factor into several individual factors
### Aliases: fac.divide
### Keywords: factor manip

### ** Examples

## generate a small completely randomized design for 6 treatments 
n <- 12
CRD.unit <- list(Unit = n)
treat <- factor(rep(1:4, each = 3))
CRD.lay <- designRandomize(allocated = treat, recipient = CRD.unit, seed=956)

## divide the treatments into two two-level factor A nd B
CRD.facs <- fac.divide(CRD.lay$treat, factor.names = list(A = 2, B = 2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.divide", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.gen")
### * fac.gen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.gen
### Title: Generate all combinations of several factors and, optionally,
###   replicate them
### Aliases: fac.gen
### Keywords: design factor datagen

### ** Examples

## generate a 2^3 factorial experiment with levels - and +, and 
## in Yates order
mp <- c("-", "+")
fnames <- list(Catal = mp, Temp = mp, Press = mp, Conc = mp)
Fac4Proc.Treats <- fac.gen(generate = fnames, order="yates")

## Generate the factors A, B and D. The basic pattern has 4 repetitions
## of the levels of D for each A and B combination and 3 repetitions of 
## the pattern of the B and D combinations for each level of A. This basic 
## pattern has each combination repeated twice, and the whole of this 
## is repeated twice. It generates 864 A, B and D combinations.
gen <- list(A = 3, 3, B = c(0,100,200), 4, D = c("0","1"))
fac.gen(gen, times=2, each=2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.gen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.match")
### * fac.match

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.match
### Title: Match, for each combination of a set of columns in 'x', the row
###   that has the same combination in 'table'
### Aliases: fac.match
### Keywords: design factor

### ** Examples
## Not run: 
##D #A single unmatched combination
##D kdata <- data.frame(Expt="D197-5", 
##D                     Row=8, 
##D                     Column=20, stringsAsFactors=FALSE)
##D index <- fac.match(kdata, D197.dat, c("Expt", "Row", "Column"))
##D 
##D # A matched and an unmatched combination
##D kdata <- data.frame(Expt=c("D197-5", "D197-4"), 
##D                     Row=c(8, 10), 
##D                     Column=c(20, 8), stringsAsFactors=FALSE)
##D index <- fac.match(kdata, D197.dat, c("Expt", "Row", "Column"))
## End(Not run)


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.match", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.meanop")
### * fac.meanop

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.meanop
### Title: computes the projection matrix that produces means
### Aliases: fac.meanop
### Keywords: array projector

### ** Examples

## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a generalized factor whose levels are the combinations of A and B
AB <- fac.combine(list(A,B))

## obtain the operator that computes the AB means from a vector of length 12
M.AB <- fac.meanop(AB)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.meanop", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.multinested")
### * fac.multinested

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.multinested
### Title: Creates several factors, one for each level of nesting.fac and
###   each of whose values are either generated within those of a level of
###   nesting.fac or using the values of nested.fac within a levels of
###   nesting.fac.
### Aliases: fac.multinested
### Keywords: factor manip

### ** Examples

  lay <- data.frame(A = factor(rep(c(1:3), c(3,6,4)), labels = letters[1:3]))
  lay$B <-fac.nested(lay$A)

  #Add factors for B within each level of A
  lay2 <- cbind(lay, fac.multinested(lay$A))
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summary(canon2)

  #Add factors for B within each level of A, but with levels and outlabel given
  lay2 <- cbind(lay, fac.multinested(lay$A, nested.levs = seq(10,60,10), outlabel = "other"))
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summary(canon2)

  #Replicate the combinations of A and B three times and index them with the factor sample
  lay3 <- rbind(lay,lay,lay)
  lay3$sample <- with(lay3, fac.nested(fac.combine(list(A,B))))
  
  #Add factors for B within each level of A
  lay4 <- cbind(lay3, fac.multinested(nesting.fac = lay$A, nested.fac = lay$B))
  canon4 <- designAnatomy(list(~(A/(a+b+c))/sample), data = lay4)
  summary(canon4)

  #Add factors for sample within each combination of A and B
  lay5 <- with(lay4, cbind(lay4, 
                           fac.multinested(nesting.fac = a, fac.prefix = "a"),
                           fac.multinested(nesting.fac = b, fac.prefix = "b"),
                           fac.multinested(nesting.fac = c, fac.prefix = "c")))

  #Add factors for sample within each level of A
  lay6 <- cbind(lay4, 
                fac.multinested(nesting.fac = lay4$A, nested.fac = lay$sample, fac.prefix = "samp"))
  canon6 <- designAnatomy(list(~A/(a/sampa+b/sampb+c/sampc)), data = lay6)
  summary(canon6)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.multinested", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.nested")
### * fac.nested

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.nested
### Title: creates a factor, the nested factor, whose values are generated
###   within those of the factor nesting.fac
### Aliases: fac.nested
### Keywords: factor manip

### ** Examples

## set up factor A
A <- factor(c(1, 1, 1, 2, 2))

## create nested factor
B <- fac.nested(A)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.nested", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.recode")
### * fac.recode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.recode
### Title: Recodes factor 'levels' using values in a vector. The values in
###   the vector do not have to be unique.
### Aliases: fac.recode
### Keywords: factor manip

### ** Examples

## set up a factor with labels
a <- factor(rep(1:4, 4), labels=c("A","B","C","D"))
 
## recode "A" and "D" to 1 and "B" and "C" to 2
b <- fac.recode(a, c(1,2,2,1), labels = c("a","b"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.recode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.sumop")
### * fac.sumop

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.sumop
### Title: computes the summation matrix that produces sums corresponding
###   to a (generalized) factor
### Aliases: fac.sumop
### Keywords: array projector

### ** Examples

## set up a two-level factoir and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a generlaized factor whose levels are the combinations of A and B
AB <- fac.combine(list(A,B))

## obtain the operator that computes the AB means from a vector of length 12
S.AB <- fac.sumop(AB)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.sumop", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac.vcmat")
### * fac.vcmat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac.vcmat
### Title: forms the variance matrix for the variance component of a
###   (generalized) factor
### Aliases: fac.vcmat
### Keywords: array

### ** Examples

## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## create a 12 x 12 ar1 matrix corrresponding to B
vc.B <- fac.vcmat(B, 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac.vcmat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitted.aovlist")
### * fitted.aovlist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitted.aovlist
### Title: Extract the fitted values for a fitted model from an aovlist
###   object
### Aliases: fitted.aovlist fitted
### Keywords: methods models htest

### ** Examples

## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## two equivalent ways of extracting the fitted values
fit  <- fitted.aovlist(RCBDPen.aov)
fit <- fitted(RCBDPen.aov, error.term = "Blend:Flask")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitted.aovlist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitted.errors")
### * fitted.errors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitted.errors
### Title: Extract the fitted values for a fitted model
### Aliases: fitted.errors
### Keywords: models htest

### ** Examples

## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## three equivalent ways of extracting the fitted values
fit  <- fitted.aovlist(RCBDPen.aov)
fit <- fitted(RCBDPen.aov, error.term = "Blend:Flask")
fit <- fitted.errors(RCBDPen.aov, error.term = "Blend:Flask")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitted.errors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get.daeTolerance")
### * get.daeTolerance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get.daeTolerance
### Title: Gets the value of daeTolerance for the package dae
### Aliases: get.daeTolerance
### Keywords: manip projector

### ** Examples

## get daeTolerance.
get.daeTolerance()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get.daeTolerance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("harmonic.mean")
### * harmonic.mean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: harmonic.mean
### Title: Calcuates the harmonic mean.
### Aliases: harmonic.mean
### Keywords: manip

### ** Examples

y <- c(seq(0.1,1,0.2))
harmonic.mean(y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("harmonic.mean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("interaction.ABC.plot")
### * interaction.ABC.plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: interaction.ABC.plot
### Title: Plots an interaction plot for three factors
### Aliases: interaction.ABC.plot
### Keywords: aplot hplot design

### ** Examples

## Not run: 
##D ## plot for Example 14.1 from Mead, R. (1990). The Design of Experiments: 
##D ## Statistical Principles for Practical Application. Cambridge, 
##D ## Cambridge University Press.  
##D ## use ?SPLGrass.dat for details
##D data(SPLGrass.dat)
##D interaction.ABC.plot(Main.Grass, x.factor=Period,
##D                      groups.factor=Spring, trace.factor=Summer,
##D                      data=SPLGrass.dat,
##D                      title="Effect of Period, Spring and Summer on Main Grass")
##D 
##D ## plot for generated data
##D ## use ?ABC.Interact.dat for data set details
##D data(ABC.Interact.dat)
##D ## Add standard errors for plotting 
##D ## - here data contains a single value for each combintion of A, B and C
##D ## - need to supply name for data twice 
##D ABC.Interact.dat$se <- rep(c(0.5,1), each=4)
##D interaction.ABC.plot(MOE, A, B, C, data=ABC.Interact.dat,
##D                      ggplotFunc=list(geom_errorbar(data=ABC.Interact.dat, 
##D                                                    aes(ymax=MOE+se, ymin=MOE-se), 
##D                                                    width=0.2)))
## End(Not run)


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("interaction.ABC.plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is.allzero")
### * is.allzero

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is.allzero
### Title: Tests whether all elements are approximately zero
### Aliases: is.allzero
### Keywords: manip

### ** Examples

## create a vector of 9 zeroes and a one
y <- c(rep(0,9), 1)

## check that vector is only zeroes is FALSE 
is.allzero(y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is.allzero", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is.projector")
### * is.projector

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is.projector
### Title: Tests whether an object is a valid object of class projector
### Aliases: is.projector
### Keywords: array projector

### ** Examples

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## check that it is a valid projector
is.projector(proj.m)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is.projector", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("marginality")
### * marginality

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: marginality
### Title: Extracts the marginality matrix (matrices) from a
###   'pstructure.object' or a 'pcanon.object'.
### Aliases: marginality.pcanon marginality.pstructure marginality
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain pstructure.object and extract marginality matrix
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
unit.marg <- marginality(unit.struct)

##obtain combined decomposition and extract marginality matrices
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), data = PBIBD2.lay)
marg <- marginality(unit.trt.canon)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("marginality", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.I")
### * mat.I

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.I
### Title: Forms a unit matrix
### Aliases: mat.I
### Keywords: array

### ** Examples

    col.I <- mat.I(order=4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.I", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.J")
### * mat.J

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.J
### Title: Forms a square matrix of ones
### Aliases: mat.J
### Keywords: array

### ** Examples

    col.J <- mat.J(order=4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.J", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.Vpred")
### * mat.Vpred

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.Vpred
### Title: Calculates the variances of a set of predicted effects from a
###   mixed model
### Aliases: mat.Vpred
### Keywords: array design

### ** Examples

## Reduced example from Smith et al. (2015)
## Generate two-phase design
mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                            levels = c("Y","W","G","M","D","E"))
start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
rownames(start.design) <- NULL

## Set up matrices
n <- nrow(start.design)
W <- model.matrix(~ -1+ Variety, start.design)
ng <- ncol(W)
Gg<- diag(1, ng)
Vu <- with(start.design, fac.vcmat(Mrep, 0.3) + 
                         fac.vcmat(fac.combine(list(Mrep, Mday)), 0.2) + 
                         fac.vcmat(Frep, 0.1) + 
                         fac.vcmat(fac.combine(list(Frep, Fplot)), 0.2))
R <- diag(1, n)
  
## Calculate the variance matrix of the predicted random Variety effects
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
designAmeasures(Vp)

## Calculate the variance matrix of the predicted fixed Variety effects, 
## elminating the grand mean
Vp.reduc <- mat.Vpred(W = W, Gg = 0, Vu = Vu, R = R, 
                      eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
designAmeasures(Vp.reduc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.Vpred", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.Vpredicts")
### * mat.Vpredicts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.Vpredicts
### Title: Calculates the variances of a set of predicted effects from a
###   mixed model, based on supplied matrices or formulae.
### Aliases: mat.Vpredicts
### Keywords: array design

### ** Examples

## Reduced example from Smith et al. (2015)
## Generate two-phase design
mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                            levels = c("Y","W","G","M","D","E"))
start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
rownames(start.design) <- NULL

## Set gammas
terms <- c("Variety", "Frep", "Frep:Fplot", "Mrep", "Mrep:Mday", "Mrep:Mday:Mord")
gammas <- c(1, 0.1, 0.2, 0.3, 0.2, 1)
names(gammas) <- terms

## Specify matrices to calculate the variance matrix of the predicted fixed Variety effects 
W <- model.matrix(~ -1 + Variety, start.design)
Vu <- with(start.design, fac.vcmat(Mrep, gammas["Mrep"]) + 
                         fac.vcmat(fac.combine(list(Mrep,Mday)), gammas["Mrep:Mday"]) + 
                         fac.vcmat(Frep, gammas["Frep"]) + 
                         fac.vcmat(fac.combine(list(Frep,Fplot)), gammas["Frep:Fplot"]))
R <- diag(1, nrow(start.design))
  
## Calculate variance matrix
Vp <- mat.Vpredicts(target = W, random=Vu, R=R)

## Calculate the variance matrix of the predicted random Variety effects using formulae
Vp <- mat.Vpredicts(target = ~ -1 + Variety, Gt = 1, 
                    fixed = ~ 1, 
                    random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                    G = as.list(gammas[c(4,5,2,3)]), 
                    R = R, design = start.design)
designAmeasures(Vp)

## Calculate the variance matrix of the predicted fixed Variety effects, 
## elminating the grand mean
n <- nrow(start.design)
Vp.reduc <- mat.Vpredicts(target = ~ -1 + Variety, 
                          random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                          G = as.list(gammas[c(4,5,2,3)]), 
                          eliminate = projector(matrix(1, nrow = n, ncol = n)/n), 
                          design = start.design)
designAmeasures(Vp.reduc)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.Vpredicts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ar1")
### * mat.ar1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ar1
### Title: Forms an ar1 correlation matrix
### Aliases: mat.ar1
### Keywords: array

### ** Examples

    corr <- mat.ar1(rho=0.4, order=4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ar1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ar2")
### * mat.ar2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ar2
### Title: Forms an ar2 correlation matrix
### Aliases: mat.ar2
### Keywords: array

### ** Examples

    corr <- mat.ar2(ARparameters = c(0.4, 0.2), order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ar2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ar3")
### * mat.ar3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ar3
### Title: Forms an ar3 correlation matrix
### Aliases: mat.ar3
### Keywords: array

### ** Examples

    corr <- mat.ar3(ARparameters = c(0.4, 0.2, 0.1), order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ar3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.arma")
### * mat.arma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.arma
### Title: Forms an arma correlation matrix
### Aliases: mat.arma
### Keywords: array

### ** Examples

    corr <- mat.arma(ARparameter = 0.4, MAparameter = -0.2, order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.arma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.banded")
### * mat.banded

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.banded
### Title: Form a banded matrix from a vector of values
### Aliases: mat.banded
### Keywords: array

### ** Examples

      m <- mat.banded(c(1,0.6,0.5), 5,5)
      m <- mat.banded(c(1,0.6,0.5), 3,4)
      m <- mat.banded(c(1,0.6,0.5), 4,3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.banded", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.dirprod")
### * mat.dirprod

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.dirprod
### Title: Forms the direct product of two matrices
### Aliases: mat.dirprod
### Keywords: array

### ** Examples

    col.I <- mat.I(order=4)
    row.I <- mat.I(order=28)
    V <- mat.dirprod(col.I, row.I)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.dirprod", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.dirsum")
### * mat.dirsum

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.dirsum
### Title: Forms the direct sum of a list of matrices
### Aliases: mat.dirsum
### Keywords: array

### ** Examples

       m1 <- matrix(1:4, nrow=2)
       m2 <- matrix(11:16, nrow=3)
       m3 <- diag(1, nrow=2, ncol=2)
       dsum <- mat.dirsum(list(m1, m2, m3))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.dirsum", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.exp")
### * mat.exp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.exp
### Title: Forms an exponential correlation matrix
### Aliases: mat.exp
### Keywords: array

### ** Examples

    corr <- mat.exp(coordinates=c(3:6, 9:12, 15:18), rho=0.1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.exp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.gau")
### * mat.gau

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.gau
### Title: Forms an exponential correlation matrix
### Aliases: mat.gau
### Keywords: array

### ** Examples

    corr <- mat.gau(coordinates=c(3:6, 9:12, 15:18), rho=0.1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.gau", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ma1")
### * mat.ma1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ma1
### Title: Forms an ma1 correlation matrix
### Aliases: mat.ma1
### Keywords: array

### ** Examples

    corr <- mat.ma1(MAparameter=0.4, order=4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ma1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ma2")
### * mat.ma2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ma2
### Title: Forms an ma2 correlation matrix
### Aliases: mat.ma2
### Keywords: array

### ** Examples

    corr <- mat.ma2(MAparameters = c(0.4, -0.2), order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ma2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.ncssvar")
### * mat.ncssvar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.ncssvar
### Title: Calculates the variance matrix of the random effects for a
###   natural cubic smoothing spline
### Aliases: mat.ncssvar
### Keywords: array design

### ** Examples

Gs <- mat.ncssvar(knot.points = 1:10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.ncssvar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.sar")
### * mat.sar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.sar
### Title: Forms an sar correlation matrix
### Aliases: mat.sar
### Keywords: array

### ** Examples

    corr <- mat.sar(SARparameter = -0.4, order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.sar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mat.sar2")
### * mat.sar2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.sar2
### Title: Forms an sar2 correlation matrix
### Aliases: mat.sar2
### Keywords: array

### ** Examples

    corr <- mat.sar2(gamma = c(-0.4, 0.2), order = 4)
    corr <- mat.sar2(gamma = c(-0.4, 0.2), order = 4, print = "ar3")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.sar2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mpone")
### * mpone

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mpone
### Title: Converts the first two levels of a factor into the numeric
###   values -1 and +1
### Aliases: mpone
### Keywords: factor manip

### ** Examples

## generate all combinations of two two-level factors
mp <- c("-", "+")
Frf3.trt <- fac.gen(list(A = mp, B = mp))

## add factor C, whose levels are the products of the levles of A and B
Frf3.trt$C <- factor(mpone(Frf3.trt$A)*mpone(Frf3.trt$B), labels = mp)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mpone", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("no.reps")
### * no.reps

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: no.reps
### Title: Computes the number of replicates for an experiment
### Aliases: no.reps
### Keywords: design

### ** Examples

## Compute the number of replicates (blocks) required for a randomized 
## complete block design with four treatments. 
no.reps(multiple = 1, df.num = 3,
        df.denom = expression(df.num * (r - 1)), delta = 5,
	        sigma = sqrt(20), print = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("no.reps", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("power.exp")
### * power.exp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: power.exp
### Title: Computes the power for an experiment
### Aliases: power.exp
### Keywords: design

### ** Examples

## Compute power for a randomized complete block design with four treatments 
## and five blocks. 
rm <- 5
power.exp(rm = rm, df.num = 3, df.denom = 3 * (rm - 1), delta = 5,
          sigma = sqrt(20),print = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("power.exp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.aliasing")
### * print.aliasing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.aliasing
### Title: Print an aliasing data.frame
### Aliases: print.aliasing
### Keywords: array projector

### ** Examples

## Generate a data.frame with 3 factors length 12
pseudo.lay <- data.frame(pl = factor(1:12),
                         ab = factor(rep(1:4, times=3)),
                         a = factor(rep(1:2, times=6)))


## create a pstructure object
trt.struct <- pstructure(~ ab+a, data = pseudo.lay)

## print the object either using the Method function, the generic function or show
print.aliasing(trt.struct$aliasing)
print(trt.struct$aliasing, which.criteria = "none")
trt.struct$aliasing



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.aliasing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.projector")
### * print.projector

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.projector
### Title: Print projectors
### Aliases: print.projector print,projector-method
### Keywords: array projector

### ** Examples

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## print the object either using the Method function, the generic function or show
print.projector(proj.m)
print(proj.m)
proj.m



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.projector", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.pstructure")
### * print.pstructure

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.pstructure
### Title: Prints a pstructure.object
### Aliases: print.pstructure
### Keywords: array projector

### ** Examples

## Generate a data.frame with 4 factors, each with three levels, in standard order
ABCD.lay <- fac.gen(list(A = 3, B = 3, C = 3, D = 3))

## create a pstructure object based on the formula ((A*B)/C)*D
ABCD.struct <- pstructure.formula(~ ((A*B)/C)*D, data =ABCD.lay)

## print the object either using the Method function, the generic function or show
print.pstructure(ABCD.struct)
print(ABCD.struct)
ABCD.struct



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.pstructure", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.summary.p2canon")
### * print.summary.p2canon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.summary.p2canon
### Title: Prints the values in an 'summary.p2canon' object
### Aliases: print.summary.p2canon
### Keywords: design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain projectors using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and print summary
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summ <- summary(unit.trt.p2canon)
print(summ)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.summary.p2canon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.summary.pcanon")
### * print.summary.pcanon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.summary.pcanon
### Title: Prints the values in an 'summary.pcanon' object
### Aliases: print.summary.pcanon
### Keywords: design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain combined decomposition and summarize
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt),
                                data = PBIBD2.lay)
summ <- summary(unit.trt.canon, which = c("aeff","eeff","order"))
print(summ)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.summary.pcanon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("proj2.combine")
### * proj2.combine

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: proj2.combine
### Title: Compute the projection and Residual operators for two, possibly
###   nonorthogonal, projectors
### Aliases: proj2.combine
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## obtain the projection operators for the interblock analysis
PBIBD2.Bops <- proj2.combine(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])
Q.B.T <- PBIBD2.Bops$Qconf
Q.B.res <- PBIBD2.Bops$Qres

## demonstrate their orthogonality
is.allzero(Q.B.T %*% Q.B.res)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("proj2.combine", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("proj2.efficiency")
### * proj2.efficiency

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: proj2.efficiency
### Title: Computes the canonical efficiency factors for the joint
###   decomposition of two projectors
### Aliases: proj2.efficiency
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## save intrablock efficiencies
eff.intra <- proj2.efficiency(unit.struct$Q[["Block"]], trt.struct$Q[["trt"]])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("proj2.efficiency", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("proj2.eigen")
### * proj2.eigen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: proj2.eigen
### Title: Canonical efficiency factors and eigenvectors in joint
###   decomposition of two projectors
### Aliases: proj2.eigen
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

## obtain intra- and inter-block decompositions
decomp.inter <- proj2.eigen(unit.struct$Q[["Block"]], trt.struct$Q[["trt"]])
decomp.intra <- proj2.eigen(unit.struct$Q[["Unit[Block]"]], trt.struct$Q[["trt"]])

#extract intrablock efficiencies
decomp.intra$efficiencies



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("proj2.eigen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("projector-class")
### * projector-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projector-class
### Title: Class projector
### Aliases: projector-class coerce,projector,matrix-method
###   coerce<-,projector,matrix-method
### Keywords: classes array projector

### ** Examples

showClass("projector")

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## check that it is a valid projector
is.projector(proj.m)

## create a projector based on the matrix m
proj.m <- new("projector", data=m)

## add its degrees of freedom and print the projector
degfree(proj.m) <- proj.m



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projector-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("projector")
### * projector

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projector
### Title: Create projectors
### Aliases: projector
### Keywords: array projector

### ** Examples

## set up a 2 x 2 mean operator that takes the mean of a vector of 2 values
m <- matrix(rep(0.5,4), nrow=2)

## create an object of class projector
proj.m <- projector(m)

## check that it is a valid projector
is.projector(proj.m)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projector", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("projs.2canon")
### * projs.2canon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projs.2canon
### Title: A canonical analysis of the relationships between two sets of
###   projectors
### Aliases: projs.2canon
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain projectors using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summary(unit.trt.p2canon)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projs.2canon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("projs.combine.p2canon")
### * projs.combine.p2canon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projs.combine.p2canon
### Title: Extract, from a p2canon object, the projectors that give the
###   combined canonical decomposition
### Aliases: projs.combine.p2canon
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

## obtain sets of projectors
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
UcombineT <- projs.combine.p2canon(unit.trt.p2canon)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projs.combine.p2canon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pstructure.formula")
### * pstructure.formula

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pstructure.formula
### Title: Takes a formula and constructs a 'pstructure.object' that
###   includes the orthogonalized projectors for the terms in a formula
### Aliases: pstructure.formula pstructure
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)
## manually obtain projectors for units
Q.G <- projector(matrix(1, nrow=24, ncol=24)/24)                         
Q.B <- projector(fac.meanop(PBIBD2.lay$Block) - Q.G)
Q.BP <- projector(diag(1, nrow=24) - Q.B - Q.G)

## manually obtain projector for trt
Q.T <- projector(fac.meanop(PBIBD2.lay$trt) - Q.G)

##compute intrablock efficiency criteria
effic <- proj2.efficiency(Q.BP, Q.T)
effic
efficiency.criteria(effic)

##obtain projectors using pstructure.formula
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summary(unit.trt.p2canon, which = c("aeff","eeff","order"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pstructure.formula", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("qqyeffects")
### * qqyeffects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: qqyeffects
### Title: Half or full normal plot of Yates effects
### Aliases: qqyeffects
### Keywords: iplot hplot design htest

### ** Examples

## analysis of 2^4 factorial experiment from Table 10.6 of Box, Hunter and 
## Hunter (1978) Statistics for Experimenters. New York, Wiley.
## use ?Fac4Proc.dat for data set details
data(Fac4Proc.dat)
Fac4Proc.aov <- aov(Conv ~ Catal * Temp * Press * Conc + Error(Runs),
                                                            Fac4Proc.dat)
qqyeffects(Fac4Proc.aov, error.term="Runs", data=Fac4Proc.dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("qqyeffects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("resid.errors")
### * resid.errors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: resid.errors
### Title: Extract the residuals for a fitted model
### Aliases: resid.errors
### Keywords: models htest

### ** Examples

## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## two equivalent ways of extracting the residuals
res  <- residuals.aovlist(RCBDPen.aov)
res <- residuals(RCBDPen.aov, error.term = "Blend:Flask")
res <- resid.errors(RCBDPen.aov)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("resid.errors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("residuals.aovlist")
### * residuals.aovlist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: residuals.aovlist
### Title: Extract the residuals from an aovlist object
### Aliases: residuals.aovlist residuals
### Keywords: methods models htest

### ** Examples

## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## two equivalent ways of extracting the residuals
res  <- residuals.aovlist(RCBDPen.aov)
res <- residuals(RCBDPen.aov, error.term = "Blend:Flask")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("residuals.aovlist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rmvnorm")
### * rmvnorm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rmvnorm
### Title: generates a vector of random values from a multivariate normal
###   distribution
### Aliases: rmvnorm
### Keywords: datagen

### ** Examples

## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## generate random values from a multivariate normal for which 
#the mean is 20 for all variables and 
#the variance matrix has random effects for factor A, ar1 pattern for B and 
#residual random variation
mean <- rep(20, 12)
V <- fac.vcmat(A, 5) + fac.ar1mat(B, 0.6) + 2*mat.I(12)
y <- rmvnorm(mean, V)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rmvnorm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("set.daeTolerance")
### * set.daeTolerance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: set.daeTolerance
### Title: Sets the values of daeTolerance for the package dae
### Aliases: set.daeTolerance
### Keywords: manip projector

### ** Examples

## set daeTolerance.
set.daeTolerance(1E-04, 1E-08)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("set.daeTolerance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("strength")
### * strength

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: strength
### Title: Generate paper strength values
### Aliases: strength
### Keywords: datagen design

### ** Examples

## Here temperature is a factor with 4*3 = 12 values whose
## first 3 values specify the temperatures to be applied in
## the 3 runs on the first day, values 4 to 6 specify the
## temperatures for the 3 runs on day 2, and so on.
temperature <- factor(rep(c(80,85,90), 4))
exp.strength <- strength(nodays = 4, noruns = 3,
                         temperature = temperature, ident = 0123456)

## In this second example, a completely randomized design is generated 
## for the same 3 temperatures replicated 4 times. The layout is stored 
## in the data.frame called Design.
Design <- designRandomize(allocated = temperature, 
                          recipient = list(runs = 12), 
                          seed = 5847123)
## eradicate the unrandomized version of temperature
remove("temperature")

## The 12 temperatures in Design are to be regarded as being assigned to 
## days and runs in the same manner as for the first example.
exp.strength <- strength(nodays = 4, noruns = 3,
                         temperature = Design$temperature, ident = 0123456)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("strength", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.p2canon")
### * summary.p2canon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.p2canon
### Title: Summarize a canonical analysis of the relationships between two
###   sets of projectors
### Aliases: summary.p2canon summary,p2canon-method
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain projectors using pstructure
unit.struct <- pstructure(~ Block/Unit, data = PBIBD2.lay)
trt.struct <- pstructure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(unit.struct$Q, trt.struct$Q)
summary(unit.trt.p2canon)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.p2canon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.pcanon")
### * summary.pcanon

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.pcanon
### Title: Summarizes the anatomy of a design, being the decomposition of
###   the sample space based on its canonical analysis, as produced by
###   designAnatomy
### Aliases: summary.pcanon summary,pcanon-method
### Keywords: array design projector

### ** Examples

## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- designRandomize(allocated = trt, 
                              recipient = PBIBD2.unit, 
                              nested.recipients = PBIBD2.nest)

##obtain combined decomposition and summarize
unit.trt.canon <- designAnatomy(list(unit=~ Block/Unit, trt=~ trt), 
                                data = PBIBD2.lay)
summary(unit.trt.canon, which = c("aeff","eeff","order"))
summary(unit.trt.canon, which = c("aeff","eeff","order"), labels.swap = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.pcanon", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tukey.1df")
### * tukey.1df

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tukey.1df
### Title: Performs Tukey's one-degree-of-freedom-test-for-nonadditivity
### Aliases: tukey.1df
### Keywords: models htest

### ** Examples

## set up data frame for randomized complete block design in Table 4.4 from 
## Box, Hunter and Hunter (2005) Statistics for Experimenters. 2nd edn 
## New York, Wiley.
RCBDPen.dat <- fac.gen(list(Blend=5, Flask=4))
RCBDPen.dat$Treat <- factor(rep(c("A","B","C","D"), times=5))
RCBDPen.dat$Yield <- c(89,88,97,94,84,77,92,79,81,87,87,
                       85,87,92,89,84,79,81,80,88)

## perform the analysis of variance
RCBDPen.aov <- aov(Yield ~ Blend + Treat + Error(Blend/Flask), RCBDPen.dat)
summary(RCBDPen.aov)

## Obtain the quantities for Tukey's test
tukey.1df(RCBDPen.aov, RCBDPen.dat, error.term = "Blend:Flask")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tukey.1df", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("yates.effects")
### * yates.effects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: yates.effects
### Title: Extract Yates effects
### Aliases: yates.effects
### Keywords: design htest

### ** Examples

## analysis of 2^4 factorial experiment from Table 10.6 of Box, Hunter and 
## Hunter (1978) Statistics for Experimenters. New York, Wiley.
## use ?Fac4Proc.dat for data set details
data(Fac4Proc.dat)
Fac4Proc.aov <- aov(Conv ~ Catal * Temp * Press * Conc + Error(Runs),
                                                            Fac4Proc.dat)
round(yates.effects(Fac4Proc.aov, error.term="Runs", data=Fac4Proc.dat), 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("yates.effects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
