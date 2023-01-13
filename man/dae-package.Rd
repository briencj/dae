\newcommand{\packageVer}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Version")}}
\newcommand{\packageDate}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1", fields="Date")}}
\name{dae-package}
\alias{dae-package}
\alias{dae}
\docType{package}
\title{
\packageTitle{dae}
}
\description{\packageDescription{dae}

\bold{Version:} {\packageVer{dae}}

\bold{Date:} \packageDate{dae}
}
\section{Index}{
\tabular{ll}{
(i) Data\cr\cr

\code{\link{ABC.Interact.dat}} 
  \tab Randomly generated set of values indexed by \cr
  \tab three factors \cr
\code{\link{BIBDWheat.dat}} 
  \tab Data for a balanced incomplete block experiment \cr
\code{\link{Casuarina.dat}}
  \tab Data for an experiment with rows and columns from \cr
  \tab Williams (2002) \cr
\code{\link{Exp249.munit.des}}
 \tab Systematic, main-plot design for an experiment to be \cr
 \tab run in a greenhouse \cr
\code{\link{Fac4Proc.dat}}
  \tab Data for a 2^4 factorial experiment\cr
\code{\link{LatticeSquare_t49.des}}
 \tab A Lattice square design for 49 treatments \cr
\code{\link{McIntyreTMV.dat}}
 \tab The design and data from McIntyre (1955) two-phase \cr
 \tab experiment \cr
\code{\link{Oats.dat}}
 \tab Data for an experiment to investigate nitrogen response \cr
 \tab of 3 oats varieties    \cr 
\code{\link{Sensory3Phase.dat}}
 \tab Data for the three-phahse sensory evaluation \cr
 \tab experiment in Brien, C.J. and Payne, R.W. (1999)\cr
\code{\link{Sensory3PhaseShort.dat}}
 \tab Data for the three-phase sensory evaluation \cr
 \tab experiment in Brien, C.J. and Payne, R.W. (1999), \cr
 \tab but with short factor names \cr
\code{\link{SPLGrass.dat}}
  \tab Data for an experiment to investigate the \cr
  \tab effects of grazing patterns on pasture \cr
  \tab composition \cr
\cr
(ii) Factor manipulation functions \cr\cr

\code{\link{as.numfac}}             
 \tab Convert a factor to a numeric vector \cr
\code{\link{fac.combine}}
  \tab Combines several factors into one \cr
\code{\link{fac.divide}}
  \tab Divides a factor into several separate factors\cr
\code{\link{fac.gen}}
  \tab Generate all combinations of several factors and, \cr
  \tab optionally, replicate them\cr
\code{\link{fac.gen}}
  \tab Generate all combinations of the levels of the supplied \cr
  \tab factors, without replication \cr
\code{\link{fac.match}}
  \tab Match, for each combination of a set of columns \cr
  \tab in 'x', the row that has the same combination\cr
  \tab in 'table'\cr
\code{\link{fac.multinested}}
  \tab Creates several factors, one for each level of a nesting.fac \cr
  \tab and each of whose values are either generated within \cr 
  \tab those of the level of nesting.fac or using the values \cr
  \tab of a nested.fac\cr
\code{\link{fac.nested}}
  \tab Creates a factor, the nested factor, whose values are \cr
  \tab generated within those of a nesting factor\cr
\code{\link{fac.recast}}
  \tab Recasts a factor by modifying the values in the factor vector\cr 
  \tab and/or the levels attribute, possibly combining \cr 
  \tab some levels into a single level. \cr
\code{\link{fac.recode}}
  \tab Recodes factor 'levels' using possibly nonunique \cr
  \tab values in a vector. (May be deprecated in future.)\cr
\code{\link{fac.split}}
  \tab Splits a factor whose levels consist of several delimited \cr
  \tab strings into several factors.   \cr
\code{\link{fac.uncombine}}
  \tab Cleaves a single factor, each of whose levels has delimited \cr 
  \tab strings, into several factors using the separated strings.  \cr
\code{\link{fac.uselogical}}
  \tab Forms a two-level factor from a logical object \cr
\code{\link{mpone}}
  \tab Converts the first two levels of a factor into\cr
  \tab the numeric values -1 and +1\cr
\cr
(iii) Design functions\cr\cr

\code{\link{blockboundaryPlot}}
  \tab This function plots a block boundary on a plot \cr
  \tab produced by 'designPlot'. It supersedes  \cr
  \tab blockboundary.plot. \cr
\code{\link{designAmeasures}}
  \tab Calculates the A-optimality measures from the \cr
  \tab variance matrix for predictions.\cr
\code{\link{designAnatomy}}
  \tab Given the layout for a design, obtain its anatomy via \cr 
  \tab the canonical analysis of its projectors to show the \cr
  \tab confounding and aliasing inherent in the design.\cr
\code{\link{designBlocksGGPlot}}
  \tab Adds block boundaries to a plot produced by \code{\link{designGGPlot}}. \cr
\code{\link{designGGPlot}}
  \tab Plots labels on a two-way grid of coloured cells using \code{ggplot2} \cr
  \tab to represent an experimental design. \cr
\code{\link{designLatinSqrSys}}
  \tab Generate a systematic plan for a Latin Square design.\cr
\code{\link{designPlot}}
  \tab A graphical representation of an experimental design \cr 
  \tab using labels stored in a matrix. \cr
  \tab It superseded design.plot.\cr
\code{\link{designPlotlabels}}
  \tab Plots labels on a two-way grid using \code{ggplot2}. \cr
\code{\link{designRandomize}}
  \tab Randomize allocated to recipient factors to produce  \cr 
  \tab a layout for an experiment. It supersedes \code{fac.layout}. \cr
\code{\link{designTwophaseAnatomies}}
  \tab Given the layout for a design and three structure formulae, \cr
  \tab obtain the anatomies for the (i) two-phase,  \cr 
  \tab (ii) first-phase, (iii) cross-phase, treatments, and \cr
  \tab (iv) combined-units designs. \cr
\code{\link{marginality.pstructure}} 
  \tab Extracts the marginality matrix from a \cr
  \tab \code{\link{pstructure.object}} \cr
\code{\link{marginality.pstructure}} 
  \tab Extracts a list containing the marginality matrices from \cr
  \tab a \code{\link{pcanon.object}} \cr
\code{\link{print.aliasing}} 
  \tab Prints an aliasing data.frame \cr
\code{\link{no.reps}}
  \tab Computes the number of replicates for an experiment\cr
\code{\link{detect.diff}}
  \tab Computes the detectable difference for an experiment\cr
\code{\link{power.exp}}
  \tab Computes the power for an experiment\cr
\code{\link{summary.pcanon}}
  \tab Summarizes the anatomy of a design, being the \cr
  \tab decomposition of the sample space based on its \cr
  \tab canonical analysis. \cr
\cr
(iv) ANOVA functions \cr\cr

\code{\link{fitted.aovlist}}
  \tab Extract the fitted values for a fitted model\cr
  \tab from an aovlist object\cr
\code{\link{fitted.errors}}
  \tab Extract the fitted values for a fitted model\cr
\code{\link{interaction.ABC.plot}}
  \tab Plots an interaction plot for three factors\cr
\code{\link{qqyeffects}}
  \tab Half or full normal plot of Yates effects\cr
\code{\link{resid.errors}}
  \tab Extract the residuals for a fitted model\cr
\code{\link{residuals.aovlist}}
  \tab Extract the residuals from an aovlist object\cr
\code{\link{strength}}
  \tab Generate paper strength values\cr
\code{\link{tukey.1df}}
  \tab Performs Tukey's\cr
  \tab one-degree-of-freedom-test-for-nonadditivity\cr
\code{\link{yates.effects}}
  \tab Extract Yates effects\cr
\cr
(v) Matrix functions\cr\cr

\code{\link{elements}}
  \tab Extract the elements of an array specified by\cr
  \tab the subscripts\cr
\code{\link{fac.ar1mat}}
  \tab Forms the ar1 correlation matrix for a\cr
  \tab (generalized) factor\cr
\code{\link{fac.sumop}}
  \tab Computes the summation matrix that produces\cr
  \tab sums corresponding to a (generalized) factor\cr
\code{\link{fac.vcmat}}
  \tab Forms the variance matrix for the variance\cr
  \tab component of a (generalized) factor\cr
\code{\link{mat.I}}
  \tab Forms a unit matrix\cr
\code{\link{mat.J}}
  \tab Forms a square matrix of ones\cr
\code{\link{mat.ar1}}
  \tab Forms an ar1 correlation matrix\cr
\code{\link{mat.ar2}}
  \tab Forms an ar2 correlation matrix\cr
\code{\link{mat.ar3}}
  \tab Forms an ar3 correlation matrix\cr
\code{\link{mat.arma}}
  \tab Forms an arma correlation matrix\cr
\code{\link{mat.banded}}
  \tab Forms a banded matrix\cr
\code{\link{mat.exp}}
  \tab Forms an exponential correlation matrix \cr
\code{\link{mat.gau}}
  \tab Forms a gaussian correlation matrix \cr
\code{\link{mat.ma1}}
  \tab Forms an ma1 correlation matrix\cr
\code{\link{mat.ma2}}
  \tab Forms an ma2 correlation matrix\cr
\code{\link{mat.ncssvar}}
  \tab Forms a variance matrix for random cubic \cr
  \tab smoothing spline effects\cr
\code{\link{mat.sar}}
  \tab Forms an sar correlation matrix\cr
\code{\link{mat.sar2}}
  \tab Forms an sar2 correlation matrix\cr
\code{\link{mat.dirprod}}
  \tab Forms the direct product of two matrices\cr
\code{\link{mat.dirsum}}
  \tab Forms the direct sum of a list of matrices\cr
\code{\link{mat.random}}
  \tab Calculates the variance matrix for the \cr
  \tab random effects from a mixed model, based  \cr
  \tab on a formula or a supplied matrix  \cr
\code{\link{mat.Vpred}}
  \tab Forms the variance matrix of predictions\cr
\code{\link{mat.Vpredicts}}
  \tab Forms the variance matrix of predictions,\cr
  \tab based  on supplied matrices or formulae. \cr
\code{\link{Zncsspline}}
  \tab Forms the design matrix for fitting the \cr
  \tab random effects for a natural cubic smoothing \cr
  \tab spline. \cr
\cr
(vi) Projector and canonical efficiency functions\cr\cr
\cr
Projector class:\cr\cr

\code{\link{projector}}
  \tab Create projectors\cr
\code{\link{projector-class}}
  \tab Class projector\cr
\code{\link{is.projector}}
  \tab Tests whether an object is a valid object of\cr
  \tab class projector\cr
\code{\link{print.projector}}
  \tab Print projectors\cr
\code{\link{correct.degfree}}
  \tab Check the degrees of freedom in an object of\cr
  \tab class projector\cr
\code{\link{degfree}}
  \tab Degrees of freedom extraction and replacement\cr
\cr
Accepts two or more formulae:\cr\cr

\code{\link{designAnatomy}}
  \tab An anatomy of a design, obtained from \cr
  \tab a canonical analysis of the relationships \cr
  \tab between sets of projectors.\cr
\code{\link{summary.pcanon}}
  \tab Summarizes the anatomy of a design, being the \cr
  \tab decomposition of the sample space based on its \cr
  \tab canonical analysis \cr
\code{\link{print.summary.pcanon}}
  \tab Prints the values in an 'summary.pcanon' object\cr
\code{\link{efficiencies.pcanon}}
  \tab Extracts the canonical efficiency factors from a \cr
  \tab list of class 'pcanon'\cr
\cr
Accepts exactly two formulae:\cr\cr

\code{\link{projs.2canon}}
  \tab A canonical analysis of the relationships between\cr
  \tab two sets of projectors\cr
\code{\link{projs.combine.p2canon}}
  \tab Extract, from a p2canon object, the projectors\cr
\code{\link{summary.p2canon}}
  \tab A summary of the results of an analysis of\cr
  \tab the relationships between two sets of projectors\cr
\code{\link{print.summary.p2canon}}
  \tab Prints the values in an 'summary.p2canon' object\cr
  \tab that give the combined decomposition\cr
\code{\link{efficiencies.p2canon}}
  \tab Extracts the canonical efficiency factors from\cr
  \tab a list of class 'p2canon'\cr
\cr
Accepts a single formula:\cr\cr

\code{\link{as.data.frame.pstructure}}
  \tab Coerces a pstructure.object to a data.frame \cr
\code{\link{print.pstructure}}
  \tab Prints a pstructure.object \cr
\code{\link{pstructure.formula}}
  \tab Takes a formula and constructs a pstructure.object \cr
  \tab that includes the orthogonalized projectors for the \cr
  \tab terms in a formula\cr

Others:\cr\cr

\code{\link{decomp.relate}}
  \tab Examines the relationship between the\cr
  \tab eigenvectors for two decompositions\cr
\code{\link{efficiency.criteria}}
  \tab Computes efficiency criteria from a set of\cr
  \tab efficiency factors\cr
\code{\link{fac.meanop}}
  \tab Computes the projection matrix that produces means\cr
\code{\link{proj2.eigen}}
  \tab Canonical efficiency factors and eigenvectors\cr
  \tab in joint decomposition of two projectors\cr
\code{\link{proj2.efficiency}}
  \tab Computes the canonical efficiency factors for\cr
  \tab the joint decomposition of two projectors\cr
\code{\link{proj2.combine}}
  \tab Compute the projection and Residual operators\cr
  \tab for two, possibly nonorthogonal, projectors\cr
\code{\link{show-methods}}
  \tab Methods for Function 'show' in Package dae\cr
\cr
(vii) Miscellaneous functions\cr\cr

\code{\link{extab}}
  \tab Expands the values in table to a vector\cr
\code{\link{get.daeTolerance}}
  \tab Gets the value of daeTolerance for the package dae\cr
\code{\link{harmonic.mean}}
  \tab Calcuates the harmonic mean.\cr
\code{\link{is.allzero}}
  \tab Tests whether all elements are approximately zero\cr
\code{\link{rep.data.frame}}
  \tab Replicate the rows of a data.frame by repeating  \cr 
  \tab each row consecutively and/or repeating all rows \cr
  \tab as a group. \cr
\code{\link{rmvnorm}}
  \tab Generates a vector of random values from a\cr
  \tab multivariate normal distribution\cr
\code{\link{set.daeTolerance}}
  \tab Sets the value of daeTolerance for the package dae\cr
}
}
\author{
\packageAuthor{dae}

Maintainer: \packageMaintainer{dae}
}
\keyword{ package }
