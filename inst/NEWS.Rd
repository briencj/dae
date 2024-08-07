\name{news}
\title{News for Package dae}
\section{Changes in dae version 3.2.29 (2024-07-23)}{
  \itemize{
  \item Add arguments \code{center} and \code{scale} to as.numfac to centre and scale the numeric.
  \item Remove links to functions in other libraries. 
}}
\section{Changes in dae version 3.2.28 (2024-06-12)}{
  \itemize{
  \item Update for submission to CRAN.
}}
\section{Changes in dae version 3.2.27 (2024-05-27)}{
  \itemize{
  \item Add the variable \code{daeRNGkind} to the \code{daeEnv}, with default \code{Mersenne-Twister}, and the functions set.daeRNGkind and get.RNGkind. This ensures tha the random number generator used by dae has to be explicitly changed by the user.
}}
\section{Changes in dae version 3.2.26 (2024-04-21)}{
  \itemize{
  \item Fixed bugs in testing whether projectors are symmetric or idempotent caused by unimportant differences in attributes of objects involved in the testing.
}}
\section{Changes in dae version 3.2.25 (2024-04-09)}{
  \itemize{
  \item Update for submission to CRAN.
}}
\section{Changes in dae version 3.2.24 (2024-03-19)}{
  \itemize{
  \item Fixed bug in pstructure, also affecting designAnatomy, that resulted aliasing of two terms going undetected, nesting of sources was not being correctly determined when aliasing was present, and computed marginality relationships not being transitive.
}}
\section{Changes in dae version 3.2.23 (2024-02-25)}{
  \itemize{
  \item Added trap for svd error so that only a warning will be given and an NA returned. 
}}
\section{Changes in dae version 3.2.22 (2024-02-13)}{
  \itemize{
  \item Add link to GitHub issues page for BugReports in DESCRIPTION file. 
  \item Fixed problem in mat.ginv for inverting matrices of rank one that led to a conformability error caused by R returning a single column matrix as a vector. 
}}
\section{Changes in dae version 3.2.21 (2024-01-14)}{
  \itemize{
    \item Corrects some minor issues in the manual.
}}
\section{Changes in dae version 3.2.20 (2024-01-14)}{
  \itemize{
    \item Update the \code{DesignNotes} vignette.
    \item Add the function porthogonalize.list to orthogonalize a sequence of \code{projectors}. It externalizes the process for othogonalizing projectors previously performed within pstructure.formula.
    \item Add the function mat.ginv for computing generalized inverses.
    \item Add the functions mat.cor and mat.corg for forming correlation matrices.
}}
\section{Changes in dae version 3.2.19 (2023-08-08)}{
  \itemize{
    \item Fix CRAN errors and increment version number.
}}
\section{Changes in dae version 3.2.18 (2023-08-07)}{
  \itemize{
    \item Correct the DesignNotes vignette.
    \item Add \code{facetstrips.switch} argument to designGGPlot and the \code{facetstrips.placement} argument to designGGPlot and designBlocksGGPlot.
}}
\section{Changes in dae version 3.2.17 (2023-03-30)}{
  \itemize{
    \item Correct the DesignNotes vignette.
    \item Add a check for a term that is equivalent to the grand mean term in designAnatomy.
}}
\section{Changes in dae version 3.2.16 (2023-03-26)}{
  \itemize{
    \item Revise the DesignNotes vignette.
}}
\section{Changes in dae version 3.2.15 (2023-03-10)}{
  \itemize{
    \item Update ggplot2 functions to use tidy evaluation.
}}
\section{Changes in dae version 3.2-14 (2023-01-13)}{
  \itemize{
    \item Update \code{DesignNotes} vignette.
    \item Add a \code{replications} argument to designAmeasures and implement the calculation of weighted average pairwise variances using the \code{replications}.
}}
\section{Changes in dae version 3.2-13 (2021-11-15)}{
  \itemize{
    \item Rebuild for R version 4.1.2.
}}
\section{Changes in dae version 3.2-12 (2021-11-05)}{
  \itemize{
    \item Fix bug in the function fac.multinested that resulted in the order of the levels for each new nested factor, when each has the same levels, depending on the order that they occurred within the level of the \code{nesting.fac}.
}}
\section{Changes in dae version 3.2-11 (2021-09-01)}{
  \itemize{
    \item Add the function fac.uncombine.
}}
\section{Changes in dae version 3.2-10 (2021-08-6)}{
  \itemize{
    \item Add the functions fac.split and fac.genfactors.
    \item Modify fac.gen so that the number of levels for each recipient factor is the number of observed levels of that factor.
    \item Fix bug in finding marginal terms in pstructure.formula when there is partial aliasing amongst terms.
}}
\section{Changes in dae version 3.1-41 (2021-06-07)}{
  \itemize{
    \item Add the function fac.recast, a function that replaces and extends fac.recode.
}}
\section{Changes in dae version 3.1-40 (2021-05-12)}{
  \itemize{
    \item Fix bug in pstructure that resulted in errors in aliasing reporting.
}}
\section{Changes in dae version 3.1-39 (2021-04-10)}{
  \itemize{
    \item Add \code{scales.free} argument to allow the axis scales to vary between facets in designGGPlot.
}}
\section{Changes in dae version 3.1-38 (2021-03-19)}{
  \itemize{
    \item Rebuild vignette
}}
\section{Changes in dae version 3.1-37 (2021-03-17)}{
  \itemize{
    \item Recompile for submission to CRAN
}}
\section{Changes in dae version 3.1-36 (2021-03-04)}{
  \itemize{
    \item Fix bug in pstructure for parenthesizing terms that consist of generalized or combined factors.
}}
\section{Changes in dae version 3.1-35 (2021-02-25)}{
  \itemize{
    \item Fix bug in the aliasing table, where "## Information remaining" not printed for pstructure.
    \item Add \code{...} argument to print.aliasing call in print.structure.
    \item Allow multiple values for \code{which} argument in print.pstructure.
}}
\section{Changes in dae version 3.1-34 (2021-02-16)}{
  \itemize{
    \item Minor refinements of the the tables of aliasing generated by pstructure.formula, projs.2canon and designAnatomy.
}}
\section{Changes in dae version 3.1-33 (2021-02-15)}{
  \itemize{
    \item Tidy up print.pstructure by removing the ineffective \code{which.criteria} option and adding a \code{which} option to control which components are printed.
    \item Change the tables of aliasing generated by pstructure.formula, projs.2canon and designAnatomy so that the amount of information that is confounded is reported, instead of the amount of information that remains. The amount that remains is reported in the decomposition table.
}}
\section{Changes in dae version 3.1-32 (2020-11-11)}{
  \itemize{
    \item Rebuild vignettes and and package for R 4.0.3 and CRAN submission.
}}
\section{Changes in dae version 3.1-31 (2020-09-21)}{
  \itemize{
    \item Fix a bug in print.summary.pcanon that resulted in incorrect consecutive repetitions of sources within a column.
}}
\section{Changes in dae version 3.1-30 (2020-09-15)}{
  \itemize{
    \item Check that at least one of \code{labels} and \code{cellfillcolour.column} is set in a call to designGGPlot.
}}
\section{Changes in dae version 3.1-29 (2020-08-05)}{
  \itemize{
    \item Add \code{label.size} and \code{labeller} arguments to designGGPlot.
}}
\section{Changes in dae version 3.1-28 (2020-07-20)}{
  \itemize{
    \item Re-implement fac.nested so that it is much more efficient.
}}
\section{Changes in dae version 3.1-27 (2020-07-09)}{
  \itemize{
    \item Rebuild vignettes and and package for R 4.0.2 and CRAN submission.
}}
\section{Changes in dae version 3.1-26 (2020-06-15)}{
  \itemize{
    \item Add \code{x.axis.position} argument to designGGPlot, making \code{top} the default.
}}
\section{Changes in dae version 3.1-25 (2020-05-12)}{
  \itemize{
    \item Rebuild for R version 4.0.0.
    \item Minor revisions to \code{DesignNotes} vignette.
}}
\section{Changes in dae version 3.1-24 (2020-03-23)}{
  \itemize{
    \item Change the Maintainer email address.
    \item Add function mat.vcomp.
}}
\section{Changes in dae version 3.1-23 (2020-03-16)}{
  \itemize{
    \item Add the function rep.data.frame.
    \item Rename the design species for which designTwoPhaseAnatomies produces anatomies. Also, allow the overriding of the default titles for each anatomy.
}}
\section{Changes in dae version 3.1-22 (2020-01-16)}{
  \itemize{
    \item Add the function fac.uselogical.
}}
\section{Changes in dae version 3.1-21 (2020-01-12)}{
  \itemize{
    \item Reduce example for fac.multinested and build for CRAN resubmission.
}}
\section{Changes in dae version 3.1-20 (2020-01-11)}{
  \itemize{
    \item Rebuild for CRAN submission.
}}
\section{Changes in dae version 3.1-19 (2019-12-30)}{
  \itemize{
    \item Add names of countries to Countries factor in Casuarina.dat.
    \item Make additions to the DesignNotes vignette.
}}
\section{Changes in dae version 3.1-18 (2019-12-27)}{
  \itemize{
    \item Add names of missing factor/covariate to error message in pstructure.
    \item Add function fac.multinested
}}
\section{Changes in dae version 3.1-17 (2019-11-23)}{
  \itemize{
    \item Revise the DesignNotes vignette.
    \item Further debug designRandomize for some \code{data.frame}s not in standard order.
}}
\section{Changes in dae version 3.1-16 (2019-11-07)}{
  \itemize{
    \item Correct typo and re-upload to CRAN.
}}
\section{Changes in dae version 3.1-15 (2019-11-07)}{
  \itemize{
    \item Prepare for upload to CRAN.
}}
\section{Changes in dae version 3.1-14 (2019-10-29)}{
  \itemize{
    \item Add the capability to mat.Vpredicts to return either the variance or the information matrix. Also, add the rank of the matrix as an attribute of the matrix.
}}
\section{Changes in dae version 3.1-13 (2019-10-27)}{
  \itemize{
    \item Revise DesignNotes vignette.
    \item Further debug designRandomize for either \code{data.frame} not in standard order or column order inconsistencies.
    \item Rename name of data set \code{Exp249.mplot.des} to \code{Exp249.munit.des}.
}}
\section{Changes in dae version 3.1-12 (2019-10-20)}{
  \itemize{
    \item Include check for names supplied in \code{nested.recipients} are in the names in  \code{recipients}.
}}
\section{Changes in dae version 3.1-11 (2019-10-10)}{
  \itemize{
    \item Fix bug in designRandomize that occurs when a supplied \code{data.frame} is not in standard order.
    \item Deprecate fac.layout
    \item Update the DesignNotes vignette.
    }}
\section{Changes in dae version 3.1-10 (2019-09-29)}{
  \itemize{
    \item Merge features of designPlot and designPlotlabels into a new function designGGPlot and replace blockboundaryPlot with designBlocksGGPlot. The functions designPlot, blockboundaryPlot and designPlotlabels have been marked for future deprecation.
    \item Update the DesignNotes vignette.
    }}
\section{Changes in dae version 3.0-32 (2019-05-15)}{
  \itemize{
    \item Rebuild under version 3.6.0 prior to resubmission to CRAN.
    }}
\section{Changes in dae version 3.0-31 (2019-04-04)}{
  \itemize{
    \item Update the Author field to include an ORCID and the Description field of the DESCRIPTION.
    \item Add new function designTwophaseAnatomies.
}}
\section{Changes in dae version 3.0-30 (2019-03-29)}{
  \itemize{
    \item Add \code{reverse.x} and \code{reverse.y} arguments to designPlotlabels.
}}
\section{Changes in dae version 3.0-29 (2019-03-23)}{
  \itemize{
    \item Dcoumentation revisions.
}}
\section{Changes in dae version 3.0-28 (2019-03-12)}{
  \itemize{
    \item Fix bug in mat.Vpredicts when there are no random effects.
}}
\section{Changes in dae version 3.0-27 (2019-02-21)}{
  \itemize{
    \item Use R.rsp to install vignettes so that hyperlinks are active.
}}
\section{Changes in dae version 3.0-26 (2019-02-18)}{
  \itemize{
    \item Deal with a mean only formula (~1) in pstructure.formula.
}}
\section{Changes in dae version 3.0-25 (2019-02-07)}{
  \itemize{
    \item Reorganize vignette and include manual as a vignette.
}}
\section{Changes in dae version 3.0-24 (2019-02-01)}{
  \itemize{
    \item Minor updates to the daeDesignNotes vignette and move to a vignettes subdirectory.
}}
\section{Changes in dae version 3.0-23 (2018-10-17)}{
  \itemize{
    \item Add diagnostic message to mat.Vpredicts.
    \item Fix a bug in setValidity of class "projector".
}}
\section{Changes in dae version 3.0-22 (2018-10-01)}{
  \itemize{
    \item Fix bug in mat.Vpredicts when have a a single random term for a single factor.
}}
\section{Changes in dae version 3.0-21 (2018-09-29)}{
  \itemize{
    \item Add the new function mat.Vpredicts that has renamed arguments and the capability of using formulae in generating the required matrices.
}}
\section{Changes in dae version 3.0-20 (2018-09-02)}{
  \itemize{
    \item Add \code{eliminate} argument to mat.Vpred.
}}
\section{Changes in dae version 3.0-19 (2018-08-25)}{
  \itemize{
    \item Fix bug in summary.pcanon for a single \code{formula} when \code{projectors} are saved.
}}
\section{Changes in dae version 3.0-18 (2018-08-09)}{
  \itemize{
    \item Small documentation edits.
}}
\section{Changes in dae version 3.0-17 (2018-6-12)}{
  \itemize{
    \item Rebuild under latest version of R 3.5.0 with requirement for R >= 3.1.0.
}}
\section{Changes in dae version 3.0-16 (2018-6-11)}{
  \itemize{
    \item Rebuild under latest version of R 3.6.0.
}}
\section{Changes in dae version 3.0-15 (2018-04-06)}{
  \itemize{
    \item Fix bugs in designAnatomy when there is only one structure.
    \item Fix bugs in print.summary.pcanon when there are no entries in a column of the table.
    \item Add a warning to print.summary.canon when the sum of  sources' dimensions is less than the number of rows in data.
}}
\section{Changes in dae version 3.0-14 (2017-12-24)}{
  \itemize{
    \item Modified checking of options for function arguments.
}}
\section{Changes in dae version 3.0-13 (2017-12-08)}{
  \itemize{
    \item Modified mat.Vpred so that the Moore-Penrose inverse, which is unique, is always used for matrix inverses.
}}
\section{Changes in dae version 3.0-12 (2017-12-04)}{
  \itemize{
    \item Rebuild under latest version of R.
}}
\section{Changes in dae version 3.0-11 (2017-10-22)}{
  \itemize{
    \item Minor editing of documentation.
}}
\section{Changes in dae version 3.0-10 (2017-10-08)}{
  \itemize{
    \item Changed designRandomize so that, if \code{allocated} is \code{NULL}, then a random permutation of the recipient factors is produced.
}}
\section{Changes in dae version 3.0-09 (2017-09-28)}{
  \itemize{
    \item Adjust caption of summary.pcanon.
    \item Added new function designPlotlabels that uses ggplot to plot labels on a grid.
}}
\section{Changes in dae version 3.0-08 (2017-09-20)}{
  \itemize{
    \item Correct bug in axis labelling in interaction.ABC.plot.
}}
\section{Changes in dae version 3.0-07 (2017-07-31)}{
  \itemize{
    \item Reinstate ability to form and print a summary from old pcannon objects (pre version 3.0-0).
}}
\section{Changes in dae version 3.0-06 (2017-07-26)}{
  \itemize{
    \item Fixed a bug in summary.pcanon so that aliasing table is printed when required.
}}
\section{Changes in dae version 3.0-05 (2017-07-13)}{
  \itemize{
    \item Put in a check that factors/covariates in a formula are in data.
    \item Formatted the aliasing table to be similar to the decomposition table.
}}
\section{Changes in dae version 3.0-04 (2017-07-05)}{
  \itemize{
    \item Apply David Stern's fix to fac.nested so that it deals with missing values in the nesting factor.
}}
\section{Changes in dae version 3.0-02 (2017-06-30)}{
  \itemize{
    \item Added links to documentation.
}}
\section{Changes in dae version 3.0-01 (2017-06-29)}{
  \itemize{
    \item Added a component named \code{aliasing} to \code{pstructure.object} that stores the information about (partial) aliasing between sources derived a single \code{formula}. Added a component named \code{aliasing} to  \code{pcanon.object} that collects all the aliasing that has occured in the decomposition.
    \item Added the function print.aliasing to print a data.frame of class \code{aliasing} that contains aliasing information.
    \item Renamed Ameasures to designAmeasures
    \item Replaced the vignette for randomzing designs with a vignette that discusses the facilities for designing experiments that includes the use of designRandomize and designAnatomy.
}}
\section{Changes in dae version 3.0-0 (2017-06-11)}{
  \itemize{
    \item Added the new function designRandomize with the intention that it will eventaully replace fac.layout. The new function replaces the arguments \code{unrandomized} and \code{randomized} with \code{recipients} and \code{allocated}, respectively, and \code{nested.factors} with \code{nested.recipients}. Also, it changes the default for \code{unit.permutation} to \code{FALSE}.
    \item Deprecated projs.structure, replacing it with pstructure.formula, creating a pstructure.object that stores the projectors, the lists of their term names and source names, and the marginality matrix for the terms. Added arguments \code{marginality} and \code{check.marginality} to allow the marginality matrix to be supplied, overwriting any that is produced by pstructure.formula. The argument \code{check.marginality} allows one to decide if the supplied and calculated marginality matrices are the same. Methods for as.data.frame, marginality and print were added for objects of class \code{pstructure}.
    \item Deprecated projs.canon, it being replaced by designAnatomy.
    \item Modified the function designAnatomy to accomodate the changes to pstructure. Added new argument \code{labels} to designAnatomy, which has option \code{sources} that will result in sources indicating whether they represent interactions and/or nested effects. Added argument \code{marginality} that allows a list of marginality matrices to be supplied. Added the argument \code{check.marginality} allows one to decide if the supplied and calculated marginality matrices are the same. 
    \item Added argument \code{labelswap} to summary.pcanon so that can swap between terms and sources in the output.
    \item Added new method marginality.pcannon to extract the list of marginality matrices corresponding to the \code{formulae} in the call to designAnatomy. Added the argument \code{labelswap} to the summary.pcanon method so that either terms or sources can be output.
    \item Added the oats split-plot design from Yates (1937).
  }
}
\section{Changes in dae version 2.7-26 (2017-05-04)}{
  \itemize{
    \item Converted projs.canon to designAnatomy, changing the position of \code{data} to be second in the argument list and changing \code{meanTerm} to \code{grandMean}. Retained projs.canon as a synonym with the same arguments as before. Also changed \code{meanTerm} to \code{grandMean} in projs.structure.
    \item Added new functions mat.ncssvar and Zspline that calculate the varaince matrix of the random spline effects and the design matrix for fitting random spline effects, respectively.
  }
}
\section{Changes in dae version 2.7-25 (2017-04-12)}{
  \itemize{
    \item In designPlot, have the whole \code{designMatrix} plotted if \code{labels} is \code{NULL}.
  }
}
\section{Changes in dae version 2.7-24 (2017-02-10)}{
  \itemize{
    \item Add the more general \code{hybrid} option to projs.structure and  projs.canon and made it the default option.
    \item Changed fac.combine so that a single factor would pass through without generating an error message.
    \item Revise interaction.ABC.plot to account for ggplot deprecation.
    \item Add functions mat.sar, to form the correlation matrix from an SAR parameter, and mat.gau, to form a correlation function based on the square of the distance between coordinates.
  }
}
\section{Changes in dae version 2.7-23 (2017-01-05)}{
  \itemize{
    \item Fix a bug in tukey.1df when there is no Error function in the model.
    \item Add tips.
  }
}
\section{Changes in dae version 2.7-22 (2016-10-03)}{
  \itemize{
    \item Fix bug in no.reps when the number of reps to achieve the desired power is less than 2.
  }
}
\section{Changes in dae version 2.7-21 (2016-09-19)}{
  \itemize{
    \item Fix bug in fac.combine that can introduce a level equal to \code{NA} when none of the factors has such a level.
  }
}
\section{Changes in dae version 2.7-20 (2016-09-16)}{
  \itemize{
    \item Rebuild under R 3.4.0.
  }
}
\section{Changes in dae version 2.7-19 (2016-09-07)}{
  \itemize{
    \item Include the argument \code{meanTerm} in projs.structure and projs.canon to allow the option to include the grand mean term.
  }
}
\section{Changes in dae version 2.7-18 (2016-08-07)}{
  \itemize{
    \item Fixed a bug for aliased terms in projs.jandw.
  }
}
\section{Changes in dae version 2.7-17 (2016-05-18)}{
  \itemize{
    \item Extended Ameasures to include a \code{groups} argument that allows explicit specification of the elements of groups for which A values are required.
  }
}
\section{Changes in dae version 2.7-16 (2016-05-16)}{
  \itemize{
    \item Fixed a bug in fac.gen when there is a factor with a single character value.
  }
}
\section{Changes in dae version 2.7-15 (2016-05-14)}{
  \itemize{
    \item Revised fac.recode and fac.divide to deal with factors with both numeric and character levels and to improve efficiency.
  }
}
\section{Changes in dae version 2.7-14 (2016-05-10)}{
  \itemize{
    \item Added designLatinSqrSys. 
    \item Revised some documentation.
  }
}
\section{Changes in dae version 2.7-13 (2016-04-25)}{
  \itemize{
    \item Added functions mat.Vpred and Ameasures. 
    \item Fixed bugs in fac.vcmat. 
  }
}
\section{Changes in dae version 2.7-12 (2016-03-10)}{
  \itemize{
    \item Reversed order of parameters in mat.ar1 and mat.exp. 
    \item Added functions mat.ar2, mat.ar3, mat.ma1, mat.ma2, mat.arma and mat.sar2. 
  }
}
\section{Changes in dae version 2.7-11 (2016-02-11)}{
  \itemize{
    \item Added \code{xefficiency}, the maximum efficiency factor projector-based functions. 
  }
}
\section{Changes in dae version 2.7-10 (2016-01-19)}{
  \itemize{
    \item Improved the code for analyzing the relationships between pairs of projectors. 
  }
}
\section{Changes in dae version 2.7-9 (2015-01-03)}{
  \itemize{
    \item Include new functions designPlot and blockboundaryPlot, that are the same as design.plot and blockboundary.plot, except that they have less cryptic argument names and some arguments have been modified. It is envisaged that the new functions will supersede the old functions with the old functions being deprecated in due course.
    \item Fixed bug in projs.structure to ensure that arguments, including \code{keep.order}, are passed to the terms function. 
  }
}
\section{Changes in dae version 2.7-8 (2015-12-16)}{
  \itemize{
    \item Add arguments rchtdiv and cchtdiv to design.plot so that can size row and column labels independetly of the cell labels. 
  }
}
\section{Changes in dae version 2.7-7 (2015-11-15)}{
  \itemize{
    \item Improved the documentation of efficiency criteria in a number of functions. 
  }
}
\section{Changes in dae version 2.7-6 (2015-11-02)}{
  \itemize{
    \item Allow a single formula to be specified in formula of projs.canon.
    \item Improved proj2.combine to equate Qconf to Q1 when all of Q1 is confounded by Q2. 
  }
}
\section{Changes in dae version 2.7-5 (2015-10-29)}{
  \itemize{
    \item Fix bug in mat.banded in which the elements were being placed starting with the first subdiagonal instead of the diagonal. 
  }
}
\section{Changes in dae version 2.7-4 (2015-10-01)}{
  \itemize{
    \item Add arrgument \code{ggplotFuncs} to interaction.ABC.plot to allow \code{ggplot} functions to be used to 
    modify the plot. 
    \item Add detect.diff to allow the calculation of the detectable difference for given replication, power, delta and alpha.
  }
}
\section{Changes in dae version 2.7-3 (2015-08-23)}{
  \itemize{
    \item Change option \code{differencingfirst} to \code{differencing} in projs.structure and projs.canon. 
          Code modified so that differencing only orthogonalizes by subtracting from a projector those 
          projectors whose factors are a subset of its factors. 
    \item Added new criteria \code{dforth} to efficiency.criteria.
    \item Added \code{order} to the default setting of \code{which.criteria} in all functions that have this argument.
    \item NEWS added.
  }
}
\section{Changes in dae version 2.7-2 (2015-06-24)}{
  \itemize{
    \item minor reprogramming.
    \item documentation of deprecated functions added.
    \item NEWS added.
  }
}
\section{Changes in dae Version 2.7-1 (2015-05-26)}{
  \itemize{
    \item An except argument added to fac.layout so that nominated unrandomized factors are not
          permuted.
  }
}
\section{Changes in dae Version 2.7-0 (2015-04-29)}{
  \itemize{
    \item Added functions to perform an eigenalysis of the projection matrices for an experimental  
          design to produce the decomposition or skeleton-anova table. These functions include 
          projs.combine.p2canon, projs.structure, projs.2canon, projs.canon, summary.p2canon, 
          summary.pcanon, print.summary.p2canon, print.summary.pcanon. The top-level functions 
          are those with the suffix 'pcanon'.
    \item This version uploaded to CRAN.
  }
}
