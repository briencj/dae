#S3 methods

efficiencies <- function(object, ...) UseMethod("efficiencies")
marginality <- function(object, ...) UseMethod("marginality")
pstructure <- function(formula, ...) UseMethod("pstructure")


#"resid.errors" = function(object, ...) UseMethod("residuals")
#fitted.errors = function(x) UseMethod("fitted.errors")
#setGeneric("fitted.errors.aovlist")

#Deprecations

Ameasures <- function(...)
{ 
  .Deprecated(new = "designAmeasures", package = "dae")
  invisible()
}

blockboundary.plot <- function(...)
{ 
  .Deprecated(new = "blockboundaryPlot", package = "dae")
  invisible()
}

design.plot <- function(...)
{ 
  .Deprecated(new = "designPlot", package = "dae")
  invisible()
}

fac.layout <- function(...)
{ 
  .Deprecated(new = "designRandomize", package = "dae")
  invisible()
}

proj2.decomp <- function(...)
{ 
  .Deprecated(new = "proj2.eigen", package = "dae")
  invisible()
}

proj2.ops <- function(...)
{ 
  .Deprecated(new = "proj2.combine", package = "dae")
  invisible()
}

projs.canon <- function(...)
{ 
  .Deprecated(new = "designAnatomy", package = "dae")
  invisible()
}

projs.structure <- function(...)
{ 
  .Deprecated(new = "pstructure.formula", package = "dae")
  invisible()
}

