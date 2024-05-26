daeEnv <- new.env()
tol <- c(.Machine$double.eps ^ 0.5, .Machine$double.eps ^ 0.5)
names(tol) <- c("element.tol", "eigen.tol")
assign("daeTolerance",  tol, envir=daeEnv)
assign("daeRNGkind", "Mersenne-Twister", envir=daeEnv)

#Set the kind of random number generator for dae
set.daeRNGkind <- function(kind = "Mersenne-Twister")
{
  assign("daeRNGkind", kind, envir = daeEnv)
  invisible(kind)
}

#Get the current setting of the random number generator for dae
get.daeRNGkind <- function()
{ 
  kind <- get("daeRNGkind", envir=daeEnv)
  return(kind)
}

#Set element and eigen tolerances for dae
set.daeTolerance <- function(element.tol=NULL, eigen.tol=NULL)
{ 
  tol <- get("daeTolerance", envir=daeEnv)
  
  #Set value for element.tol, if one supplied
  if (!is.null(element.tol))
  { 
    if (!is.numeric(element.tol))
      stop("non-numeric argument")
    if (length(element.tol) != 1)
    { 
      element.tol <- element.tol[1]
      warning("More than one element.tol value supplied - only first used")
    } 
    tol[1] <- element.tol
  }
    
  #Set value for eigen.tol, if one supplied
  if (!is.null(eigen.tol))
  { 
    if (!is.numeric(eigen.tol))
      stop("non-numeric argument")
    if (length(eigen.tol) != 1)
    { 
      eigen.tol <- eigen.tol[1]
      warning("More than one eigen.tol value supplied - only first used")
    } 
    tol[2] <- eigen.tol
  }
  
  names(tol) <- c("element.tol", "eigen.tol")
  assign("daeTolerance", tol, envir = daeEnv)
  invisible(tol)
}

get.daeTolerance <- function()
{ 
  tol <- get("daeTolerance", envir=daeEnv)
  return(tol)
}

#function to test whether all elements are zero
#was all(abs(x) < daeTolerance)
"is.allzero" <- function(x)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
  mean(abs(x), na.rm=TRUE) < daeTolerance[["element.tol"]]
}


remove.repeats <- function(x, tolerance = 1E-06)
  # function to remove repeated values that differ by no more than daeTolerance  
{ #daeTolerance <- get("daeTolerance", envir=daeEnv)
  n <- length(x)
  if (n > 1)
  { repeats <-   c(FALSE, abs(x[2:n] - x[1:(n-1)]) < tolerance)
    x <- x[!repeats]
  }
  return(x)
}

"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ 
  kopt <- pmatch(arg.val, options)
  if (any(is.na(kopt)))
    stop("Value ",paste(arg.val, collapse = ","), " is either not unique or is not an allowed option for its argument")
  if (length(kopt) > 1)
  {
    warning(paste("Only one value allowed for argument where", 
                  paste(arg.val, collapse = ","), "have been supplied", 
                  sep = " "))
    kopt <- kopt[1]
  }
  return(kopt)
}

"separateFunction" <- function(var)
  #A function to separate the name of a function and the argument to the function
{ 
  #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("(", var, fixed=TRUE))!=0) 
  { 
    var <- (strsplit(var, "(", fixed=TRUE) )[[1]]
    var[2] <- (strsplit(var[2], ")", fixed=TRUE) )[[1]][1]
  }
  return(var)
}

"rmFunction" <- function(var, asreml.obj)
  #A function that returns the variable without any function
{ 
  var <- separateFunction(var)
  if (length(var)==2)
  { 
    var <- var[2]
    #Check for further arguments and strip, if found
    if (length(grep(",", var, fixed=TRUE))!=0) 
    { 
      var <- (strsplit(var, ",", fixed=TRUE) )[[1]]  
      var <- var[1]
    } 
  }  
  return(var)
}

"fac.getinTerm" <- function(term, rmfunction=FALSE)
  #function to return the set of factors/variables in a term separated by ':"
{ 
  if (length(term) != 1)
    stop("Multiple terms supplied where only one allowed")
  vars <- unlist(strsplit(term, ":", fixed=TRUE))
  if (rmfunction)
    vars <- unlist(lapply(vars, rmFunction))
  return(vars)
}
