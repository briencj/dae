#"as.numfac" <- function(factor) {as.numeric(as.vector(factor))}

"as.numfac" <- function(factor)
{ if (!is.numeric(factor))
  { levs <- levels(factor)
    if (any(is.na(suppressWarnings(as.numeric(levs[!is.na(levs)])))))
      warning("Some levels do not have values that are numeric in nature")
    #see factor help
    factor <- as.numeric(levels(factor))[factor]
  }
  return(factor)
}

"mpone" <- function(factor) {2*as.numeric(factor)-3}

"fac.recode" <- function(factor, newlevels, ...)
#function to form a new factor by changing the levels of factor
{ 
  nlev <- length(levels(factor))
  if (nlev != length(newlevels))
  { stop("Must supply a new level for every level in the supplied factor")}
  new.fac <- factor(newlevels[factor], ...)
  return(new.fac)
}

"fac.recast" <- function(factor, newlevels = NULL, levels.order = NULL, newlabels = NULL, ...)
  #function to form a new factor by changing the levels of factor
{ 
  new.fac <- factor
  oldlevs <- levels(new.fac)
  if (!is.null(newlevels))
  {
    if (length(oldlevs) != length(newlevels))
      stop("Must supply a level for every level in the supplied factor")
    new.fac <- factor(newlevels[new.fac], ...)
    oldlevs <- levels(new.fac)
  }

  if (!(is.null(levels.order) & is.null(newlabels)))
  {
    if (is.null(levels.order))
      levels.order <- oldlevs
    if (!all(oldlevs %in% levels.order) | length(unique(levels.order)) < length(levels.order))
      stop("The set of levels.order must be unique and contain all the levels in the factor being recast")
    if (!is.null(newlabels))
      if (length(newlabels) != length(oldlevs) | length(unique(newlabels)) < length(newlabels))
        stop("The newlabels must be a set of unique values equal in length to the number of levels in the factor being recast")
    
    #Recast the factor
    if (is.null(newlabels))
      new.fac <- factor(new.fac, levels = levels.order, ...)
    else
      new.fac <- factor(new.fac, levels = levels.order, labels = newlabels, ...)
  }
  
  return(new.fac)
}

"fac.uselogical" <- function(x, levels = c(TRUE, FALSE), labels = c("yes", "no"), ...)
#function to form a two-level factor from a logical object
{
  if (length(levels) != 2)
    warning("length of supplied levels is not 2")
  if (length(labels) != 2)
    warning("length of supplied labels is not 2")
  if (!is.logical(x))
    x <- as.logical(x)
  fac <- factor(x, levels = levels, labels = labels, ...)
  return(fac)
}

"fac.combine" <- function(factors, order="standard", combine.levels=FALSE, sep=",", ...)
{
  #
  # test arguments
  #
  if (mode(factors) != "list") stop("Must supply a list")
  if (!all(sapply(factors, inherits, what="factor"))) 
    stop("All elements of list must be factors or ordereds")
  new.fac <- factors[[1]]
  nfac <- length(factors)
  if (nfac > 1)
  {
    if (var(sapply(factors, length)) != 0) stop("All factors must be of the same length")
    which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
    if (which.ord == "")	stop("order must be either standard or yates")
    # standard order
    if (which.ord == "1") counter <- nfac:1
    # Yates order
    else if (which.ord== "2") counter <- 1:nfac
    #
    # compute new factor
    #		
    radix <- 1
    new.fac <- rep(1, length(factors[[1]])) 
    na.lev <- FALSE
    for (i in counter)
      #reassign factor so unused levels removed
    { 
      f <- factor(factors[[i]])
      nlev <- length(levels(f))
      new.fac <- (as.numeric(f)-1) * radix + new.fac
      if (combine.levels)
      {
        if (any(is.na(levels(f))))
          na.lev <- TRUE
        if (i == counter[1])
          radix.lev <- paste(levels(f))
        else
        { if (which.ord == 1)
          radix.lev <- paste(rep(levels(f), each=radix), rep(radix.lev, times=nlev), sep=sep)
        else
          radix.lev <- paste(rep(radix.lev, times=nlev), rep(levels(f), each=radix), sep=sep)
        }
      }
      radix <- radix * nlev
    }
    if (combine.levels)
    { 
      obslevs <- unique(new.fac)
      obslevs <- obslevs[order(obslevs)]
      if (!na.lev)
        obslevs <- obslevs[!is.na(obslevs)]
      new.fac <- factor(new.fac, labels=radix.lev[obslevs], ...)
    }	
    else
      new.fac <- factor(new.fac, ...)
  }
  return(new.fac)
}

#Function the takes a factor whose levels include a separator that can be used to cleave 
#each level into the levels for several new factors.
fac.uncombine <- function(factor, new.factors, sep = ",", ...)
{
  if (mode(new.factors) != "list") stop("Must supply a list for new.factors")
  if (!inherits(factor, what = "character"))
  {
    if (inherits(factor, what = "factor"))
      factor <- as.character(factor)
    else
      stop("factor is not a factor or character")
  }
  fac.list <- strsplit(factor, split = sep)
  #name the values in each elements of fac.list 
  fac.list <- lapply(fac.list, 
                     function(obs, new.factors)
                     {
                       names(obs) <- names(new.factors)
                       if (length(obs) != length(new.factors))
                         stop("Not all levels of factor have separated into ",length(new.factors)," levels")
                       return(obs)
                     }, new.factors = new.factors)
  #form separate columns for the values in each element of fac.list
  new.facs <- lapply(names(new.factors), 
                     function(elem)
                     {
                       unlist(lapply(fac.list, function(obs, elem) obs[elem], 
                                     elem = elem))
                     })
  names(new.facs) <- names(new.factors)
  df <- as.data.frame(new.facs, stringsAsFactors = FALSE)
  #Make factors
  for ( fac in names(new.factors))
  {
    if (!is.null(new.factors[[fac]]))
      df[fac] <- factor(df[[fac]], levels = new.factors[[fac]])
    df[fac] <- factor(df[[fac]], ...) #make sure no unused levels
  }
  return(df)
}

"fac.split" <- function(combined.factor, factor.names, sep = ",", ...)
{
  #
  # test arguments
  #
  if (!inherits(combined.factor, what = "factor")) 
    stop("combined.factor must be a factor")
  if (mode(factor.names) != "list") 
    stop("Must supply a list of factor names")
  
  vals <- as.character(combined.factor)
  nfac <- length(factor.names)
  vals <- strsplit(vals, sep, fixed = TRUE)
  kfac <- length(vals[[1]])
  if (nfac != kfac)
    stop("The number of names in factor.names does no match the number of comma-separated strings in combined.factor")
  
  #Extract the values for each factor
  dat <- as.data.frame(do.call(cbind,lapply(1:nfac,
                                            function(k)
                                            {
                                              unlist(lapply(vals,
                                                            function(val,k )val[k],
                                                            k = k))
                                            })))
  
  #Form the factors
  dat <- mapply(dat, factor.names, 
                FUN = function(fac, nam) 
                {
                  if (!is.null(nam))
                    fac <- factor(fac, levels = nam, ...)
                  else
                    fac <- factor(fac, ...)
                  return(fac)
                }, SIMPLIFY = FALSE)
  dat <- as.data.frame(dat)
  names(dat) <- names(factor.names)
  return(dat)
}
  

"fac.divide" <- function(combined.factor, factor.names, order="standard")
{
#
# test arguments
#
#	if (mode(factor.names) != "list" | length(factor.names) == 1) 
#    stop("Must supply a list of more than one component")
	if (mode(factor.names) != "list") 
    stop("Must supply a list of factor names")
	nfac <- length(factor.names)
	which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1") counter <- nfac:1
# Yates order
	else if (which.ord== "2") counter <- 1:nfac
#
# convert factor or numeric to numeric with values 1:nlev and compute new factors
#		
	if (is.factor(combined.factor))
	  nlev <- length(levels(combined.factor))
	else
	  nlev <- length(unique(combined.factor))
	kombined.factor <- c(1:nlev)[combined.factor]
	for (i in counter)
	{ 
	  klev <- length(factor.names[[i]])
    if (is.numeric(factor.names[[i]]))
    { 
      if (klev != 1)
      { 
        lev <- factor.names[[i]]
        val <- (kombined.factor - 1) %% klev + 1
        factor.names[[i]] <- factor(lev[val])
      }
      else
      {	
        klev <- factor.names[[i]]
		    factor.names[[i]] <- factor((kombined.factor - 1) %% klev + 1)
	    }
    }
    else
	  {	
	    factor.names[[i]] <- factor((kombined.factor - 1) %% klev + 1, 
                                          labels = factor.names[[i]])
	  }
		kombined.factor <- (kombined.factor - 1) %/% klev + 1
	}
  new.factors <- data.frame(factor.names)
}

"fac.nested" <- function(nesting.fac, nested.levs = NA, nested.labs = NA, ...)
{
  #Deal with deprecated levels and labels function arguments
  tempcall <- list(...)
  if (length(tempcall))
  {
    if ("levels" %in% names(tempcall))
      stop("levels has been deprecated in fac.nested - use nested.levs")
    if ("labels" %in% names(tempcall))
      stop("labels has been deprecated in fac.nested - use nested.labs")
  }
  nested.fac <- split(nesting.fac, nesting.fac)
  nested.fac <- lapply(nested.fac, function(nesting.lev) (1:length(nesting.lev)))
  no.lev.within <- max(unlist(lapply(nested.fac, length)))
  nested.fac <- unsplit(nested.fac, nesting.fac)
  if (length(nested.levs) == 1 && is.na(nested.levs)) nested.levs <- 1:no.lev.within
	if (length(nested.labs) == 1 && is.na(nested.labs)) nested.labs <- as.character(nested.levs)
	nested.fac <- factor(nested.fac, levels=nested.levs, labels=nested.labs, ...)
	return(nested.fac)
}

"fac.multinested" <- function(nesting.fac, nested.fac = NULL, fac.prefix = NULL, 
                              nested.levs = NA, nested.labs = NA, outlevel = 0, outlabel = "rest", ...)
{
  n <- length(nesting.fac)
  levs <- na.omit(levels(nesting.fac))
  levs <- levs[levs != outlabel]
  no.lev.between <- length(levs)

  if (is.null(nested.fac))
  {
    no.lev.within <- max(table(nesting.fac))
    reps <- table(nesting.fac)[levs]
    names(reps) <- levs
    nest.facs <- lapply(levs, 
                        function(lev, nesting.fac, reps, nested.levs, nested.labs, outlevel, outlabel) 
                        {
                          nest.fac <- rep(outlevel, n)
                          nest.fac[is.na(nesting.fac)] <-NA
                          #Set levels if nested factor for each levels of nested factor, avoiding missing values
                          if (length(nested.levs) == 1 && is.na(nested.levs)) 
                            klevs <- 0:reps[lev]
                          else
                            klevs <- c(outlevel,nested.levs[1:reps[lev]])
                          if (length(nested.labs) == 1 && is.na(nested.labs)) 
                            klabs <- as.character(c(outlabel, klevs[-1]))
                          else
                            klabs <- as.character(c(outlabel, nested.labs[1:reps[lev]]))
                          nest.fac[!is.na(nesting.fac)][nesting.fac[!is.na(nesting.fac)] == lev] <- 
                            klevs[-1][1:reps[lev]]
                          nest.fac <- factor(nest.fac, levels=klevs, labels=klabs, ...)
                          return(nest.fac)
                        }, nesting.fac = nesting.fac, reps = reps, nested.levs=nested.levs, nested.labs=nested.labs, 
                           outlevel = outlevel, outlabel = outlabel)
  } else
  {
    nested.fac <- as.character(nested.fac)
    nest.facs <- lapply(levs, 
                        function(lev, nesting.fac, nested.fac, nested.labs, outlevel, outlabel) 
                        {
                          nest.fac <- rep(outlevel, n)
                          nest.fac[is.na(nesting.fac) | is.na(nested.fac)] <-NA
                          #Set levels if nested factor for each levels of nested factor, avoiding missing values
                          nest.fac[!is.na(nesting.fac)][nesting.fac[!is.na(nesting.fac)] == lev] <- 
                            nested.fac[!is.na(nesting.fac)][nesting.fac[!is.na(nesting.fac)] == lev]
                          klevs <- nested.fac[!is.na(nesting.fac)][nesting.fac[!is.na(nesting.fac)] == lev]
                          klevs <- c(outlevel,klevs)
                          if (length(nested.labs) == 1 && is.na(nested.labs)) 
                            klabs <- as.character(c(outlabel, klevs[-1]))
                          else
                            klabs <- as.character(c(outlabel, nested.labs))
                          nest.fac <- factor(nest.fac, levels=klevs, labels=klabs, ...)
                          return(nest.fac)
                        }, nesting.fac = nesting.fac, nested.fac = nested.fac, nested.labs=nested.labs, 
                           outlevel = outlevel, outlabel = outlabel)
  }
  nest.facs <- data.frame(nest.facs)
  names(nest.facs) <- paste0(fac.prefix, levs)
  return(nest.facs)
}

