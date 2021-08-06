"factor.list" <- function(generate, order="standard")
  #takes a generate list and creates a list of factor names, with levels 
  #information, and a list of factor relative replications, both of which are 
  #returned as a list of the two parallel lists called factors and reps.
{ 
  n <- length(generate)
  which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
  if (which.ord == "")	stop("order must be either standard or yates")
  if (which.ord == "1") # standard order
    counter <- 1:n
  else                  # Yates order
    counter <- n:1
  kfac <- 0
  for(i in counter) 
  { 
    if (!(names(generate[i]) == ""))
    { 
      kfac=kfac+1
      if (kfac == 1)
      { 
        fnames <- list(generate[[i]])
        names(fnames) <- names(generate)[i]
        freps <- 1
      }
      else
      { 
        knames <- list(generate[[i]])
        names(knames) <- names(generate)[i]
        fnames <- c(fnames, knames)
        freps <- c(freps, 1)
      }
    }
    else
    { 
      if (kfac == 0)
        if (which.ord == "1")
          stop("Must start with a factor name - set times argument instead")
      else
        stop("Must end with a factor name - set each argument instead")
      freps[kfac] <- generate[[i]]
    }
  }
  if (which.ord == "2") #reverse lists for Yates order
  { 
    fnames <- fnames[kfac:1]
    freps <- freps[kfac:1]
  }
  return(list(factors = fnames,reps = freps))
}

#Supply recipient in standard order
#Returns a permutation of the standard order of the recipient factors
"fac.recip.cross" <- 
  function(recipient, rec.names, rec.levels, nested.recipients=NULL, except=NULL, 
           allocated, seed=NULL)
  { 
    if (is.null(allocated))
    {
      n <- prod(rec.levels)
    } else
    {
      if (is.data.frame(allocated))
        n <- nrow(allocated)
      else
        n <- length(allocated)
    }
    which.unr <- 0
    if (is.data.frame(recipient)) #for data.frame
    { 
      which.unr <- 1
      facgen <- recipient
    }
    else 
    { 
      facgen <- fac.gen(generate=rec.names)
    }
    #generate random factor values
    nunr <- ncol(facgen)
    for(i in 1:nunr)
    { 
      if (!(names(rec.names)[i] %in% except))
        rno <- runif(rec.levels[i])[as.integer(facgen[[i]])]
      else
        rno <- as.integer(facgen[[i]])
      if (i == 1)
        facrecip <- data.frame(rno)
      else
        facrecip <- data.frame(facrecip,rno)
    }
    names(facrecip) <- names(rec.names)
    facrecip.ord <- 1:n
    facrecip.ord[do.call(order, facrecip)] <- facrecip.ord
    facrecip <- fac.divide(facrecip.ord, rec.names)
    return(facrecip)
  }

#Supply recipient in standard order
#Returns a permutation of the standard order of the recipient factors
"fac.recip.nest" <- function(recipient, rec.names, rec.levels, nested.recipients=NULL, 
                             except=NULL, allocated, seed=NULL)
{ 
  nnested <- length(nested.recipients)
  names.nested <- names(nested.recipients)
  if (!all(names.nested %in% names(rec.names)))
    stop("The following nested recipients are not in the set of recipient factors: ",
         paste(names.nested[!(names.nested %in% names(rec.names))], 
               collapse = ", "),"\n\n")
  #check nested factors for transitivity
  for (i in 1:nnested)
  { 
    if (!all(nested.recipients[[i]] %in% names(rec.names)))
      stop("The following nesting recipients are not in the set of recipient factors: ",
           paste(nested.recipients[[i]][!(nested.recipients[[i]] %in% names(rec.names))], 
                 collapse = ", "),"\n\n")
    nnest <- length(nested.recipients[[i]])
    
    for (j in 1:nnest)
    { 
      knested <- match(nested.recipients[[i]][j], names.nested)
      if (!is.na(knested))
      { 
        if (!all(nested.recipients[[knested]] %in% nested.recipients[[i]]))
        {  
          stop(names.nested[i]," is nested in ",nested.recipients[[i]][j],
               " but is not nested in all those that ",nested.recipients[[i]][j]," is.")
        }
      }
    }
  }
  #order the recipient factors so that any subsets of the generalized factor
  #formed from the nesting factors of a nested factor are to the left of the 
  #nested factor. Also sort numbers of levels in levels and form vector 
  #rec.denestord giving location of recipient factors in nestord.
  #first a list of all non-nested factors
  if (is.null(allocated))
  {
    n <- prod(rec.levels)
  } else
  {
    if (is.data.frame(allocated))
      n <- nrow(allocated)
    else
      n <- length(allocated)
  }
  which.unr <- 0
  if (is.data.frame(recipient)) #for data.frame
  { 
    which.unr <- 1
    nunr <- ncol(recipient)
  }
  else 
    nunr <- length(rec.names)
  kunr <- 0
  rec.levels.nestord <- rep(1, length=nunr)
  rec.denestord <- rep(NA, length=nunr)
  for(i in 1:nunr)
  { 
    if (is.na(match(names(rec.names)[i], names.nested)))
    { 
      kunr <- kunr+1
      rec.denestord[i] <- kunr
      rec.levels.nestord[kunr] <- rec.levels[i]
      if (kunr == 1)
      { 
        rec.nestord <- list(rec.names[[i]])
        names(rec.nestord) <- names(rec.names)[i]
      }
      else
      { 
        knames <- list(rec.names[[i]])
        names(knames) <- names(rec.names)[i]
        rec.nestord <- c(rec.nestord, knames)
      }
    }
  } 
  #now sort nested factors for number of nesting factors and add to list in this order
  if (!is.null(nested.recipients))
  { 
    nested.sort <- sort(sapply(nested.recipients, FUN=length))
    for (i in 1:nnested)
    { 
      kunr <- kunr+1
      krec.name <- names(nested.sort)[i]
      kno <- match(krec.name, names(rec.names))
      rec.denestord[kno] <- kunr
      rec.levels.nestord[kunr] <- rec.levels[kno]
      if (kunr == 1)
      { 
        rec.nestord <- list(rec.names[[kno]])
        names(rec.nestord) <- names(rec.names)[kno]
      }
      else
      { 
        knames <- list(rec.names[[kno]])
        names(knames) <- names(rec.names)[kno]
        rec.nestord <- c(rec.nestord, knames)
      }
    }
  }
  #generate random factor values with the factors in nested order
  facgen <- fac.gen(generate=rec.nestord)
  for(i in 1:nunr)
  { 
    if (names(rec.nestord)[i] %in% except)
      rno <- as.integer(facgen[[i]])
    else
    { 
      knest <- match(names(rec.nestord)[i], names.nested)
      if (is.na(knest))   #nonnested factor
      { 
        rno <- runif(rec.levels.nestord[i])[as.integer(facgen[[i]])]
      }
      else     #nested factor
      { 
        kfac <- length(nested.recipients[[knest]])+1
        kfacnos <- rep(1, length=kfac)
        kfacnos[1] <- i
        for (j in 1:(kfac-1))
        { 
          kfacnos[j+1] <- match(nested.recipients[[knest]][j], names(rec.nestord))
          if (is.na(kfacnos[j+1]))
            stop("Nesting factor not in list of recipient factors.")
        }
        sort(kfacnos)
        #determine number of random nos required and 
        #generate radix to expand to n-vector
        radix <- rep(1, length=n)
        each <- 1
        for (j in kfac:1)
        { 
          radix <- radix + (as.integer(facgen[[kfacnos[j]]])-1)*each
          each <- each*rec.levels.nestord[kfacnos[j]]
        }
        rno <- runif(each)[radix]
      }
    }
    if (i == 1)
      facrecip <- data.frame(rno)
    else
      facrecip <- data.frame(facrecip,rno)
  }
  names(facrecip) <- names(rec.nestord)
  facrecip.ord <- 1:n
  facrecip.ord[do.call(order, facrecip)] <- facrecip.ord
  facrecip <- fac.divide(facrecip.ord, rec.nestord)
  #reorder facrecip so that the factors and their values are in an appropriate 
  #order for the original factor order i.e. in rec.names
  facgen <- fac.gen(generate=rec.names)
  #form radix from facgen, numbers in radix labelling the levels of the factors
  #in order as for rec.nestord
  radix <- rep(1, length=n)
  each <- 1
  #loop over factors in rec.nestord
  for (j in nunr:1)
  { 
    kno <- match(names(rec.nestord)[j], names(rec.names))
    radix <- radix + (as.integer(facgen[[kno]])-1)*each
    each <- each*rec.levels[kno]
  }
  facrecip <- facrecip[radix,]
  facrecip <- facrecip[, rec.denestord]
  return(facrecip)
}

"designRandomize" <- function(allocated = NULL, recipient, nested.recipients=NULL, 
                              except=NULL, seed=NULL, unit.permutation = FALSE, ...)
{
  #generate a layout for a design consisting of allocated factors that are 
  #allocated to the recipient factors, taking into account the nesting between
  #the recipient factors.
  #recipient is a data.frame or list of factors (no numbers), along with their 
  #levels. If it is a list, each component of the list is a factor name that 
  #contains either a single numeric value that is the number of levels, a numeric 
  #vector that contains the levels of the factor or a character vector that 
  #contains the labels of the factor.
  #nested.factor is a list of the factors in recipient that are nested in other
  #factors in recipient. The factors within which they are nested are the 
  #elements of the component.
  #allocated is a factor or a data frame containing the generated factor or 
  #factors to be allocated.
  
  #Deal with deprecated recipient and randomized function arguments
  tempcall <- list(...)
  if (length(tempcall))
  {
    if ("unrandomized" %in% names(tempcall))
      stop("unrandomized has been deprecated in designRandomize - use recipient")
    if ("randomized" %in% names(tempcall))
      stop("randomized has been deprecated in designRandomize - use allocated")
  }
  
  #process allocated argument
  if (!is.null(allocated))
  {
    if (!is.data.frame(allocated) & !is.factor(allocated))
      stop("allocated must be a factor or data frame.")
    if (is.factor(allocated))
    {
      allocname <- deparse(substitute(allocated))
      allocated <- data.frame(allocated)
      names(allocated) <- allocname
    }
    n <- nrow(allocated)
  } else
    n <- 0
  #process seed argument
  if (!is.null(seed))
    set.seed(seed)
  #process recipient argument, form recipient factor list and, 
  #if a data.frame, ensure in standard order  
  if (is.data.frame(recipient)) #for data.frame
  { 
    if (!all(sapply(recipient, FUN=is.factor)))
      stop("All columns in the recipient data.frame must be factors")
    nunr <- ncol(recipient)
    if (!is.null(allocated) & nrow(recipient) != n)
      stop("The number of rows in the recipient data frame and the length of the factor(s) must be equal.")
    rec.names <- vector("list", length = nunr)
    names(rec.names) <- as.list(names(recipient))
    for (i in 1:nunr)
    {
      rec.names[[i]] <- levels(recipient[[i]])
      rec.names[[i]] <- rec.names[[i]][rec.names[[i]] %in% unique(recipient[[i]])]
      
    }
    #Order recipient into standard order
    perm.recip <- do.call(order, recipient)
    recipient.ord <- recipient[perm.recip,]
  }
  else
  { 
    if (!is.list(recipient))  #for (generate) list
      stop("recipient must be a list or a data.frame.")
    if (any(names(recipient) == ""))
      stop("all components of recipient list must be named.")
    facs.reps <- factor.list(recipient, order="standard")
    rec.names <- facs.reps$factors
    nunr <- length(rec.names)
    #recipient must be in standard order
    perm.recip <- 1:n
    recipient.ord <- recipient
  }
  #if except is not NULL, check that it contains only recipient factors
  if (!is.null(except))
  { 
    if (!is.character(except))
      stop("except must be a character vector")
    if (!all(except %in% names(rec.names)))
      stop("except must contain the names of only recipient factors")
  }
  #form vector of numbers of levels
  rec.levels <- rep(1, times=nunr)
  for (i in 1:nunr)
  { 
    if (is.numeric(rec.names[[i]]) | is.character(rec.names[[i]]))
    {
      if (is.numeric(rec.names[[i]]))
        rec.levels[i] <- rec.names[[i]]
      else
        rec.levels[i] <- length(rec.names[[i]])
    } else
      stop("Levels of factors must be specified using either numeric or character vectors")
  }
  if (is.null(allocated)) 
  {
    n <- prod(rec.levels)
  } else
  {
    if (n != prod(rec.levels))
      stop("The product of the numbers of levels of the recipient factors ", 
           "must equal the length of the allocated factors.")
  }
  #process nested.factor argument to get randomized recipient factors
  if (is.null(nested.recipients))
    facrecip <- fac.recip.cross(recipient=recipient.ord, rec.names=rec.names, 
                                rec.levels=rec.levels, nested.recipients=nested.recipients, 
                                except=except, seed=seed, allocated=allocated)
  else
  { 
    if (!is.list(nested.recipients))
      stop("nested.recipients must be a list.")
    names.nested <- names(nested.recipients)
    if (any(names.nested == ""))
      stop("all components of nested.recipients must be named.")
    if (!all(sapply(nested.recipients, FUN=is.character)) )
      stop("All elements of the nested.recipients list must be of class character")
    facrecip <- fac.recip.nest(recipient=recipient.ord, rec.names=rec.names, 
                               rec.levels=rec.levels, nested.recipients=nested.recipients, 
                               except=except, seed=seed, allocated=allocated)
  }
  
  #Form layout which is to be in standard or data frame order for the 
  #recipient factors supplied in recipient i.e. in rec.names order
  if (is.data.frame(recipient)) #get into recipient data.frame order
  { 
    for (i in 1:nunr)
      attributes(facrecip[[i]]) <- attributes(recipient[[i]])

    #perm.derand, on the right, puts facrecip order into standard order 
    # and, on the left,  puts standard order into facrecip order
    #perm.recip, on the right, puts recipient order into standard order 
    # and, on the left,  puts standard order into recipient order
    #Note: order(perm.xxxx) on the right is equivalent to perm.xxx on the left
    perm.derand <- do.call(order, facrecip)

    faclay <- facrecip
    #change faclay from standard to recipient order (to correspond with allocated)
    faclay[perm.recip, ] <- facrecip 
 
    # #order to get from standard order to permuted data frame (perm.dat on the right)
    # perm.dat <- vector("numeric", length=n)
    # perm.dat[perm.recip] <- perm.derand
    # #Order to unpermute the data (perm.derand.dat on the right)
    # perm.derand.dat <- vector("numeric", length=n)
    # perm.derand.dat[perm.recip] <- order(perm.dat)
    
    #order to get from recipient data frame to permuted data frame
    #(perm.dat on the right)
    perm.dat <- vector("numeric", length=n)
    perm.dat <- perm.recip
    perm.dat <- perm.dat[perm.derand]
    perm.dat[perm.recip] <- perm.dat
    #Order to unpermute the data (perm.derand.dat on the right)
    perm.derand.dat <- vector("numeric", length=n)
    perm.derand.dat <- order(perm.dat)

    #join faclay with allocated factors and unpermute so that recipient factors are in data.frame order
    if (is.null(allocated))
    {
      if (unit.permutation)
        faclay <- data.frame(.Units = 1:n, .Permutation = perm.derand.dat, faclay)
    }
    else
    {
      faclay <- data.frame(faclay, allocated)
      if (unit.permutation)
              faclay <- data.frame(.Units = 1:n, .Permutation = perm.derand.dat,
                             faclay[perm.dat, ])
      else
        faclay <- faclay[perm.dat, ]
    }
  }
  else  #get in standard order for rec.names
  { 
    perm.derand <- do.call(order, facrecip)
    if (is.null(allocated))
    {
      faclay <- facrecip
      if (unit.permutation)
        faclay <- data.frame(.Units = 1:n, .Permutation = perm.derand, 
                             faclay)
    }
    else
    {
      faclay <- data.frame(facrecip, allocated)
      if (unit.permutation)
        faclay <- data.frame(.Units = 1:n, .Permutation = order(perm.derand), 
                             faclay[perm.derand, ])
      else
        faclay <- faclay[perm.derand, ]
    }
  }
  rownames(faclay) <- 1:n
  return(faclay)
}

