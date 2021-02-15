#"fac.getinTerm" <- function(term)
#  #function to return the set of factors/variables in a term separated by ':"
#{ unlist(strsplit(term, ":", fixed=TRUE))}


as.data.frame.pstructure <- function(x, row.names = NULL, optional = FALSE, ..., 
                                     omit.marginality = FALSE)
{
  if (!inherits(x, what = "pstructure"))
    stop("Must supply an object of class pstructure")
  #Remove aliasing
  x <- x[-match("aliasing", names(x))]
  
  #Convert projectors to df
  if (inherits(x$Q[[1]], what = "projector"))
    x$Q <- unlist(lapply(x$Q, degfree))
  
  #Form data.frame
  if (omit.marginality)
    x <- as.data.frame.list(x[-4], stringsAsFactors = FALSE)
  else
  {
    if (x$Q[1] == "Mean")
      Q <- Q[-1]
    x <- as.data.frame.list(x, row.names = row.names, optional = optional, 
                            stringsAsFactors = FALSE)
  }
  names(x)[1] <- "df"
  return(x)
}

print.pstructure <- function(x, which = "all", ...)
{
  if (!inherits(x, "pstructure"))
    stop("Must supply an object of class pstructure")
  options <- c("projectors", "marginality", "aliasing", "all")
  print.opt <- options[unlist(lapply(which, check.arg.values, options=options))]
  if (print.opt %in% c("projectors", "all"))
  {
    cat("\nProjection matrix summary\n\n")
    print(as.data.frame(x, omit.marginality = TRUE))
  }
  if (print.opt %in% c("marginality", "all"))
  { 
    cat("\nMarginality matrix\n\n")
    print(x$marginality)
  }
  if (print.opt %in% c("aliasing", "all"))
  { 
    if (is.null(x$aliasing))
      cat("\n\nNo aliasing between sources in this pstructure object\n\n")
    else
    {
      #cat("\nTable of information (partially) aliased with previous sources derived from the same formula\n\n")
      print(x$aliasing)
    }
  }
  invisible()
}

"marginality.pstructure" <- function(object, ...)
  #function to extract the marginality matrices from a pcanon object
{ if (!inherits(object, "pstructure"))
  stop("Must supply an object of class pstructure as produced by designAnatomy")

  #Get marginality matrices
  if (is.null(object$marginality))
    marginality <- NULL
  else
  {
    marginality <- object$marginality
  }
  return(marginality)
}  

print.aliasing <- function(x, which.criteria = c("aefficiency","eefficiency","order"), 
                           ...)
{
  if (!is.null(x))
  {
    if (!inherits(x, "aliasing"))
      stop("Must supply an object of class aliasing")
    
    if (nrow(x) > 0)
    {
      #check which.criteria arguments
      criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                    "order", "dforthog")
      options <- c(criteria, "none", "all")
      kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                         options=options))]
      if ("all" %in% kcriteria)
        kcriteria <- criteria
      anycriteria <- !("none" %in% kcriteria)
      
      #Print aliasing
      cols <- c("Source","df")
      if (!all(is.na(x$Alias)))
        cols <- c(cols, "Alias")
      if ("In" %in% names(x))
        cols <- c(cols, "In")
      if (anycriteria)
      {
        cols <- c(cols, kcriteria)
        if (!all(kcriteria %in% names(x)))
          stop("Not all requested efficiency criteria are available in aliasing for printing")
      }
      cat(attr(x, which = "title"))
      y <- x[cols]
      for (kcrit in kcriteria)
      {
        if (kcrit == "order")
        {
          y[kcrit] <- formatC(y[[kcrit]], format="f", digits=0, width=5)
        } else
        {
          if (kcrit == "dforthog")
          {
            y[kcrit] <- formatC(y[[kcrit]], format="f", digits=0, width=8)
          } else
          {
            y[kcrit] <- formatC(y[[kcrit]], format="f", digits=4, width=11)
          }
        }
        y[kcrit] <- gsub("NA", "  ", y[[kcrit]])
      }
      
      print.data.frame(y, na.print="  ", right=FALSE, row.names=FALSE)
    }
  }
  invisible()
}

#Function to convert term names into source names
#term.names is an unnamed character vector with the nterms names of the terms to be converted
#marginality is a nterms x nterms matrices of 0 and 1s with the column entries specifying 
#which of the row terms is marginal to the column term.
formSources <- function(term.names, marginality, grandMean = FALSE)
{
  if (grandMean)
  {
    gmterm <- term.names[1]
    term.names <- term.names[-1]
  }
  nterms <- length(term.names)
  names(term.names) <- term.names
  if (nterms == 1)
    sources <- term.names
  else
  {
    sources <- vector(mode = "character", length = nterms)
    names(sources) <- term.names
    fac.list <- lapply(term.names, fac.getinTerm)
    names(fac.list) <- term.names
    #Find all the factors
    facs <- unique(unlist(fac.list))
    #Do an incidence matrix of factors in terms
    fac.incidence <- do.call(rbind, lapply(fac.list, 
                                           function(kfac.list, facs)
                                           {
                                             ind <- rep(FALSE, length(facs))
                                             names(ind) <- facs
                                             ind[kfac.list] <- TRUE
                                             return(ind)
                                           }, facs = facs))
    #Protect generalized or compound or joint factors by changing ":" to "."
    #Find factors with identical incidence patterns
    nfac <- length(facs)
    k1 <- 0
    newterms <- list()
    for (k in 1:(nfac-1))
    {
      k1 <- k1 + 1
      facs.identical <- unlist(lapply(facs[(k1+1):nfac], 
                                      function(jfac, kfac, fac.incidence)
                                      {
                                        is.same <- all(fac.incidence[,jfac] == fac.incidence[,kfac])
                                        return(is.same)
                                      }, kfac = facs[k1], fac.incidence = fac.incidence))
      #Have a potential compound factor?
      if (any(facs.identical))
      {
        facs.identical <- c(facs[k1],
                            c(facs[(k1+1):nfac])[facs.identical])
        nposns <- length(facs.identical)
        #Change all the terms with these identical factors to have a "."
        newterms <- lapply(term.names[fac.incidence[,facs[k1]]],
                           function(term, fac.list)
                           {
                             newterm <- term
                             kfac.list <- fac.list[[term]]
                             fac.posns <- sort(match(facs.identical, kfac.list))
                             if (!any(is.na(fac.posns)))
                             {
                               nposns <- length(fac.posns)
                               #Can only be a compound factor if equal factors are consecutive
                               if ((nposns-1) == fac.posns[nposns] - fac.posns[1])
                               {
                                 if (fac.posns[1] == 1 & nposns == length(kfac.list))
                                   newterm <- paste("(",
                                                    paste(kfac.list[fac.posns[1]:fac.posns[nposns]], 
                                                          collapse = "."), 
                                                    ")", sep = "")
                                 else
                                 {
                                   newterm <- paste("(",
                                                    paste(kfac.list[fac.posns[1]:fac.posns[nposns]], 
                                                          collapse = "."), 
                                                    ")", sep = "")
                                   if (fac.posns[1] != 1)
                                     newterm <- paste(c(kfac.list[1:(fac.posns[1] - 1)], newterm), 
                                                      collapse = ":")
                                   if (fac.posns[nposns] < length(kfac.list))
                                     newterm <- paste(c(newterm, 
                                                        kfac.list[(fac.posns[nposns]+1):length(kfac.list)]), 
                                                      collapse = ":")
                                 }
                                 
                               }
                             }
                             newterm <- list(term = term, newterm = newterm)
                             return(newterm)
                           }, fac.list = fac.list)
        #Copy new names to term.names and reform fac.list
        if (length(newterms) > 0)
        {
          term.names <- unlist(lapply(term.names, 
                                      function(term, newterms)
                                      {
                                        for (k in 1:length(newterms))
                                        {
                                          if (newterms[[k]]$term == term)
                                            term <- newterms[[k]]$newterm
                                        }
                                        return(term)
                                      }, newterms = newterms))
          fac.list <- lapply(term.names, fac.getinTerm)
          names(fac.list) <- term.names
          colnames(marginality) <- unlist(term.names)
          rownames(marginality) <- unlist(term.names)
        }
        k1 <- k1 + nposns - 1
      }
      if (k1 >= (nfac-1))
        break
    }
    
    ###Recompute fac.incidence
    facs <- unique(unlist(fac.list))
    fac.incidence <- do.call(rbind, lapply(fac.list, 
                                           function(kfac.list, facs)
                                           {
                                             ind <- rep(FALSE, length(facs))
                                             names(ind) <- facs
                                             ind[kfac.list] <- TRUE
                                             return(ind)
                                           }, facs = facs))
    #Is there any implicit nesting
    implicit <- vector(mode = "list", length = length(term.names))
    names(implicit) <- term.names
    for (term in term.names)
    {
      k <- match(term, rownames(fac.incidence))
      kmarg <- as.logical(marginality[, k])
      kmarg[k] <- FALSE
      if (any(kmarg))
      {
        kmarg <- rownames(fac.incidence)[kmarg]
        implicit[[term]] <- !unlist(lapply(kmarg, 
                                           function(kterm, term, fac.list)
                                           {
                                             is.subset <- (length(intersect(fac.list[[term]], 
                                                                            fac.list[[kterm]])) > 0)
                                             return(is.subset)
                                           }, term = term, fac.list = fac.list))
        if (any(implicit[[term]]))
          for (kimpl in kmarg[implicit[[term]]])
            fac.incidence[term, ] <- fac.incidence[term, ] | fac.incidence[kimpl, ]
      }
    }
    
    ###Loop over terms
    marg.compliant <- marginality == 1
    kmarg.terms <- vector(mode = "list", length = nterms)
    names(kmarg.terms) <- term.names
    for (term in term.names)
    {
      #reduce terms to the marginality-compliant set
      j <- k <- match(term, term.names)
      #      marg.compliant <- as.logical(marginality[,k])
      if (j <= 1)
      {
        kmarg.terms[[term]] <- character()
      } else
      {
        if (!any(marg.compliant[1:(j-1), k]))
          kmarg.terms[[term]] <- character()
        else
        {
          repeat
          {
            if (j <= 0)
              break
            #get highest (row number) marginal term (j) left below k
            if (j == 1 | !any(marg.compliant[1:(j-1), k]))
              j <- numeric()
            else
              j <- max(c(1:(j-1))[marg.compliant[1:(j-1), k]])
            #repeat until no further ones below current position
            if (length(j) == 0 || is.infinite(j))
              break
            #remove any terms marginal to j from list of terms marginal to k
            if (j > 1 & any(marginality[1:(j-1), j] == 1))
            {
              marg.compliant[1:(j-1), k][marginality[1:(j-1), j]  == 1] <- FALSE
            }
          }
          kmarg.terms[[term]] <- term.names[marg.compliant[,k]]
          kmarg.terms[[term]] <- kmarg.terms[[term]][-match(term, kmarg.terms[[term]])]
        }
      }
    }
    for (term in term.names)
    {
      k <- match(term, term.names)
      #Use marginality and the marginality-compliant set to determine the source name
      if (length(kmarg.terms[[term]]) == 0) #no nesting terms
      {
        sources[k] <- gsub(".", ":", term, fixed = TRUE)
        nch <- nchar(sources[k])
        if (substr(sources[k], 1, 1) == "(" & substr(sources[k], nch, nch) == ")")
          sources[k] <- substr(sources[k], 2, nch-1)
      } else
      {
        if (length(kmarg.terms[[term]]) == 1) #kmarg.term must be must be nesting
        {
          facs.nesting <- fac.list[kmarg.terms[[term]]][[1]]
          kfacs <- colnames(fac.incidence)[fac.incidence[kmarg.terms[[term]],]]#do we need this to get implicit?
          #re-order the factors in fac.nesting
          facs.nesting <- c(na.omit(setdiff(kfacs, facs.nesting)),
                            facs.nesting)
          #get the factors in the term other than nesting
          facs.nested <- setdiff(fac.list[term][[1]], facs.nesting)
        } else #more than 1
        {
          facs.nesting <- NULL
          facs.nested <- NULL
          #find the supremum for the term, which gives the nesting factors
          common.marginals <- TRUE
          for (c in match(kmarg.terms[[term]], term.names))
            common.marginals <- common.marginals & marginality[,c]  == 1 #marg.compliant[,c]
          common.marginals <- term.names[common.marginals]
          facs.nesting <- character()
          for (cterm in common.marginals)
          {
            facs.nesting <- c(facs.nesting, fac.list[cterm][[1]])
          }
          facs.nesting <- unique(facs.nesting)
          #Check factor names for any extra suprema
          kfacs <- unlist(lapply(kmarg.terms[[term]],
                                 function(kmarg, fac.incidence)
                                 {
                                   fac <- colnames(fac.incidence)[fac.incidence[kmarg,]]
                                 }, fac.incidence = fac.incidence))
          kfacs <- c(na.omit(setdiff(kfacs, fac.list[term][[1]])),
                     fac.list[term][[1]])
          kfacs <- kfacs[!(kfacs %in% facs.nesting)]
          if (length(kfacs) > 0)
          {
            for (kfac in kfacs)
            {
              #Is this factor in all marginality-compliant terms
              if (all(fac.incidence[kmarg.terms[[term]], kfac]))
                facs.nesting <- c(facs.nesting, kfac)
            } 
          }
          #nested factors are those in the term other than the nesting factornesting
          facs.nested <- setdiff(fac.list[term][[1]], facs.nesting)
        }
        
        #Paste the factors and punctuation
        sources[k] <- paste(facs.nested, collapse = "#")
        nch <- nchar(sources[k])
        if (substr(sources[k], 1, 1) == "(" & substr(sources[k], nch, nch) == ")")
          sources[k] <- substr(sources[k], 2, nch-1)
        if (length(facs.nesting) > 0)
        {
          facs.nesting <- paste(facs.nesting, collapse = ":")
          facs.nesting <- gsub("(", "", facs.nesting, fixed = TRUE)
          facs.nesting <- gsub(")", "", facs.nesting, fixed = TRUE)
          sources[k] <- paste(sources[k], "[", facs.nesting,"]", sep = "")
        }
        sources[k] <- gsub(".", ":", sources[k], fixed = TRUE)
      }
    }
  }
  if (grandMean)
  {
    sources <- c(gmterm, sources)
    names(sources)[1] <- gmterm    
  }
  return(sources)
}

"projs.jandw" <- function(R, Q, which.criteria = c("aefficiency","eefficiency","order"), 
                          aliasing.print = TRUE)
  #A function to use J&W to orthogonalise the set of Q to R and previous Q
{ 
  if (!is.list(Q))
    stop("The matrices to orthogonalize to must be in a list")
  
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  #set up aliasing summary
  nc <- 3 + length(criteria)
  aliasing <- data.frame(matrix(nrow = 0, ncol=nc))
  colnames(aliasing) <- c("Source", "df", "Alias", criteria)
  class(aliasing) <- c("aliasing", "data.frame")
  eff.crit <- vector(mode="list",length=length(criteria))
  names(eff.crit) <- criteria

  #do orthogonalization
  terms <- names(Q)
  aliased <- rep(FALSE, length(terms))
  names(aliased) <- terms
  for (i in 1:length(terms))
  { 
    decomp <- proj2.combine(R, Q[[terms[i]]])
#    Q[[terms[i]]] <- decomp$Qconf
#    R <- decomp$Qres
    df <- degfree(decomp$Qconf)
    if (df == 0)
    { 
      warning(paste(terms[[i]],"is aliased with previous terms in the formula",
                    "and has been removed", sep=" "))
      eff.crit[criteria] <- 1 #0
      aliasing <- rbind(aliasing, 
                        data.frame(c(list(Source = terms[[i]], df = degfree(Q[[terms[i]]]), 
                                          Alias = "unknown"), 
                                     eff.crit), 
                                   stringsAsFactors = FALSE))
    } else
    { 
      keff.crit <- efficiency.criteria(decomp$efficiencies)
      if ((df - keff.crit[["dforthog"]]) != 0)
      { 
        P <- projector(mat.I(nrow(R)) - R)
        decompP <- proj2.combine(P, Q[[terms[i]]]) 
        keff.crit <- efficiency.criteria(decompP$efficiencies)
        warning(paste(terms[[i]],"is partially aliased with previous terms in the formula", sep=" "))
        aliasing <- rbind(aliasing, 
                          data.frame(c(list(Source = terms[[i]], df = degfree(decompP$Qconf), Alias = "unknown"), 
                                       keff.crit), 
                                     stringsAsFactors = FALSE))
      }
    }
    Q[[terms[i]]] <- decomp$Qconf
    R <- decomp$Qres
  }
  
  #Print out the efficiency criteria if any aliasing
  if (nrow(aliasing) > 0)
  { 
    rownames(aliasing) <- NULL
    class(aliasing) <- c("aliasing", "data.frame")
    attr(aliasing, which = "title") <- 
      "\nTable of information (partially) aliased with previous sources derived from the same formula\n\n"
    if (aliasing.print)
    {
      if (anycriteria)
        print(aliasing, which = kcriteria)
      else
        print(aliasing, which = "none")
    }
  } else
  {
    aliasing <- NULL
  }
  Q <- Q[!aliased]
  return(list(Q = Q, aliasing = aliasing))
}

"pstructure.formula" <- function(formula, keep.order = TRUE, grandMean = FALSE, 
                                 orthogonalize = "hybrid", labels = "sources", 
                                 marginality = NULL, check.marginality = TRUE, 
                                 omit.projectors = FALSE, 
                                 which.criteria = c("aefficiency","eefficiency","order"), 
                                 aliasing.print = TRUE, data = NULL, ...)
{ #generate a set of mutually orthogonal projection matrices, one for each term in the formula
  options <- c("differencing", "eigenmethods", "hybrid")
  opt <- options[check.arg.values(orthogonalize, options)]
  options <- c("terms", "sources")
  lab.opt <- options[check.arg.values(labels, options)]
  if (orthogonalize == "eigenmethods" & lab.opt == "sources" & is.null(marginality))
    warning("When orthogonalize = eigenmethods and marginality not supplied, labels = sources is ignored")
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  #Initialize
  if (is.null(data) | !is.data.frame(data))
    stop("Must supply a data.frame for data")
  n <- nrow(data)
  Q.G = projector(matrix(1, nrow=n, ncol=n)/n)
  aliasing <- NULL
  
  #check for ~1 and create pstructure
  if (formula == ~1)
  {
    if (!grandMean) 
      warning("formula is for grand mean, yet grandMean is FALSE; it will be saved")
    terms <- "Mean"
    nterms <- 1
    Q <- list(Mean = Q.G)
    marginality <- diag(1, nrow = nterms, ncol = nterms)
    rownames(marginality) <- terms
    colnames(marginality) <- terms
    sources <- terms
  } else #formula is not ~1
  {
    #get terms, check all factors/variates in data
    terms <- attr(terms(formula, keep.order = keep.order, ...), which="term.labels")
    nterms <- length(terms)
    fac.list <- lapply(terms, fac.getinTerm)
    names(fac.list) <- terms
    #Find all the factors
    facs <- unique(unlist(fac.list))
    if (!all(facs %in% names(data)))
      stop(paste("Some factor/covariates missing from data:", 
                 paste0(facs[!(facs %in% names(data))], collapse = ", ")))
    
    #form mean operators for each term
    fac.modl <- model.frame(formula, data=data)
    if (grandMean)  #add grand mean term if grandMean is TRUE
    {
      Q <- vector("list", length=nterms+1)
      names(Q) <- c("Mean", terms)
      Q[["Mean"]] <- Q.G
    } else
    {
      Q <- vector("list", length=nterms)
      names(Q) <- terms
    }
    for (k in 1:nterms)
    { 
      Q[[terms[k]]] <- model.matrix(as.formula(paste("~ ",terms[k])), data=fac.modl)
      Q[[terms[k]]] <- Q[[terms[k]]] %*% ginv(t(Q[[terms[k]]]) %*% Q[[terms[k]]]) %*% t(Q[[terms[k]]])
      Q[[terms[k]]] <- projector(Q[[terms[k]]] - Q.G)
    }
    
    #check marginality matrix when it is supplied
    if (!is.null(marginality))
    {
      if (nrow(marginality) != nterms | ncol(marginality) != nterms)
        stop(paste("The number of rows and columns in the supplied marginality matrix must be ",
                   "the same as the number of terms in the supplied formula"))
      if (!all(rownames(marginality) == terms) | !all(colnames(marginality) == terms))
        warning("Not all row and column names for the supplied marginality are the same as the expanded set of term names")
      if (!all(marginality == 1 | marginality ==  0))
        stop("All elements of the supplied marginality matrix must be 0 or 1")
    }
    
    #initialize marginality matrix
    marg.mat <- diag(1, nrow = nterms, ncol = nterms)
    rownames(marg.mat) <- terms
    colnames(marg.mat) <- terms
    
    #form projection matrices of the structure
    if (length(terms) == 1)
    {
      marginality <- as.matrix(c(1))
      rownames(marginality) <- colnames(marginality) <- terms
      sources <- terms <- names(Q)
      aliasing <- NULL
    } else
    { 
      if (orthogonalize == "hybrid") #by difference and eigenmethods
      { 
        #set up aliasing summary
        nc <- 3 + length(criteria)
        aliasing <- data.frame(matrix(nrow = 0, ncol=nc))
        colnames(aliasing) <- c("Source", "df", "Alias", criteria)
        eff.crit <- vector(mode="list",length=length(criteria))
        names(eff.crit) <- criteria
        
        #Loop over terms
        for (i in 2:length(terms))
        { 
          Q.work <- Q[[terms[i]]]
          aliasstatus <- "none"
          #loop over terms to see if any are marginal to term i
          for (j in 1:(i-1))
          { 
            if (marg.mat[i, i] == 1)
            {
              Qji <- Q[[terms[j]]] %*% Q.work
              if (!is.allzero(Qji)) # if not orthogonal then have to orthogonalize
              {
                if (is.allzero(Qji-Q[[terms[j]]])) #test marginality
                {
                  if (is.allzero(Q[[terms[j]]] - Q.work)) #are terms equal?
                  {
                    aliasstatus <- "full"
                    marg.mat[terms[i], terms[i]] <- 0
                    warning(paste(terms[[i]],"is aliased with previous terms in the formula",
                                  "and has been removed", sep=" "))
                    eff.crit[criteria] <- 1 #0
                    aliasing <- rbind(aliasing, 
                                      data.frame(c(list(Source = terms[[i]], 
                                                        df = degfree(Q[[terms[i]]]), #0, 
                                                        Alias = terms[[j]]), 
                                                   eff.crit), 
                                                 stringsAsFactors = FALSE))
                  } else # i is a marginal term to j - difference
                  {
                    marg.mat[terms[j], terms[i]] <- 1
                    Q.work <- projector(Q.work - Q[[terms[j]]])
                  }
                } else
                {
                  if (is.allzero(Qji - Q.work)) #i is aliased with j
                  {
                    aliasstatus <- "full"
                    marg.mat[terms[i], terms[i]] <- 0
                    warning(paste(terms[[i]],"is aliased with previous terms in the formula",
                                  "and has been removed", sep=" "))
                    eff.crit[criteria] <- 1 #0
                    aliasing <- rbind(aliasing, 
                                      data.frame(c(list(Source = terms[[i]], 
                                                        df = degfree(Q.work), #0, 
                                                        Alias = terms[[j]]), 
                                                   eff.crit), 
                                                 stringsAsFactors = FALSE))
                  } else #partial aliasing of j with i - orthogonalize
                  {
                    aliasstatus <- "partial"
                    decompP <- proj2.combine(Q.work, Q[[terms[j]]])
                    keff.crit <- efficiency.criteria(decompP$efficiencies)
                    aliasing <- rbind(aliasing, 
                                      data.frame(c(list(Source = terms[[i]], 
                                                        df = degfree(decompP$Qconf),
                                                        Alias = terms[[j]]), 
                                                   keff.crit[criteria]), 
                                                 stringsAsFactors = FALSE))
                    R <- projector(diag(1, nrow = n, ncol = n) - Q[[terms[j]]])
                    decomp <- proj2.combine(R, Q.work)
                    Q.work <- decomp$Qconf
                  }
                }
              }
            }
          }
          
          #Check that this term is orthogonal to previous projectors - necessary?
          # if (degfree(Q.work) != 0) 
          # { 
          #   i1 <- i - 1 
          #   if (i1 > 0)
          #     for (j in 1:i1)
          #       if (!is.allzero(Q.work %*% Q[[terms[j]]]))
          #       { 
          #         warning(paste("** Projection matrices for ",terms[i], " and ", terms[j], 
          #                       " are not orthogonal", sep=""))
          #       }
          # }
          if (aliasstatus == "partial")
          {
            decompR <- proj2.combine(Q[[terms[i]]], Q.work)
            keff.crit <- efficiency.criteria(decompR$efficiencies)
            aliasing <- rbind(aliasing, 
                              data.frame(c(list(Source = terms[[i]], 
                                                df = degfree(decompR$Qconf),
                                                Alias = "## Information remaining"), 
                                           keff.crit[criteria]), 
                                         stringsAsFactors = FALSE))
          }
          if (aliasstatus == "full")
          {
            keff.crit <- rep(0, length(criteria))
            names(keff.crit) <- criteria
            aliasing <- rbind(aliasing, 
                              data.frame(c(list(Source = terms[[i]], 
                                                df = degfree(Q[[terms[i]]]),
                                                Alias = "## Aliased"), 
                                           keff.crit[criteria]), 
                                         stringsAsFactors = FALSE))
          }
          Q[[terms[i]]] <- Q.work
        }
        
        #Print out the aliasing summary if any aliasing
        if (nrow(aliasing) > 0)
        { 
          rownames(aliasing) <- NULL
          class(aliasing) <- c("aliasing", "data.frame")
          #Change terms in aliasing to sources, if appropriate 
          if (lab.opt == "sources")
          {
            tmp <- marg.mat
            diag(tmp) <- 1
            sources <- formSources(terms, tmp, grandMean = grandMean)
            which.sources <- aliasing$Source %in% names(sources)
            if (any(which.sources))
              aliasing$Source[which.sources] <- sources[aliasing$Source[which.sources]]
            which.sources <- aliasing$Alias %in% names(sources)
            if (any(which.sources))
              aliasing$Alias <- sources[aliasing$Alias]
          }
          attr(aliasing, which = "title") <- 
            "\nTable of information (partially) aliased with previous sources derived from the same formula\n\n"
          
          #Print aliasing
          if (aliasing.print)
          {
            if (anycriteria)
              print(aliasing, which = kcriteria)
            else
              print(aliasing, which = "none")
          }
        } else
        {
          aliasing <- NULL
        }
        
        #Remove aliased terms
        aliased <- (diag(marg.mat) == 0)
        if (any(aliased))
        {
          marg.mat <- marg.mat[!aliased, !aliased]
          Q <- Q[!aliased]
        }
        terms <- names(Q)
        if (is.null(marginality))
        {
          sources <- formSources(terms, marg.mat, grandMean = grandMean)
          marginality <- marg.mat
        }
        else
        {
          if (check.marginality & all(marg.mat == marginality))
            warning("Supplied marginality matrix differs from that computed by pstructure formula")
          sources <- formSources(terms, marginality, grandMean = grandMean)
        }
        if (lab.opt == "sources")
          names(Q) <- sources
      } else
      {
        if (orthogonalize == "differencing") #by difference
        { 
          #set up aliasing summary
          nc <- 3
          aliasing <- data.frame(matrix(nrow = 0, ncol=nc))
          colnames(aliasing) <- c("Source", "df", "Alias")
          class(aliasing) <- c("aliasing", "data.frame")
          eff.crit <- vector(mode="list",length=length(criteria))
          names(eff.crit) <- criteria
          
          orthogonal <- TRUE
          fac.mat <- attr(terms(formula, keep.order = keep.order, ...), which="factors")
          for (i in 2:length(terms))
          { 
            Q.work <- Q[[terms[i]]]
            for (j in 1:length(terms))
            { 
              if (i != j)
              { 
                if (all((fac.mat[,j] != 0) == (fac.mat[,j] & fac.mat[,i])))
                {
                  Q.work <- Q.work - Q[[terms[j]]]
                  marg.mat[terms[j], terms[i]] <- 1
                }
              }
            }
            Q.work <- projector(Q.work)
            if (degfree(Q.work) == 0)
            {
              marg.mat[terms[i], terms[i]] <- 0
              warning(paste(terms[[i]],"is aliased with previous terms in the formula",
                            "and has been removed", sep=" "))
              eff.crit[criteria] <- 1 #0
              aliasing <- rbind(aliasing, 
                                data.frame(list(Source = terms[[i]], df = 0, 
                                                Alias = "unknown"), 
                                           eff.crit,
                                           stringsAsFactors = FALSE))
            }
            else #Check that this term is orthogonal to previous projectors
            { 
              i1 <- i - 1 
              if (i1 > 0)
                for (j in 1:i1)
                  if (!is.allzero(Q.work %*% Q[[terms[j]]]))
                  { 
                    warning(paste("** Projection matrices for ",terms[i], " and ", terms[j], 
                                  ' are not orthogonal\n',
                                  '   Try orthogonalize = "hybrid" or "eigenmethods"', sep=""))
                    eff.crit[criteria] <- NA
                    aliasing <- rbind(aliasing, 
                                      data.frame(list(Source = terms[[i]], df = NA, 
                                                      Alias = terms[[j]]),
                                                 eff.crit, 
                                                 stringsAsFactors = FALSE))
                  }
            }
            Q[[terms[i]]] <- Q.work
          }
          
          #Print out the aliasing summary if any aliasing
          if (nrow(aliasing) > 0)
          { 
            rownames(aliasing) <- NULL
            class(aliasing) <- c("aliasing", "data.frame")
            #Change terms in aliasing to sources, if appropriate 
            if (lab.opt == "sources")
            {
              tmp <- marg.mat
              diag(tmp) <- 1
              sources <- formSources(terms, tmp, grandMean = grandMean)
              which.sources <- aliasing$Source %in% names(sources)
              if (any(which.sources))
                aliasing$Source[which.sources] <- sources[aliasing$Source[which.sources]]
              which.sources <- aliasing$Alias %in% names(sources)
              if (any(which.sources))
                aliasing$Alias <- sources[aliasing$Alias]
            }
            attr(aliasing, which = "title") <- 
              "\nTable of information (partially) aliased with previous sources derived from the same formula\n\n"
            
            #Print aliasing
            if (aliasing.print)
              print(aliasing, which = "none")
          } else
          {
            aliasing <- NULL
          }
          
          #Remove aliased terms
          aliased <- (diag(marg.mat) == 0)
          if (any(aliased))
          {
            marg.mat <- marg.mat[!aliased, !aliased]
            Q <- Q[!aliased]
          }
          terms <- names(Q)
          if (is.null(marginality))
          {
            sources <- formSources(terms, marg.mat, grandMean = grandMean)
            marginality <- marg.mat
          }
          else
          {
            if (check.marginality & all(marg.mat == marginality))
              warning("Supplied marginality matrix differs from that computed by pstructure formula")
            sources <- formSources(terms, marginality, grandMean = grandMean)
          }
          if (lab.opt == "sources")
            names(Q) <- sources
        } else  #by recursive orthogonalization
        { 
          R <- projector(diag(1, nrow = n, ncol = n) - Q[[terms[1]]] - Q.G)
          if (length(terms) > 1)
          {
            Qorth <- projs.jandw(R, Q[terms[2:length(terms)]], which.criteria = kcriteria,
                                 aliasing.print = aliasing.print)
            Q <- c(Q[1:match(terms[1], names(Q))], Qorth$Q)
            aliasing <- Qorth$aliasing
          }
          terms <- names(Q)
          if (is.null(marginality))
          {
            sources <- terms
            lab.opt <- "terms"
          }
          else
          {
            aliased <- !(terms %in% names(Q))
            if (grandMean)
              aliased <- aliased[-1]
            marginality <- marginality[!aliased, !aliased]
            sources <- formSources(terms, marginality, grandMean = grandMean)
            names(Q) <- sources
            lab.opt <- "sources"
          }
        }
      }
    }
  }

  #Create pstructure.object and return it
  struct <- list(Q = Q, terms = terms, sources = sources, 
                 marginality = marginality, aliasing = aliasing)
  class(struct) <- c("pstructure", "list")
  attr(struct, which = "labels") <- lab.opt
  if (omit.projectors)
    struct$Q <- unlist(lapply(struct$Q, degfree))

  return(struct)
}    

