"combined.split" <- function(x)
{ 
  terms <- strsplit(x, "[& ]")[[1]]
  return(terms)
}

"findDecomposition" <- function(icomb, term, source, object, ncomb)
#Looks for the combination of term confounded with source in object 
#that consists of ncomb decompositions. Starts with decompistion icomb
{ 
  while (!(term %in% names(object[[icomb]][[source]])))
  { 
    icomb <- icomb + 1
    if (icomb > ncomb)
      stop(paste("Could not find",term,"when confounded with",source," "))
  }
  return(icomb)
}

"findLastwithSource" <- function(icomb, source, object, ncomb)
  #Looks for the last combination that contains source
{ 
  while (!(source %in% names(object[[icomb]])))
  { 
    icomb <- icomb + 1
    if (icomb > ncomb)
      stop(paste("Could not find",source," "))
  }
  #Check for later sources
  repeat
  { 
    if (!(source %in% names(object[[icomb+1]])))
      break
    icomb <- icomb + 1
  }
  return(icomb)
}

"makeCombinedSource" <- function(terms)
{ 
  source <- paste(terms, collapse="&")
  return(source)     
}

"replace.proj" <- function(object)
  #Replace projectors with their degrees of freedom
{ 
#  if (!inherits(object, "p2canon"))
#    stop("object must be of class p2canon as produced by proj.2canon")
  
  #Loop through p2canon object
  Q1labels <- names(object)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { 
    object[[i]][["Q1res"]] <- degfree(object[[i]][["Q1res"]])
    Q2labels <- names(object[[i]])[-1]
    if (length(Q2labels) > 0)
    { 
      for (j in Q2labels)
        object[[i]][[j]][["Qproj"]] <- degfree(object[[i]][[j]][["Qproj"]])
    }
  }
  return(object)
}

# #a synonym for designAnatomy.data.frame, except that the data argument is in 
# #different positions and the default for labels is terms
# "projs.canon" <- function(formulae, orthogonalize = "hybrid", labels = "terms", 
#                           keep.order = TRUE, meanTerm = FALSE, 
#                           which.criteria = c("aefficiency", "eefficiency", "order"), 
#                           omit.projectors = c("pcanon", "combined"), data = NULL, ...)
# { 
#   warning("projs.canon may be deprecated in future versions, its synonym designAnatomy being preferred")
#   CombinedSets <- designAnatomy(formulae = formulae, data = data, 
#                                 orthogonalize = orthogonalize, 
#                                 labels = labels, 
#                                 keep.order = keep.order, 
#                                 grandMean = meanTerm, 
#                                 which.criteria = which.criteria, 
#                                 omit.projectors = omit.projectors, ...)
#   return(CombinedSets)
# }

"designAnatomy" <- function(formulae, data, keep.order = TRUE, grandMean = FALSE, 
                            orthogonalize = "hybrid", labels = "sources", 
                            marginality = NULL, check.marginality = TRUE, 
                            which.criteria = c("aefficiency", "eefficiency", 
                                               "order"), 
                            aliasing.print = FALSE, 
                            omit.projectors = c("pcanon", "combined"), ...)
{ #examine the relationships between the sets of mutually orthogonal projection matrices for the supplied formulae
  
  #Check arguments and intitialize
  if (!is.list(formulae))
    stop("formulae must be a list")
  form.named <- !is.null(names(formulae))
  if (form.named)
    if (any(is.na(names(formulae))))
      stop("If name one formula then must name all")
  ntiers <- length(formulae)
  ncomb <- ntiers - 1 
  if (!all(unlist(lapply(formulae, inherits, what="formula"))))
    stop("formulae must contain a list of formulae")

  #Check marginality argument
  if (!is.null(marginality))
  {
    if (!is.list(marginality))
      stop("marginality must be a list")
    if (!all(unlist(lapply(marginality, inherits, what="matrix")) || 
             unlist(lapply(marginality, is.null))))
      stop("marginality must contain a list of matrices or NULL components")
    #deal with the case that not all marginality matrices are supplied 
    if (length(marginality) != ntiers)
      if (is.null(names(marginality)) | is.null(names(formulae)))
      {
        stop(paste("if the marginality list is not the same length as the formulae list, ",
                   "these two lists must be named", sep = ""))
      } else #construct full length marginality list
      {
        tmp <- vector(mode = "list", length = ntiers)
        names(tmp) <- names(formulae)
        for (f in names(formulae))
          tmp[f] <- marginality[f]
        marginality <- tmp
      }
  }

  if (!is.data.frame(data))
    stop("data must be a data.frame ")
  n <- nrow(data)
  
  options <- c("hybrid", "differencing", "eigenmethods")
  if (length(orthogonalize) == 1)
    orthogonalize <- rep(orthogonalize, ntiers)
  else
    if (length(orthogonalize) != ntiers)
    { warning("Length of orthogonalize is not equal to 1 or the number of formulae - only using first value")
      orthogonalize <- rep(orthogonalize[1], ntiers)
    }
  orthog <- options[unlist(lapply(orthogonalize, check.arg.values, options=options))]
  
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  options <- c("pcanon", "combined", "none")
  omit.proj <- options[unlist(lapply(omit.projectors, check.arg.values, 
                                     options=options))]
  if ("none" %in% omit.proj)
    omit.proj <- "none"
  
  if (is.null(data) | !is.data.frame(data))
    stop("Must supply a data.frame for data")
  struct <- vector(mode = "list", length=ntiers)
  
  #Set up pcanon object
  CombinedSets <- vector(mode="list", length=5)
  names(CombinedSets) <- c("Q", "terms", "sources", "marginality", "aliasing")
  CombinedSets$Q <- CombinedSets$terms <- CombinedSets$sources <- 
    CombinedSets$marginality <- vector(mode="list", length=ntiers)
  if (form.named)
    names(CombinedSets$terms) <- names(CombinedSets$sources) <- 
       names(CombinedSets$marginality) <- names(formulae)
  
  #set up aliasing summary
  aliasing.cols <- c("Source", "df", "Alias", "In", criteria)
  nc <- length(aliasing.cols)
  aliasing <- data.frame(matrix(nrow = 0, ncol=nc))
  colnames(aliasing) <- aliasing.cols
  eff.crit <- vector(mode="list",length=length(criteria))
  names(eff.crit) <- criteria
  
  #Get set of projectors for the first formula
  if (ntiers == 1)
  {
    struct[[1]] <- pstructure(formulae[[1]], keep.order = keep.order, 
                              grandMean = grandMean,  orthogonalize = orthog, 
                              labels = labels, marginality = marginality[[1]], 
                              check.marginality = check.marginality, 
                              which.criteria = kcriteria, aliasing.print = aliasing.print,
                              data=data, ...)
    if ("pcanon" %in% omit.proj)
    {
      Qlabels <- names(struct[[1]]$Q)
      efficiencies <- vector(mode = "list", length = 0)
      for (i in Qlabels)
        struct[[1]]$Q[[i]] <- degfree(struct[[1]]$Q[[i]])
    }  
  }
  else 
    struct[[1]] <- pstructure(formulae[[1]], keep.order = keep.order, 
                              grandMean = grandMean,  orthogonalize = orthog[1], 
                              labels = labels, marginality = marginality[[1]], 
                              check.marginality = check.marginality, 
                              which.criteria = kcriteria, aliasing.print = aliasing.print, 
                              data=data, ...)

  #Load first structure into pcanon object
  CombinedSets$Q[[ntiers]] <- struct[[1]]$Q
  CombinedSets$terms[[1]] <- struct[[1]]$terms
  CombinedSets$sources[[1]] <- struct[[1]]$sources
  if (!is.null(struct[[1]]$marginality))
    CombinedSets$marginality[[1]] <- struct[[1]]$marginality
  if (!is.null(struct[[1]]$aliasing))
  {
    tmp <- struct[[1]]$aliasing
    if (form.named)
      tmp$In <- names(formulae)[1]
    else
      tmp$In <- "-"
    tmp <- tmp[aliasing.cols]
    aliasing <- rbind(aliasing, tmp)
  }
  lab <- attr(struct[[1]], which = "labels")
  #Loop over remaining  formulae
  if (ncomb > 0)
  { 
    for (k in 1:ncomb)
    { 
      ktier <- k + 1
      struct[[ktier]] <- pstructure(formulae[[ktier]], keep.order = keep.order, 
                                    grandMean = grandMean,  orthogonalize = orthog[ktier], 
                                    labels = labels, marginality = marginality[[ktier]], 
                                    check.marginality = check.marginality, 
                                    which.criteria = kcriteria, aliasing.print = aliasing.print,
                                    data=data, ...)
      if (!is.null(struct[[ktier]]$marginality))
        CombinedSets$marginality[[ktier]] <- struct[[ktier]]$marginality
      if (!is.null(struct[[ktier]]$aliasing))
      {
        tmp <- struct[[ktier]]$aliasing
        if (form.named)
          tmp$In <- names(formulae)[ktier]
        else
          tmp$In <- "-"
        tmp <- tmp[aliasing.cols]
        aliasing <- rbind(aliasing, tmp)
      }
      comb <- projs.2canon(CombinedSets$Q[[ntiers]], struct[[ktier]]$Q)
      CombinedSets$Q[[k]] <- comb$decomp
      if (!is.null(comb$aliasing))
        aliasing <- rbind(aliasing, comb$aliasing)
      CombinedSets$Q[[ntiers]] <- projs.combine.p2canon(comb)
      if ("pcanon" %in% omit.proj)
        CombinedSets$Q[[k]] <-  replace.proj(CombinedSets$Q[[k]])
      CombinedSets$terms[[ktier]] <- struct[[ktier]]$terms
      CombinedSets$sources[[ktier]] <- struct[[ktier]]$sources
      if (!is.null(struct[[ktier]]$marginality))
        CombinedSets$marginality[[ktier]] <- struct[[ktier]]$marginality
    }
    
    if ("combined" %in% omit.proj)
    { 
      comb.labels <- names(CombinedSets$Q[[ntiers]])
      for (i in comb.labels)
        CombinedSets$Q[[ntiers]][[i]] <- degfree(CombinedSets$Q[[ntiers]][[i]])
    }
    
  }
  #Set attribute for the number of units in the first tier
  attr(CombinedSets, which = "n") <- n
  
  if (nrow(aliasing) > 0)
  {
    CombinedSets$aliasing <- aliasing
    attr(CombinedSets$aliasing, which = "title") <- 
      "\nTable of information (partially) aliased with previous sources derived from the same formula\n\n"
    class(CombinedSets$aliasing) <- c("aliasing", "data.frame")
  }
  
  class(CombinedSets) <- "pcanon"
  attr(CombinedSets, which = "labels") <- lab
  
  return(CombinedSets)
}

#Function to form the summary from a pcanon.object
"summary.pcanon" <- function(object, labels.swap = FALSE, 
                             which.criteria = c("aefficiency", "eefficiency", "order"), 
                             ...)
  #Routine to output a summary of the projector analysis in a pcanon object
{ 
  if (!inherits(object, "pcanon"))
    stop("object must be of class pcanon as produced by projs.pcanon")
  if (!is.logical(labels.swap))
    stop("labels.swap must be a logical")
  if (is.null(names(object))) #Assume is an old pcannon
  {
    summ <- summ.legacy(object = object, which.criteria = which.criteria)
  } else
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
    
    ntiers <- length(object$Q)
    ncomb <- ntiers - 1
    #Check if have projector or df in pcanon object
    have.proj <- FALSE
    if (ntiers == 1)
    { 
      if (inherits(object$Q[[1]][[1]], "projector"))
        have.proj <- TRUE
    }  
    else 
      if (inherits(object$Q[[1]][[1]][["Q1res"]], "projector"))
        have.proj <- TRUE
    
    #Form data frame for summary table
    #Determine whether labels need changing
    old.labs <- attr(object, which = "labels")
    if (is.null(old.labs))
      new.labs <- "terms"
    else
      new.labs <- old.labs
    if (labels.swap & !is.null(old.labs) & all(c("terms","sources") %in% names(object)))
    {
      
      new.labs <- c("terms", "sources")[!(c("terms", "sources") %in% old.labs)]
      for (icomp in 1:ntiers)
      {
        if (icomp == 1)
        {
          #Convert uncombined first tier
          jold.labs <- object[[old.labs]][[1]]
          jnew.labs <- object[[new.labs]][[1]]
          ilabs <- names(object$Q[[1]])
          which.oldlabs <- match(ilabs, jold.labs)
          names(object$Q[[1]]) <- jnew.labs[which.oldlabs]
        } else
        {
          #Convert names indexing primary components
          ilabs.split <-lapply(names(object$Q[[icomp]]),  combined.split)
          for (jtier in 1:icomp)
          {
            #interchange labels for jtier
            jold.labs <- object[[old.labs]][[jtier]]
            jnew.labs <- object[[new.labs]][[jtier]]
            ilabs.split <- lapply(ilabs.split, 
                                  function(ilabs, jold.labs, jnew.labs)
                                  {
                                    which.oldlabs <- match(ilabs, jold.labs)
                                    which.ilabs <- which(!(is.na(which.oldlabs)))
                                    which.oldlabs <- na.omit(which.oldlabs)
                                    if (length(which.ilabs) > 0)
                                      ilabs[which.ilabs] <- jnew.labs[which.oldlabs]
                                    return(ilabs)
                                  }, jold.labs = jold.labs, jnew.labs = jnew.labs)
          }
          names(object$Q[[icomp]]) <- unlist(lapply(ilabs.split, makeCombinedSource))
        }
        
        #Convert the sublist names
        if (icomp < ntiers)
        {
          jtier <- icomp + 1
          jold.labs <- object[[old.labs]][[jtier]]
          jnew.labs <- object[[new.labs]][[jtier]]
          for (ks in 1:length(names(object$Q[[icomp]])))
          {
            ilabs <- names(object$Q[[icomp]][[ks]])
            which.oldlabs <- match(ilabs, jold.labs)
            which.ilabs <- which(!(is.na(which.oldlabs)))
            which.oldlabs <- na.omit(which.oldlabs)
            if (length(which.ilabs) > 0)
              names(object$Q[[icomp]][[ks]])[which.ilabs] <- jnew.labs[which.oldlabs]
          }
        }
      }
    }
    #Initialize
    tier.names <- names(object$terms)
    sources <- names(object$Q[[ntiers]])
    nlines <- length(sources)
    terms <- lapply(sources, combined.split)
    kcomb <- rep(1, nlines)
    ktiers <- lapply(terms, length)
    nc <- ntiers*2
    src.name <- "Source"
    if (new.labs == "terms")
      src.name <- "Term"
    if (ntiers == 1)
    {
      if (!is.null(tier.names))
        src.name <- paste(src.name, tier.names[1], sep = ".")
      srcdf.names <- c(src.name, " df")
    }
    else
    {
      if (is.null(tier.names))
        srcdf.names <- as.vector(outer(c(src.name, "df"), as.character(1:ntiers), paste, sep=""))
      else
        srcdf.names <- as.vector(sapply(1:ntiers,
                                        function(k, src.name, tier.names)
                                        {
                                          nam <- c(paste(src.name, tier.names[k], sep = "."),
                                                   paste("df", k, sep = ""))
                                        }, src.name = src.name, tier.names = tier.names))
    }
    orthogonaldesign <- TRUE
    if (anycriteria & ntiers > 1)
    { 
      nc <- nc + length(kcriteria)
      res.criteria <- vector(mode = "list", length = length(kcriteria))
      names(res.criteria) <- kcriteria
      res.criteria[kcriteria] <- NA
      summary <- data.frame(matrix(nrow = nlines, ncol=nc))
      colnames(summary) <- c(srcdf.names, kcriteria)
    }
    else
    {   
      summary <- data.frame(matrix(nrow = nlines, ncol=nc))
      colnames(summary) <- srcdf.names
    }
    #Loop over sources
    for (kl in 1:nlines)
      if (ntiers == 1)
      { #Get first tier source and df
        summary[[1]][kl] <- terms[[kl]][1]
        if (have.proj)
          summary[[2]][kl] <- degfree(object$Q[[1]][[terms[[kl]][1]]])
        else
          summary[[2]][kl] <- object$Q[[1]][[terms[[kl]][1]]]
      } else
      { #If next term is a residual, check it is for the previous term and output analogously
        if (terms[[kl]][ktiers[[kl]]] == "Residual")
        { #Check same as last term
          if (any(terms[[kl]][1:(ktiers[[kl]]-1)] != terms[[kl-1]][1:(ktiers[[kl]]-1)]))
            stop("Residual term and previous source are not confounded with the same source")
          label <- makeCombinedSource(terms[[kl]][1:(ktiers[[kl]]-1)])
          kcomb[kl] <- findLastwithSource(kcomb[kl], label, object$Q, ncomb)
          summary[kl, 1:(kcomb[kl]*2)] <- summary[(kl-1), 1:(kcomb[kl]*2)]
          summary[kl, (kcomb[kl]*2+1)] <- "Residual"
          if (have.proj)
            summary[[kcomb[kl]*2+2]][kl] <- degfree(object$Q[[kcomb[kl]]][[label]][["Q1res"]])
          else
            summary[[kcomb[kl]*2+2]][kl] <- object$Q[[kcomb[kl]]][[label]][["Q1res"]]
          if (ktiers[[kl]]  == 1)
            stop("Have a Residual in the first formula - include a factor for the units")
          else
          { 
            lastNonResidual <- max(c(1:ktiers[[kl]])["Residual" != terms[[kl]]])
            if (lastNonResidual  > 1)
            { 
              pcomb <- 1
              label <- makeCombinedSource(terms[[kl]][1:(lastNonResidual-1)])
              term2 <- terms[[kl]][lastNonResidual]
              pcomb <- findDecomposition(pcomb, term2, label, object$Q, ncomb) 
              if (abs(1 - unlist(object$Q[[pcomb]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
                orthogonaldesign <- FALSE
              if (anycriteria)
                summary[kl, kcriteria] <- unlist(object$Q[[pcomb]][[label]][[term2]][["adjusted"]][kcriteria])
            }
          }
        }
        else  #Not a Residual
        { #Get first tier source and df
          summary[[1]][kl] <- terms[[kl]][1]
          kconf.sources <- names(object$Q[[1]][[terms[[kl]][1]]])[-1]
          if (have.proj)
          { 
            totdf <- sum(unlist(lapply(kconf.sources, 
                                       function(ksrc, obj){ degfree(obj[[ksrc]][["Qproj"]])}, 
                                       obj = object$Q[[1]][[terms[[kl]][1]]])))
            totdf <- totdf + degfree(object$Q[[1]][[terms[[kl]][1]]][["Q1res"]])
          }
          else
          { 
            totdf <- sum(unlist(lapply(kconf.sources, 
                                       function(ksrc, obj){ obj[[ksrc]][["Qproj"]]}, 
                                       obj = object$Q[[1]][[terms[[kl]][1]]])))
            totdf <- totdf + object$Q[[1]][[terms[[kl]][1]]][["Q1res"]]
          }
          summary[[2]][kl] <- totdf
          label <- summary[[1]][kl]
          #Add remaining tiers
          if (ktiers[[kl]] > 1)
          { 
            for (kterm in 2:ktiers[[kl]])
            { 
              term2 <- terms[[kl]][kterm]
              if (term2 == "Residual")
              { 
                if (have.proj)
                  summary[[kcomb[kl]*2+2]][kl] <- degfree(object$Q[[kcomb[kl]]][[label]][["Q1res"]])
                else
                  summary[[kcomb[kl]*2+2]][kl] <- object$Q[[kcomb[kl]]][[label]][["Q1res"]]
              }
              else
              { 
                kcomb[kl] <- findDecomposition(kcomb[kl], term2, label, object$Q, ncomb)
                if (have.proj)
                  summary[[kcomb[kl]*2+2]][kl] <- degfree(object$Q[[kcomb[kl]]][[label]][[term2]][["Qproj"]])
                else
                  summary[[kcomb[kl]*2+2]][kl] <- object$Q[[kcomb[kl]]][[label]][[term2]][["Qproj"]]
              }
              summary[[kcomb[kl]*2+1]][kl] <- term2
              if (kterm != ktiers[[kl]])
              { 
                label <- paste(label,term2,sep="&")
                kcomb[kl] <- kcomb[kl] + 1
              }
            }
            if (abs(1 - unlist(object$Q[[kcomb[kl]]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
              orthogonaldesign <- FALSE
          }
          if (anycriteria & ktiers[[kl]] > 1)
            summary[kl, kcriteria] <- unlist(object$Q[[kcomb[kl]]][[label]][[term2]][["adjusted"]][kcriteria])
        }
      }
    titl <- "Summary table of the decomposition"
    if (!is.null(tier.names))
    {
      if (ntiers == 1)
        titl <- paste(titl, " for ", tier.names, sep = "")
      else
        titl <- paste(titl, " for ", 
                      paste(tier.names[1:(ntiers-1)], collapse = ", "),
                      " & ", tier.names[ntiers], sep = "")
    }
    if (!is.null(object$aliasing) | (ntiers > 1 & !orthogonaldesign))
      titl <- paste(titl, " (based on adjusted quantities)", sep = "")
    attr(summary, which = "title") <- paste("\n\n", titl, "\n\n", sep = "")
    attr(summary, which = "ntiers") <- ntiers
    attr(summary, which = "n") <- attr(object, which = "n", exact = "TRUE")
    attr(summary, which = "orthogonal") <- orthogonaldesign
    attr(summary, which = "labels") <- new.labs
    
    #Return summary and aliasing
    summ <- vector(mode = "list", length = 2)
    names(summ) <- c("decomp", "aliasing")
    summ$decomp <- summary
    if (!is.null(object$aliasing))
    {
      cols <- names(object$aliasing)[!(names(object$aliasing) %in% criteria)]
      if (anycriteria)
        cols <- c(cols, kcriteria)
      summ$aliasing <- object$aliasing[cols]
      attr(summ$aliasing, which = "title") <- attr(object$aliasing, which = "title")
    }
    class(summ) <- c("summary.pcanon", "list")
  }

  return(summ)
}

"summ.legacy" <- function(object, which.criteria = c("aefficiency", "eefficiency", "order"), ...)
  #Routine to output a summary of the projector analysis in a pcanon object
{ 
  if (!inherits(object, "pcanon"))
    stop("object must be of class pcanon as produced by projs.pcanon")
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  #Check if have projector or df in pcanon object
  have.proj <- FALSE
  if (inherits(object[[1]][[1]][["Q1res"]], "projector"))
    have.proj <- TRUE
  
  #Form data frame for summary table
  #Initialize
  ntiers <- length(object)
  ncomb <- ntiers - 1
  sources <- names(object[[ntiers]])
  nlines <- length(sources)
  terms <- lapply(sources, combined.split)
  kcomb <- rep(1, nlines)
  ktiers <- lapply(terms, length)
  nc <- ntiers*2
  srcdf.names <- as.vector(outer(c("Source", "df"), as.character(1:ntiers), paste, sep=""))
  orthogonaldesign <- TRUE
  if (anycriteria)
  { 
    nc <- nc + length(kcriteria)
    res.criteria <- vector(mode = "list", length = length(kcriteria))
    names(res.criteria) <- kcriteria
    res.criteria[kcriteria] <- NA
    summary <- data.frame(matrix(nrow = nlines, ncol=nc))
    colnames(summary) <- c(srcdf.names, kcriteria)
  }
  else
  {   
    summary <- data.frame(matrix(nrow = nlines, ncol=nc))
    colnames(summary) <- srcdf.names
  }
  #Loop over sources
  for (kl in 1:nlines)
  { #IF next term is a residual, check it is for the previous term and output analogously
    if (terms[[kl]][ktiers[[kl]]] == "Residual")
    { #Check same as last term
      if (any(terms[[kl]][1:(ktiers[[kl]]-1)] != terms[[kl-1]][1:(ktiers[[kl]]-1)]))
        stop("Residual term and previous source are not confounded with the same source")
      label <- makeCombinedSource(terms[[kl]][1:(ktiers[[kl]]-1)])
      kcomb[kl] <- findLastwithSource(kcomb[kl], label, object, ncomb)
      summary[kl, 1:(kcomb[kl]*2)] <- summary[(kl-1), 1:(kcomb[kl]*2)]
      summary[kl, (kcomb[kl]*2+1)] <- "Residual"
      if (have.proj)
        summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][["Q1res"]])
      else
        summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][["Q1res"]]
      if (ktiers[[kl]]  == 1)
        stop("Have a Residual in the first formula - include a factor for the units")
      else
      { 
        lastNonResidual <- max(c(1:ktiers[[kl]])["Residual" != terms[[kl]]])
        if (lastNonResidual  > 1)
        { 
          pcomb <- 1
          label <- makeCombinedSource(terms[[kl]][1:(lastNonResidual-1)])
          term2 <- terms[[kl]][lastNonResidual]
          pcomb <- findDecomposition(pcomb, term2, label, object, ncomb) 
          if (abs(1 - unlist(object[[pcomb]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
            orthogonaldesign <- FALSE
          if (anycriteria)
            summary[kl, kcriteria] <- unlist(object[[pcomb]][[label]][[term2]][["adjusted"]][kcriteria])
        }
      }
      #      if (anycriteria)
      #      { if (ktiers[[kl]]  == 1)
      #          stop("Have a Residual in the first formula - include a factor for the units")
      #        else
      #        { lastNonResidual <- max(c(1:ktiers[[kl]])["Residual" != terms[[kl]]])
      #          if (lastNonResidual  > 1)
      #          { pcomb <- 1
      #            label <- makeCombinedSource(terms[[kl]][1:(lastNonResidual-1)])
      #            term2 <- terms[[kl]][lastNonResidual]
      #            pcomb <- findDecomposition(pcomb, term2, label, object, ncomb) 
      #            summary[kl, kcriteria] <- unlist(object[[pcomb]][[label]][[term2]][["adjusted"]][kcriteria])
      #          }
      #        }
      #      }
    }
    else  #Not a Residual
    { #Get first tier source and df
      summary[[1]][kl] <- terms[[kl]][1]
      kconf.sources <- names(object[[1]][[terms[[kl]][1]]])[-1]
      if (have.proj)
      { 
        totdf <- sum(unlist(lapply(kconf.sources, 
                                   function(ksrc, obj){ degfree(obj[[ksrc]][["Qproj"]])}, 
                                   obj = object[[1]][[terms[[kl]][1]]])))
        totdf <- totdf + degfree(object[[1]][[terms[[kl]][1]]][["Q1res"]])
      }
      else
      { 
        totdf <- sum(unlist(lapply(kconf.sources, 
                                   function(ksrc, obj){ obj[[ksrc]][["Qproj"]]}, 
                                   obj = object[[1]][[terms[[kl]][1]]])))
        totdf <- totdf + object[[1]][[terms[[kl]][1]]][["Q1res"]]
      }
      summary[[2]][kl] <- totdf
      label <- summary[[1]][kl]
      #Add remaining tiers
      if (ktiers[[kl]] > 1)
      { 
        for (kterm in 2:ktiers[[kl]])
        { 
          term2 <- terms[[kl]][kterm]
          if (term2 == "Residual")
          { 
            if (have.proj)
              summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][["Q1res"]])
            else
              summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][["Q1res"]]
          }
          else
        { 
          kcomb[kl] <- findDecomposition(kcomb[kl], term2, label, object, ncomb)
          if (have.proj)
            summary[[kcomb[kl]*2+2]][kl] <- degfree(object[[kcomb[kl]]][[label]][[term2]][["Qproj"]])
          else
            summary[[kcomb[kl]*2+2]][kl] <- object[[kcomb[kl]]][[label]][[term2]][["Qproj"]]
        }
          summary[[kcomb[kl]*2+1]][kl] <- term2
          if (kterm != ktiers[[kl]])
          { 
            label <- paste(label,term2,sep="&")
            kcomb[kl] <- kcomb[kl] + 1
          }
        }
        if (abs(1 - unlist(object[[kcomb[kl]]][[label]][[term2]][["adjusted"]]["aefficiency"])) > 1e-04)
          orthogonaldesign <- FALSE
      }
      if (anycriteria & ktiers[[kl]] > 1)
        summary[kl, kcriteria] <- unlist(object[[kcomb[kl]]][[label]][[term2]][["adjusted"]][kcriteria])
    }
  }
  class(summary) <- c("summary.pcanon", "data.frame")
  attr(summary, which = "title") <- 
    "\n\nSummary table of the decomposition (based on adjusted quantities)\n\n"
  attr(summary, which = "ntiers") <- ntiers
  attr(summary, which = "orthogonal") <- orthogonaldesign
  attr(summary, which = "labels") <- "terms"
  return(summary)
}

#Function to print a summary from a pcanon object
print.summary.pcanon <- function(x, aliasing.print = TRUE, ...)
{ 
  if (!inherits(x, "summary.pcanon"))
    stop("Must supply an object of class summary.pcanon")
  
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                "order", "dforthog")
  
  if (inherits(x, "data.frame"))
    y <- x
  else
    y <- x$decomp
  cat(attr(y, which="title"))
  ntiers <- attr(y, which="ntiers")
  nlines <- nrow(y)
  if (ntiers ==1)
  { 
    dffw = max(3, floor(log10(max(y$' df'))) + 1)
    y$' df' <- formatC(y$' df', format="f", digits=0, width=dffw)
  } else
    dffw = max(3, floor(log10(max(y$df1))) + 1)
  repeats <- vector(mode="list", length=ntiers)
  labs <- attr(x, which = "labels")
  labs <- paste(toupper(substring(labs, 1, 1)), substring(labs, 2, nchar(labs)-1), sep = "")
  if (ntiers > 1)
  { 
    for (ktier in 1:ntiers)
    { 
      src.name <- colnames(y)[(ktier-1)*2+1] #paste(labs, ktier, sep="")
      if (all(is.na(y[src.name])))
        y[src.name] <- "  "
      else
      {
        if (ktier > 1)
          y[src.name][is.na(y[src.name])] <- "  "
        repeats[[ktier]] <- c(FALSE, y[2:nlines, src.name] == y[1:(nlines-1),src.name])
        repeats[[ktier]][is.na(repeats[[ktier]])] <- FALSE
        if (ktier > 1)
          repeats[[ktier]] <- (repeats[[ktier]] & repeats[[(ktier-1)]])
        y[repeats[[ktier]], src.name] <- "  "
      }
      if (ntiers == 1)
        df.name <- " df"
      else
        df.name <-paste("df", ktier, sep="")
      if (all(is.na(y[df.name])))
        y[df.name] <- "  "
      else
      {
        y[[df.name]] <- formatC(y[[df.name]], format="f", digits=0, width=dffw)
        y[repeats[[ktier]], df.name] <- "  "
        y[[df.name]] <- gsub("NA", "  ", y[[df.name]])
      }
    }
    for (kcrit in names(y)[names(y) %in% criteria])
    {
      if (all(is.na(y[kcrit])))
        y[kcrit] <- "  "
      else
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
    }
  }
  #Need to use print.data.frame for legacy summary objects
  print.data.frame(y, na.print="  ", right=FALSE, row.names=FALSE)

  if (aliasing.print & inherits(x, "list"))
    if (!is.null(x$aliasing))
      print(x$aliasing, which.criteria = names(x$aliasing)[names(x$aliasing) %in% criteria])

  n <- attr(y, which = "n", exact = "TRUE")
  if (!is.null(n))
  {
    df.sum <- sum(as.numeric(y[,2]), na.rm = TRUE)
    if (y[1,1] != "Mean")
      df.sum <- df.sum + 1
    if (n != df.sum)
      warning("The combined dimensions of the sources from the first formula", 
              " are less than the number of rows in data")
  }
  
  if (!attr(y, which="orthogonal"))
    cat("\nThe design is not orthogonal\n\n")
  invisible(x)
}

"efficiencies.pcanon" <- function(object, which = "adjusted", ...)
  #function to extract the efficiency factors from a pcanon object
{ if (!inherits(object, "pcanon"))
    stop("Must supply an object of class pcanon as produced by designAnatomy")
  options <- c("adjusted", "pairwise")
  opt <- options[check.arg.values(which, options)]
  
  #Get efficiencies
  ntiers <- length(object$Q)
  efficiencies <- vector(mode="list", length=(ntiers-1))
  for (k in 1:(ntiers-1))
     efficiencies[[k]] <- efficiencies.decomp(object$Q[[k]], which=opt)
  return(efficiencies)
}  

"marginality.pcanon" <- function(object, ...)
  #function to extract the marginality matrices from a pcanon object
{ if (!inherits(object, "pcanon"))
  stop("Must supply an object of class pcanon as produced by designAnatomy")

  #Get marginality matrices
  if (!("marginality" %in% names(object)) || is.null(object$marginality))
    marginality <- NULL
  else
  {
    marginality <- object$marginality
  }
  return(marginality)
}  

