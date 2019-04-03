
"efficiency.criteria" <- function(efficiencies)
{ 
  daeTolerance <- get("daeTolerance", envir=daeEnv)
  criteria <- vector(mode="list", length = 7)
  names(criteria) <- c('aefficiency','mefficiency','sefficiency','eefficiency',
                       "xefficiency",'order',"dforthog")
  df.orthog <- sum(abs(1-efficiencies) <  daeTolerance[["eigen.tol"]])
  eff.unique <- remove.repeats(efficiencies, daeTolerance[["eigen.tol"]])
  K <- length(eff.unique)
  if (K == 1)
  { if (eff.unique == 0)
    { 
    criteria["aefficiency"] <- 0
      criteria["mefficiency"] <- 0
      criteria["sefficiency"] <- 0
      criteria["eefficiency"] <- 0
      criteria["xefficiency"] <- 0
      criteria["order"] <- 0
      criteria["dforthog"] <- 0
  } else
    { 
      criteria["aefficiency"] <- eff.unique
      criteria["mefficiency"] <- eff.unique
      criteria["sefficiency"] <- 0
      criteria["eefficiency"] <- eff.unique
      criteria["xefficiency"] <- eff.unique
      criteria["order"] <- K
      criteria["dforthog"] <- df.orthog
    }
  } else
  { 
    criteria["aefficiency"] <- harmonic.mean(efficiencies)
    criteria["mefficiency"] <- mean(efficiencies)
    criteria["sefficiency"] <- var(efficiencies)
    criteria["eefficiency"] <- min(efficiencies)
    criteria["xefficiency"] <- max(efficiencies)
    criteria["order"] <- K
    criteria["dforthog"] <- df.orthog
  }
  return(criteria)
}

print.summary.p2canon <- function(x, ...)
{ if (!inherits(x, "summary.p2canon"))
    stop("Must supply an object of class summary.p2canon")
  cat(attr(x, which="title"))
  y <- x
  nlines <- nrow(y)
  repeats <- c(FALSE, y[2:nlines,"Source"] == y[1:(nlines-1),"Source"])
  y[repeats, "Source"] <- "  "
  if ("aefficiency" %in% names(y)) { 
    y$aefficiency <- formatC(y$aefficiency, format="f", digits=4, width=11)
    y$aefficiency <- gsub("NA", "  ", y$aefficiency)
  }
  if ("mefficiency" %in% names(y)) { 
    y$mefficiency <- formatC(y$mefficiency, format="f", digits=4, width=11)
    y$mefficiency <- gsub("NA", "  ", y$mefficiency)
  }
  if ("eefficiency" %in% names(y)) { 
    y$eefficiency <- formatC(y$eefficiency, format="f", digits=4, width=11)
    y$eefficiency <- gsub("NA", "  ", y$eefficiency)
  }
  if ("xefficiency" %in% names(y)) { 
    y$xefficiency <- formatC(y$xefficiency, format="f", digits=4, width=11)
    y$xefficiency <- gsub("NA", "  ", y$xefficiency)
  }
  if ("sefficiency" %in% names(y)) { 
    y$sefficiency <- formatC(y$sefficiency, format="f", digits=4, width=11)
    y$sefficiency <- gsub("NA", "  ", y$sefficiency)
  }
  if ("order" %in% names(y)) { 
    y$order <- formatC(y$order, format="f", digits=0, width=5)
    y$order <- gsub("NA", "  ", y$order)
  }
  if ("dforthog" %in% names(y)) { 
    y$dforthog <- formatC(y$dforthog, format="f", digits=0, width=8)
    y$dforthog <- gsub("NA", "  ", y$dforthog)
  }
  print.data.frame(y, na.print="  ", right=FALSE, row.names=FALSE)
  if (!attr(x, which="orthogonal"))
    cat("\nThe design is not orthogonal\n\n")
  invisible(x)
}

"proj2.sweep" <- function(Q1, Q2, Eff.Q1.Q2)
{ #A procedure to compute the Residual operator for P remove Q when P and Q are nonorthogonal.
  #  Corresponding projection operator for Q in P is also obtained.
  #This version, which allows the supply of efficiency factors, is not being used
  n <- nrow(Q1)
  if (n != nrow(Q2))
    stop("Matrices not conformable.")
  isproj <- is.projector(Q1) & is.projector(Q2)

  if (length(Eff.Q1.Q2) == 1 & Eff.Q1.Q2[1]==0) #check matrices are orthogonal
  { 
    Qconf <- projector(matrix(0, nrow = n, ncol = n))
    Qres <- Q1
    Eff.Q1.Q2 <- 0
    warning("Matrices are orthogonal.")
  } else
  { 
    daeTolerance <- get("daeTolerance", envir=daeEnv)
    EffUnique.Q1.Q2 <-remove.repeats(Eff.Q1.Q2, daeTolerance[["eigen.tol"]])
    K <- length(EffUnique.Q1.Q2)
    #check for all confounded (i.e. eff = 1)
    if (K == 1 & EffUnique.Q1.Q2[1] == 1 & length(Eff.Q1.Q2) == degfree(Q2))
    { 
      Qconf <- projector(Q2)
      Qres <- projector(Q1 - Q2)
    } else      #compute projection operators for partially confounded case
    { 
      I <- diag(1, nrow = n, ncol = n)
      Qres <- I
      Q121 <- Q1 %*% Q2 %*% Q1
      for(i in 1:K)
        Qres <- Qres %*% (Q1 - (Q121/EffUnique.Q1.Q2[i]))
      Qres <- projector(Qres)
      Qconf <- projector(Q1 - Qres)
    }
  }
  list(Qconf = Qconf, Qres = Qres)
}


"projs.2canon" <- function(Q1, Q2)
  #Function to do an eigenanalysis of the relationship between two sets of projection matrices
{ 
  if (!is.list(Q1) | !is.list(Q2))
    stop("Both Q1 and Q2 must be lists")
  daeTolerance <- get("daeTolerance", envir=daeEnv)
  
  #Get sizes and set up labels
  nQ1 <- length(Q1)
  nQ2 <- length(Q2)
  n <- nrow(Q1[1])
  Q1labels <- names(Q1)
  if (is.null(Q1labels))
    Q1labels <- as.character(1:nQ1)
  Q2labels <- names(Q2)
  if (is.null(Q2labels))
    Q2labels <- as.character(1:nQ2)
  criteria <- c('aefficiency','mefficiency','sefficiency','eefficiency',"xefficiency",
                'order',"dforthog")
  
  #set up aliasing summary
  nc <- 4 + length(criteria)
  aliasing <- data.frame(matrix(nrow = 0, ncol=nc))
  colnames(aliasing) <- c("Source", "df", "Alias", "In", criteria)

  #Perform the analysis
  kQ1Q2 <- 0
  multieffic <- FALSE
  results <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { 
    results[[i]][["Q1res"]] <- Q1[[i]]
    rdf <- degfree(Q1[[i]])
    for (j in Q2labels)
    { 
      if (rdf >0) #only do this if there are df left in Q1
      { #Get unadjusted criteria and store only if confounded
        Q1Q2.eff <- suppressWarnings(proj2.efficiency(Q1[[i]], Q2[[j]]))
        if (Q1Q2.eff[1] > 0)
        { 
          #Get adjusted efficiencies and matrices
          adj.Q1Q2 <- suppressWarnings(proj2.combine(results[[i]][["Q1res"]], Q2[[j]]))
          if (degfree(adj.Q1Q2$Qconf) > 0)
          { 
            #Check for adjusted orthogonality
            Qfitlab <- names(results[[i]])[-1]
            if (length(Qfitlab) > 0)
            { 
              for (k in Qfitlab)
              { 
                Qjik <- Q2[[k]] %*% Q1[[i]] %*% Q2[[j]]
                if(!is.allzero(Qjik))
                {
                  warning(paste(j,"and",k,"are partially aliased in",i, sep=" "))
                  eff.crit <- efficiency.criteria(adj.Q1Q2$efficiencies)
                  names(eff.crit) <- criteria
                  aliasing <- rbind(aliasing, 
                                    data.frame(c(list(Source = j, 
                                                      df = degfree(adj.Q1Q2$Qconf),
                                                      Alias = k,
                                                      In = i),
                                                 eff.crit[criteria]), 
                                               stringsAsFactors = FALSE))
                }
              }
            }  
            
            results[[i]][[j]] <- vector(mode = "list", length = 0)
            #store pairwise efficiencies
            results[[i]][[j]][["pairwise"]][["efficiencies"]] <- Q1Q2.eff
            results[[i]][[j]][["pairwise"]][criteria] <- efficiency.criteria(Q1Q2.eff)
            
            if (degfree(adj.Q1Q2$Qres) == 0)
            { 
              adj.Q1Q2$Qres <- matrix(0, nrow=nrow(adj.Q1Q2$Qres), ncol=ncol(adj.Q1Q2$Qres))
              adj.Q1Q2$Qres <- projector(adj.Q1Q2$Qres)
            }
            results[[i]][["Q1res"]] <- adj.Q1Q2$Qres
          
            #store adjusted efficiencies
            results[[i]][[j]][["adjusted"]][["efficiencies"]] <- adj.Q1Q2$efficiencies
            results[[i]][[j]][["adjusted"]][criteria] <- efficiency.criteria(adj.Q1Q2$efficiencies)
            if (results[[i]][[j]][["adjusted"]][["aefficiency"]] - 
                  results[[i]][[j]][["adjusted"]][["eefficiency"]] > daeTolerance[["eigen.tol"]])
              multieffic <- TRUE
            
            #Store adjusted projector
            results[[i]][[j]][["Qproj"]] <- adj.Q1Q2$Qconf
            rdf <- degfree(adj.Q1Q2$Qres)
          } else
          {
            warning(paste(j,"is aliased with previous terms in",i, sep=" "))
            eff.crit <- rep(0, length(criteria))
            names(eff.crit) <- criteria
            aliasing <- rbind(aliasing, 
                              data.frame(c(list(Source = j, 
                                                df = 0,
                                                Alias = "unknown",
                                                In = i),
                                          eff.crit), 
                                         stringsAsFactors = FALSE))
          }
        }  
      }
    }
  }
  if (nrow(aliasing) == 0 )
    aliasing <- NULL
  p2can <- list(decomp = results, aliasing = aliasing)
  class(p2can) <- "p2canon"
  return(p2can)
}

"summary.p2canon" <- function(object, which.criteria = c("aefficiency", "eefficiency", "order"), ...)
  #Routine to output a summary of the projector analysis
{ 
  if (!inherits(object, "p2canon"))
    stop("object must be of class p2canon as produced by projs.2canon")
  #check which.criteria arguments
  criteria <- c("aefficiency", "eefficiency", "mefficiency", "sefficiency", "xefficiency", 
                "order", "dforthog")
  options <- c(criteria, "none", "all")
  kcriteria <- options[unlist(lapply(which.criteria, check.arg.values, 
                                     options=options))]
  if ("all" %in% kcriteria)
    kcriteria <- criteria
  anycriteria <- !("none" %in% kcriteria)
  
  #Form data frame for summary table
  orthogonaldesign <- TRUE
  nc <- 3
  if (anycriteria)
  { 
    nc <- nc + length(kcriteria)
    res.criteria <- vector(mode = "list", length = length(kcriteria))
    names(res.criteria) <- kcriteria
    res.criteria[kcriteria] <- NA
    summary <- data.frame(matrix(nrow = 0, ncol=nc))
    colnames(summary) <- c("Source", "Confounded.source", "df", kcriteria)
  }
  else
  {   summary <- data.frame(matrix(nrow = 0, ncol=nc))
      colnames(summary) <- c("Source", "Confounded.source", "df")
  }  
  Q1labels <- names(object$decomp)
  for (i in Q1labels)
  { 
    Q2labels <- names(object$decomp[[i]])[-1]
    nconf.terms <- 0
    if (length(Q2labels) > 0)
    { 
      for (j in Q2labels)
      { 
      kdf <- degfree(object$decomp[[i]][[j]]$Qproj)
        if (kdf > 0)
        { 
          nconf.terms <- nconf.terms + 1
          if (abs(1 - unlist(object$decomp[[i]][[j]][["adjusted"]]["aefficiency"])) > 1e-04)
            orthogonaldesign <- FALSE
          if (anycriteria)
            summary <- rbind(summary,
                             data.frame(Source = i,
                                        Confounded.source = j, 
                                        df = kdf, 
                                        object$decomp[[i]][[j]][["adjusted"]][kcriteria], 
                                        stringsAsFactors = FALSE))
          else
            summary <- rbind(summary,
                             data.frame(Source = i,
                                        Confounded.source = j, 
                                        df = kdf, 
                                        stringsAsFactors = FALSE))
        }
      }
    }
    kdf <- degfree(object$decomp[[i]]$Q1res)
    if (kdf > 0)
      if (nconf.terms > 0)
      { 
        if (anycriteria)
        summary <- rbind(summary,
                         data.frame(Source = i,
                                    Confounded.source = "Residual", 
                                    df = kdf, 
                                    res.criteria[kcriteria], 
                                    stringsAsFactors = FALSE))
        else
          summary <- rbind(summary,
                           data.frame(Source = i,
                                      Confounded.source = "Residual", 
                                      df = kdf, 
                                      stringsAsFactors = FALSE))
      } else
      { 
        if (anycriteria)
        summary <- rbind(summary,
                         data.frame(Source = i,
                                    Confounded.source = "   ", 
                                    df = kdf, 
                                    res.criteria[kcriteria], 
                                    stringsAsFactors = FALSE))
        else
          summary <- rbind(summary,
                           data.frame(Source = i,
                                      Confounded.source = "   ", 
                                      df = kdf, 
                                      stringsAsFactors = FALSE))
      }
  }
  summary <- as.data.frame(summary)
  class(summary) <- c("summary.p2canon", "data.frame")
  titl <- "\n\nSummary table of the decomposition"
  if (!orthogonaldesign)
    titl <- paste(titl, " (based on adjusted quantities)\n\n", sep = "")
  else
    titl <- paste(titl, "\n\n", sep = "")
  attr(summary, which = "title") <- titl
  attr(summary, which = "orthogonal") <- orthogonaldesign
  return(summary)
}


"efficiencies.decomp" <- function(decomp, which = "adjusted", ...)
  #function to extract the efficiency factors from a p2canon object 
{ 
  options <- c("adjusted", "pairwise")
  opt <- options[check.arg.values(which, options)]
  
  #Get efficiencies
  Q1labels <- names(decomp)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { 
    Q2labels <- names(decomp[[i]])[-1]
    if (length(Q2labels) > 0)
    { 
      efficiencies[[i]] <- vector(mode = "list", length = 0)
      for (j in Q2labels)
        efficiencies[[i]][[j]] <- decomp[[i]][[j]][[opt]][["efficiencies"]]
    }  
  }
  return(efficiencies)
}  

"efficiencies.p2canon" <- function(object, which = "adjusted", ...)
#function to extract the efficiency factors from a p2canon object 
{ 
  if (!inherits(object, "p2canon"))
    stop("object must be of class p2canon as produced by projs.2canon")
  options <- c("adjusted", "pairwise")
  opt <- options[check.arg.values(which, options)]
  
  #Get efficiencies
  efficiencies <- efficiencies.decomp(object$decomp, which = which)
  return(efficiencies)
}  

"projs.combine.p2canon" <- function(object)
#function to extract, from a p2canon object, the projectors that give the combined decomposition  
{ 
  #Initialize list
  Q1combineQ2 <- vector(mode = "list", length = 0)
  
  #Loop through p2canon object
  Q1labels <- names(object$decomp)
  efficiencies <- vector(mode = "list", length = 0)
  for (i in Q1labels)
  { 
    Q2labels <- names(object$decomp[[i]])[-1]
    if (length(Q2labels) > 0)
    { 
      for (j in Q2labels)
      { 
        Q1Q2label <- paste(Q1labels[[match(i, Q1labels)]], Q2labels[[match(j, Q2labels)]], sep="&")
        Q1combineQ2[[Q1Q2label]] <-object$decomp[[i]][[j]][["Qproj"]]
      }
      #Get the residual if any
      if (degfree(object$decomp[[i]][["Q1res"]]) > 0)
      { 
        Q1Q2label <- paste(Q1labels[[match(i, Q1labels)]],"Residual", sep="&")
        Q1combineQ2[[Q1Q2label]] <- object$decomp[[i]][["Q1res"]]
      }
    } else
      Q1combineQ2[[i]] <- object$decomp[[i]][["Q1res"]]
  }
  return(Q1combineQ2)
}  
