designTwophaseAnatomies <- function(formulae, data, which.designs = "all", 
                                    printAnatomies = TRUE, titles,
                                    orthogonalize = "hybrid", 
                                    marginality = NULL, 
                                    which.criteria = c("aefficiency", "eefficiency", 
                                                       "order"), ...)
{
  ntiers <- length(formulae)
  if (ntiers != 3 || !all(unlist(lapply(formulae, plyr::is.formula))))
    stop("must supply a list with three formulas")  
  if (length(orthogonalize) == 1)
      orthogonalize <- rep(orthogonalize, ntiers)
    else
      if (length(orthogonalize) != ntiers)
      { warning("Length of orthogonalize is not equal to 1 or the number of formulae - only using first value")
        orthogonalize <- rep(orthogonalize[1], ntiers)
      }
  
  #Deal with titles argument
  anat.titls <- c("Anatomy for the full two-phase design",
                  "Anatomy for the first-phase design",
                  "Anatomy for the cross-phase, treatments design",
                  "Anatomy for the combined-units design")
  if (!missing(titles))
  {
    if (length(titles) != 4 | !is.character(titles))
      stop("titles must be a character of length fours")
    anat.titls <- mapply(FUN = function(titl, anat.titl) 
    {
      if (is.na(titl)) titl <- anat.titl
      invisible(titl)
    }, 
    titles, anat.titls, USE.NAMES = FALSE)
  }

  #check which.criteria arguments
  designs <- c("two-phase", "first-phase", "cross-phase", "combined-units")
  options <- c(designs, "all")
  kdesigns <- options[unlist(lapply(which.designs, check.arg.values, 
                                     options=options))]
  if ("all" %in% which.designs)
    kdesigns <- designs

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
  
  twoph.lay.canon <- twoph1.lay.canon <- ph12.lay.canon <- twoph2.lay.canon <- NULL
  
  #'### Anatomy for full two-phase design
  if ("two-phase" %in% kdesigns)
  {
    twoph.lay.canon <- designAnatomy(formulae = formulae, data = data, 
                                     orthogonalize = orthogonalize, 
                                     marginality = marginality, 
                                     which.criteria = which.criteria, ...)
    if (printAnatomies)
    { 
      cat("\n###", anat.titls[1],"\n")
      print(summary(twoph.lay.canon, which.criteria = which.criteria))
    }
  }
  
  #'### Anatomy for first-phase design
  if ("first-phase" %in% kdesigns)
  {
    twoph1.lay.canon <- designAnatomy(formulae = c(formulae[2], formulae[3]),
                                      data = data, 
                                      orthogonalize = orthogonalize[c(2,3)], 
                                      marginality = c(marginality[2],marginality[3]), 
                                      which.criteria = which.criteria, ...)
    if (printAnatomies)
    {
      cat("\n###", anat.titls[2],"\n")
      print(summary(twoph1.lay.canon, which.criteria = which.criteria))
    }
  }
  
  #'### Anatomy for the cross-phase, treatments design
  if ("cross-phase" %in% kdesigns)
  {
    ph12.lay.canon <- designAnatomy(formulae = c(formulae[1], formulae[3]),
                                    data = data, 
                                    orthogonalize = orthogonalize[c(1,3)], 
                                    marginality = c(marginality[1],marginality[3]), 
                                    which.criteria = which.criteria, ...)
    if (printAnatomies)
    {
      cat("\n###", anat.titls[3],"\n")
      print(summary(ph12.lay.canon, which.criteria = which.criteria))
    }
  }
  
  #'### Anatomy for the cpmbined-units design
  if ("combined-units" %in% kdesigns)
  {
    twoph2.lay.canon <- designAnatomy(formulae = c(formulae[1], formulae[2]),
                                      data = data, 
                                      orthogonalize = orthogonalize[c(1,2)], 
                                      marginality = c(marginality[1],marginality[2]), 
                                      which.criteria = which.criteria, ...)
    if (printAnatomies)
    {
      cat("\n###", anat.titls[4],"\n")
      print(summary(twoph2.lay.canon, which.criteria = which.criteria))
    }
  }  
  anats <- list(twophase = twoph.lay.canon, first = twoph1.lay.canon, cross = ph12.lay.canon, 
                units = twoph2.lay.canon)
  attr(anats, which = "titles") <- anat.titls
  invisible(anats)
}
