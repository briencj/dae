"elements" <- function(x, subscripts)
#a function that returns, in a vector, the elements of x specified by the subscripts.
#x is an array
#subscripts is a two dimensional object interpreted as elements by dimensions
{ if (!(class(subscripts) == "matrix" | class(subscripts) == "data.frame"))
    stop("subscripts must be in a matrix or data.frame")
  n <- nrow(subscripts)
  nsub <- ncol(subscripts)
  if (nsub == 2)
  { sapply(1:n, function(i)x[subscripts[i,1],subscripts[i,2]])
  }
  else
    sapply(1:n, function(i)eval(parse(text=paste("x[", paste(subscripts[i,], collapse=","), "]"))))
}

"mat.I" <- function(order)
{ diag(rep(1, order))
}

"mat.J" <- function(order)
{ n <- order*order
  matrix(rep(1, n), nrow=order, ncol=order)
}

"mat.banded" <- function(x, nrow, ncol)
  #Function to form a banded matrix from a vector of values
  #- band 1 is the diagonal, band 2 the first subdiagonal and so on
{ nband <- length(x)
  if (nband > min(nrow, ncol))
    stop("Have supplied values for more than ",min(nrow, ncol) ," bands")
  matrix <- matrix(0, nrow=nrow, ncol=ncol)
  for (i in 1:nband)
    matrix[row(matrix)==col(matrix)+i-1 | row(matrix)+i-1 == col(matrix)] <- x[i]
  return(matrix)
}

"mat.ar1" <- function(rho, order)
#function to form the correlation matrix of size order with an ar1 pattern for
#correlation parameter rho
{ if (abs(rho) > 1)
    stop(paste("abs(rho) must be less than 1 \n",
                "Note order of mat.ar1 arguments changed to (rho, order)",sep=""))
  n <- order*order
  if (order < 2)
    stop("The order must be 2 or more")
  ar1 <- matrix(rep(rho, n), nrow=order, ncol=order)
  power <- abs(outer(1:order, 1:order, "-"))
  ar1 <- ar1^power
  return(ar1)
}

"mat.ar2" <- function(ARparameters, order)
{ #check ARparameters values
  if (length(ARparameters) <2)
    stop("Must supply two values in ARparameters")
  if (abs(ARparameters[1]) >= 1)
    stop("ARparameters[1] must be less than 1")
  if (abs(ARparameters[2]) >= 1)
    stop("ARparameters[2] must be less than 1")
  if (order < 3)
    stop("The order must be 3 or more")
  
  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- ARparameters[1]/(1-ARparameters[2])
  for (k in 3:order)
    corrs[k] <- ARparameters[1]*corrs[k-1] + ARparameters[2]*corrs[k-2]
  
  #Form correlation matrix
  ar2 <- mat.banded(corrs, order, order)  
  return(ar2)
}

"mat.ar3" <- function(ARparameters, order)
{ #check ARparameters values
  if (abs(ARparameters[1]) >= abs(1-ARparameters[2]))
    stop("ARparameters[1] must be less than 1 - ARparameters[2]")
  if (abs(ARparameters[2]) >= 1)
    stop("ARparameters[2] must be less than 1")
  if (abs(ARparameters[3]) >= 1)
    stop("ARparameters[3] must be less than 1")
  if (order < 4)
    stop("The order must be 4 or more")
  
  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  omega <- 1 - ARparameters[2] - ARparameters[3] * (ARparameters[1] + ARparameters[3])
  corrs[1] <- 1
  corrs[2] <- (ARparameters[1]+ARparameters[2]*ARparameters[3])/omega
  corrs[3] <- (ARparameters[1] * (ARparameters[1] + ARparameters[3]) + ARparameters[2] * (1 - ARparameters[2])) / omega
  for (k in 4:order)
    corrs[k] <- ARparameters[1]*corrs[k-1] + ARparameters[2]*corrs[k-2] + ARparameters[3]*corrs[k-3]
  
  #Form correlation matrix
  ar3 <- mat.banded(corrs, order, order)  
  return(ar3)
}

"mat.exp" <- function(rho, coordinates) 
  { if (abs(rho) > 1)
      stop(paste("abs(rho) must be less than 1 \n",
               "Note order of mat.exp arguments changed to (rho, coordinates)",sep=""))
    order <- length(coordinates)
    n <- order * order
    mat <- matrix(rep(rho, n), nrow = order, ncol = order)
    rownames(mat) <- colnames(mat) <- coordinates
    power  <- abs(outer(coordinates,coordinates,'-'))
    mat <- mat^power
    return(mat)
  }

"mat.gau" <- function(rho, coordinates) 
{ if (abs(rho) > 1)
  stop(paste("abs(rho) must be less than 1 \n",
             "Note order of mat.exp arguments changed to (rho, coordinates)",sep=""))
  order <- length(coordinates)
  n <- order * order
  mat <- matrix(rep(rho, n), nrow = order, ncol = order)
  rownames(mat) <- colnames(mat) <- coordinates
  power  <- (outer(coordinates,coordinates,'-'))^2
  mat <- mat^power
  return(mat)
}

"mat.sar" <- function(SARparameter, order)
{ #check SARparameter values
  if (abs(SARparameter) >= 1)
    stop("SARparameter must be less than 1")

  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- SARparameter/(1 + (SARparameter*SARparameter/4))
  for (k in 3:order)
    corrs[k] <- SARparameter*corrs[k-1] - (SARparameter*SARparameter/4)*corrs[k-2]
  
  #Form correlation matrix
  sar <- mat.banded(corrs, order, order)  
  return(sar)
}

"mat.sar2" <- function(gamma, order, print = NULL)
{ options <- c("ar3parameters")
  if (is.null(print))
    opt <- "none"
  else  
    opt <- options[unlist(lapply(print, check.arg.values, options=options))]
  
  #Calculate phi values
  phi <- vector(mode = "numeric", length = 3)
  phi[1] <- gamma[1] + 2*gamma[2]
  phi[2] <- -gamma[2] * (2*gamma[1] + gamma[2])
  phi[3] <- gamma[1] * gamma[2] * gamma[2]
  if (opt == "ar3parameters")
    cat(paste("Calculated  parameter values for the associated AR3 process\n", 
              paste(phi, collapse=", "),"\n"))
  
  #Calculate autocorrelations
  sar2 <- mat.ar3(phi, order)
  return(sar2)
}

"mat.ma1" <- function(MAparameter, order)
  #function to form the correlation matrix of size order with an ma1 pattern for
  #MAparameter
{ #check MAparameters value
  if (abs(MAparameter) > 1)
  stop("abs(MAparameter) must be less than 1 \n")
  n <- order*order
  if (order < 2)
    stop("The order must be 2 or more")

  #Calculate autocorrelations
  rho <- -MAparameter / (1 + MAparameter*MAparameter)
  
  #Form correlation matrix
  ma1 <- mat.banded(x = c(1, rho), nrow = order, ncol = order)
  return(ma1)
}

"mat.ma2" <- function(MAparameters, order)
  #function to form the correlation matrix of size order with an ma1 pattern for
  #MAparameter
{ #check MAparameters values
  if (length(MAparameters) <2)
    stop("Must supply two values in MAparameters")
  if (abs(MAparameters[1]) >= 1)
    stop("MAparameters[1] must be less than 1")
  if (abs(MAparameters[2]) >= 1)
    stop("MAparameters[2] must be less than 1")
  if ((MAparameters[1] + MAparameters[2] >= 1) | 
      (MAparameters[1] - MAparameters[2] >= 1))
    stop("The sum and difference of MAparameters must be less than 1")
  if (order < 3)
    stop("The order must be 3 or more")
  
  #Calculate autocorrelations
  div <- 1 + MAparameters[1]*MAparameters[1] + MAparameters[2]*MAparameters[2]
  rho1 <- -MAparameters[1]*(1 - MAparameters[2]) / div
  rho2 <- -MAparameters[2] / div

  #Form correlation matrix
  ma2 <- mat.banded(x = c(1, rho1, rho2), nrow = order, ncol = order)
  return(ma2)
}

"mat.arma" <- function(ARparameter, MAparameter, order)
  #function to form the correlation matrix of size order with an arma(1,1) pattern for
  #given ARparameter and MAparameter
{ #check ARparameter and MAparameter values
  if (abs(MAparameter) > 1)
    stop("abs(MAparameter) must be less than 1 \n")
  if (abs(ARparameter) > 1)
    stop("abs(ARparameter) must be less than 1 \n")
  n <- order*order
  if (order < 3)
    stop("The order must be 3 or more")
  
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- (MAparameter - ARparameter)*(1 - MAparameter*ARparameter) / 
                 (1 + MAparameter*MAparameter - 2*MAparameter*ARparameter)
  for (k in 3:order)
    corrs[k] <- ARparameter * corrs[k-1]

  #Form correlation matrix
  arma <- mat.banded(corrs, order, order)  
  return(arma)
}

"mat.dirprod" <- function(A, B)
#function create the direct product of A and B
{ rA <- nrow(A); cA <- ncol(A)
  rB <- nrow(B); cB <- ncol(B)
  Aexp <- A[rep(1:rA, each=rB), rep(1:cA, each=cB)]
  Bexp <- eval(parse(text=paste("cbind(", paste(rep("B", cA), collapse=","), ")")))
  Bexp <- eval(parse(text=paste("rbind(", paste(rep("Bexp", rA), collapse=","), ")")))
  Aexp*Bexp
}


"mat.dirsum" <- function(matrices)
  #Function to form the direct sum of a list of matrices
{ if (!is.list(matrices))
    stop("Must supply a list for matrices")
  if (!all(unlist(lapply(matrices, is.matrix))))
    stop("All elements of the matrices list must be matrices")
  nr <- lapply(matrices, nrow)
  nc <- lapply(matrices, ncol)
  m <- sum(unlist(nr))
  n <- sum(unlist(nc))
  dsum <- matrix(0, nrow=m, ncol=n)
  r1 <- r2 <- c1 <- c2 <- 0
  for (i in 1:length(matrices))
  { r1 <- r2 + 1
    c1 <- c2 + 1
    r2 <- r2 + nr[[i]]
    c2 <- c2 + nc[[i]]
    dsum[r1:r2, c1:c2] <- matrices[[i]]
  }
  return(dsum)  
}

#Function form variance matrix for random effects for a natural cubic smoothing spline
mat.ncssvar <- function(sigma2s = 1, knot.points, print = FALSE)
{
  #Check knot.points  
  if (is.matrix(knot.points))
    if (ncol(knot.points) != 1)
      stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  krange <- (h[r] - h[1]) / (r - 1)
  h <- diff(h, lag=1) / krange
  
  TD <- matrix(0, nrow = r, ncol = r)
  diag(TD)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
  TD[row(TD)==col(TD)-1] <- c(0,h[2:(r-2)], 0)/6
  TD[row(TD)==col(TD)+1] <- c(0,h[2:(r-2)], 0)/6
  TD <- TD[2:(r-1), 2:(r-1)]
  TD <- sigma2s * TD
  if (print)
  {
    cat("\n\n#### TD \n\n")
    print(TD)
  }
  
  return(TD)
}

#design matrix for a natural cubic smoothing splines
Zncsspline <- function(knot.points, Gpower = 0, print = FALSE)
{
  #Check knot.points  
  if (is.matrix(knot.points))
    if (ncol(knot.points) != 1)
      stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  krange <- (h[r] - h[1]) / (r - 1)
  h <- diff(h, lag=1) / krange
  
  delta <- diag(as.vector(1/h), nrow = r, ncol = (r-2))
  delta[row(delta)==col(delta)+1] <- -(1/h[1:(r-2)] + 1/h[2:(r-1)])
  delta[row(delta)==col(delta)+2] <- 1/h[2:(r-1)]
  Z <- delta %*% ginv(t(delta) %*% delta)
  if (print)
  {
    cat("\n\n#### delta\n\n")
    print(delta)
  }
  
  TD <- matrix(0, nrow = r, ncol = r)
  diag(TD)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
  TD[row(TD)==col(TD)-1] <- c(0,h[2:(r-2)], 0)/6
  TD[row(TD)==col(TD)+1] <- c(0,h[2:(r-2)], 0)/6
  TD <- TD[2:(r-1), 2:(r-1)]
  if (print)
  {
    cat("\n\n#### TD \n\n")
    print(TD)
  }
  
  #operator to get power of a matrix
  "%^%" <- function(x, n) 
    with(eigen(x), vectors %*% (values^n * t(vectors))) 
  
  if (Gpower != 0)  
    Z <- Z %*% (TD %^% Gpower)
  return(Z)
}

### Function to calculate the variance of predictions for Genotypes 
### based on Hookes (2009, Equation 17)
"mat.Vpred" <- function(W, Gg = 0, X = matrix(1, nrow = nrow(W), ncol = 1), Vu = 0, R, eliminate)
{ #set W or X to a column vector of 0s and Gg or Vu to a matrix of 0s if not effects for them not random
  warning("mat.Vpred is superseded by mat.Vpredicts, being retained for backwards compatibility; it may be deprecated in future versions")
  Gg.zero <- all(Gg < 1e-08)
  if (Gg.zero)
    Gginv <- Gg
  else
    Gginv <- ginv(Gg)
  Vinv <- ginv(Vu + R)
  if (!missing(eliminate))
  {
    if (!inherits(eliminate, "projector"))
      stop("Must supply an object of class projector")
    if (!Gg.zero)
      stop("Can only eliminate effects when the effects to be predicted are fixed")
    eliminate <- projector(diag(1, nrow = nrow(Vinv), ncol = nrow(Vinv)) - 
                             eliminate)
    Vinv <- eliminate %*% Vinv %*% eliminate
  }
  A <- t(W)%*%Vinv 
  Vpred <- A%*%W + Gginv
  if (!all(X < 1e-08))
  {
    AX <- A%*%X
    Vpred <- Vpred - AX%*%ginv(t(X)%*%Vinv%*%X)%*%t(AX)
  }
  Vpred <- ginv(Vpred)
  return(Vpred)
}

"mat.Vpredicts" <- function(target, Gt = 0, fixed = ~ 1, random, G, R, design, 
                            eliminate, keep.order = TRUE)
  #Gt is the component for the target factor; if zero the target treated as fixed, otherwise it is random
  #G is a list with a component for each random term; 
  #it can be a single value for the component value of a diagonal matrix, or
  #a matrix that is the same size as the number of levels in the model term;
  #the order in the list must correspond to the order of terms in the expanded formula
{
  method <- "onestep" #else twostep (as in Butler, 2013)
  #Generate the target matrix
  if (inherits(target, what = "matrix"))
    W <- target
  else
  {
    if (!inherits(target, what = "formula"))
      stop("target must be a matrix or a formula")
    W <- model.matrix(target, design, keep.order = keep.order)
  }
  #Set up Gt matrix
  if (inherits(Gt, what = "matrix"))
  {
    if (!all(dim(Gt) == ncol(W)))
      stop("Gt is a matrix that is not conformable with the target design matrix")
  } else
  {
    if (length(Gt) != 1)
      stop("Gt is not a matrix or a scalar")
    Gt <- diag(Gt, nrow = ncol(W), ncol = ncol(W))
  }
  #Get X from fixed
  if (inherits(fixed, what = "matrix"))
  {
    X <- fixed
  } else
  {
    if (!inherits(fixed, what = "formula"))
      stop("fixed must be a matrix or a formula")
    fix.attribs <- attributes(terms(fixed, keep.order = keep.order))
    terms <- fix.attribs$term.labels
    if (length(terms) != 0)
    {
      X <- do.call(cbind, lapply(terms, 
                                 function(term, design) 
                                   model.matrix(as.formula(paste("~ -1 +", term, sep = " ")), 
                                                design),
                                 design = design))
      if (fix.attribs$intercept != 0)
      {
        X <- cbind(matrix(1, nrow = nrow(W), ncol = 1), X)
        colnames(X)[1] <- "(Intercept)"    
      }
    } else
    {
      if (fix.attribs$intercept != 0)
        X <- matrix(1, nrow = nrow(W), ncol = 1)
      else
        X <- NULL
    }
  }
  fix.cols <- ncol(X)
  #Generate Z and G for the random terms when a formula is supplied
  if (missing(random))
  {
    if (!missing(G))
      stop("A G matrix has been specified without random having been set; perhaps it is the random matrix")
    else
      Vu <- matrix(0, nrow = nrow(W), ncol = nrow(W))
  } else #process the random argument
  {
    if (inherits(random, what = "matrix"))
    {
      if (method != "onestep")
        stop("If supply a matrix for random then method must be onestep")
      Vu <- random
    } else
    {
      if (!inherits(random, what = "formula"))
        stop("random must be a matrix or a formula")
      ran.attribs <- attributes(terms(random, keep.order = keep.order))
      if (missing(G))
        stop("Have specified a random formula without specifying components")
      if (!inherits(G, what = "list"))
        stop("G should be a list")
      if (ran.attribs$intercept != 0)
        stop("An intercept has been included in the random formula")
      nranterm <- length(ran.attribs$term.labels) + ran.attribs$intercept
      if (length(G) != nranterm)
        stop("The number of supplied components is not equal to the number of random terms")
      names(G) <- ran.attribs$term.labels
      terms <- ran.attribs$term.labels
      Z <- do.call(cbind, lapply(terms, 
                                 function(term, design) 
                                   model.matrix(as.formula(paste("~ - 1 +", term, sep = " ")), 
                                                design),
                                 design = design))
      
      #Generate G from component values
      for (term in ran.attribs$term.labels)
      {
        if (length(ran.attribs$factors) == 1)
        {
          facs <- rownames(ran.attribs$factors)
          nlev <- length(levels(design[[facs]]))
        } else
        {
          facs <- names(ran.attribs$factors[,term])[ran.attribs$factors[,term]  != 0]
          nlev <- prod(unlist(lapply(facs, 
                                     function(fac, design) length(levels(design[[fac]])), 
                                     design = design)))
        }
        if (length(G[[term]]) == 1)
        {
          G[[term]] <- diag(G[[term]], nrow = nlev, ncol = nlev)
        } else
        {
          if (!all(dim(G[[term]]) == nlev))
            stop("The dimensions of the G component for ",term, " is not correct")
        }
      }
      G <- mat.dirsum(G)
      if (!all(dim(G) == ncol(Z)))
        stop("The design matrix for the combined random terms is not conformable with G")
      if (method == "onestep")
        Vu <- Z %*% G %*% t(Z)
      else
        Vu = 0
    }
  }
  #Generate R if missing
  if (missing(R))
    R <- diag(1, nrow = nrow(W), ncol = nrow(W))
  
  if (method == "onestep")
  {
    #Now compute the variance matrix of the predictions
    Gt.zero <- all(Gt < 1e-08)
    if (Gt.zero)
      Gtinv <- Gt
    else
      Gtinv <- ginv(Gt)
    if (any(dim(Vu) != dim(R)))
      stop("The variance matrix for random effects has dimensions ",nrow(Vu),", ",ncol(Vu),
           " which is not conformable with R that has dimensions ", nrow(R),", ",ncol(R),
           "; also check target")
    Vinv <- ginv(Vu + R)
    if (!missing(eliminate))
    {
      if (!inherits(eliminate, "projector"))
        stop("Must supply an object of class projector")
      if (!Gt.zero)
        stop("Can only eliminate effects when the effects to be predicted are fixed")
      eliminate <- projector(diag(1, nrow = nrow(Vinv), ncol = nrow(Vinv)) - 
                               eliminate)
      Vinv <- eliminate %*% Vinv %*% eliminate
    }
    A <- t(W)%*%Vinv 
    Vpred <- A%*%W + Gtinv
    if (!all(X < 1e-08))
    {
      AX <- A%*%X
      Vpred <- Vpred - AX%*%ginv(t(X)%*%Vinv%*%X)%*%t(AX)
    }
    Vpred <- ginv(Vpred)
  } else #twostep
  {
    target.cols <- ncol(W)
    
    if (Gt == 0) #add to fixed model
    {
      X <- cbind(X, W)
    }
    else #add to random model
    {
      G <- mat.dirsum(list(Gt, G))
      if (!missing(random))
        Z <- cbind(W, Z)
      else
        Z <- W
    }
    if (!missing(random) || Gt !=0)
      Ginv <- mat.dirsum(list(diag(0, nrow = ncol(X), ncol = ncol(X)), 
                              ginv(G)))
    
    #Form the information matrix
    C <- cbind(X, Z)
    C <- t(C) %*% ginv(R) %*% C + Ginv
    V.p <- ginv(C)
    
    #Extract target
    subcols <- fix.cols + 1:target.cols
    #  V.p <- V.p[subcols, subcols]
    
    cols.nonT <- c(1:fix.cols, (fix.cols+target.cols+1):ncol(C))
    C11 <- C[subcols, subcols]
    C22 <- C[cols.nonT, cols.nonT]
    C12 <- C[subcols, cols.nonT]
    Cadj <- C11 - (C12 %*% ginv(C22) %*% t(C12))
    Vpred <- ginv(Cadj)
  }
  return(Vpred)
}

### Function to calculate 
"designAmeasures" <- function(Vpred, groupsizes = NULL, groups = NULL)
{ 
  #determine any groupings of the variances
  n <- nrow(Vpred)
  if (n != ncol(Vpred))
    stop("Vpred must be square")  
  if (is.null(groupsizes) & is.null(groups))
  {
    groupsizes <- n
    groups <- list(all = 1:n)
    ngrp <- 1
  } else
  {
    if (!is.null(groups)) #working on groups
    {
      if (!is.list(groups))
        stop("the groups argument should be a list")
      if (!is.null(groupsizes))
        warning("Both groups and grouspsizes are set - using groups")
      ngrp <- length(groups)
      groupsizes <- unlist(lapply(groups, length))
    } else #working on group sizes
    {
      ngrp <- length(groupsizes)
      end <- cumsum(groupsizes)
      if (end[ngrp] > n)
        stop("The sum of the group sizes must less than or equal to the size of Vpred")
      if(ngrp == 1)
        start <- 1
      else
        start <- c(0, end[-ngrp]) +1
      groups <- mapply(function(start,end)
                             {start:end},
                       start, end, SIMPLIFY = FALSE)
      if (!is.null(names(groupsizes)))
        names(groups) <- names(groupsizes)
    }
  }
  
  # Calculate the variance matrix for differences between predictions
  varDiff <- matrix(rep(diag(Vpred), each = n), nrow = n) + 
             matrix(rep(diag(Vpred), times = n), nrow = n) - 2 * Vpred
  
  # Calculate all within and between group A-measure values
  A <- matrix(0, nrow = ngrp, ncol = ngrp)
  for (i in 1:ngrp)
  {
    for(j in i:ngrp)
    {
      A[i, j] <- sum(varDiff[groups[[i]], groups[[j]]])
      if (i==j)
        A[i, i] <- A[i, i]/(groupsizes[i]*(groupsizes[i]-1))
      else
      {
        A[i, j] <- A[i, j]/(groupsizes[i]*groupsizes[j])
        A[j, i] <- A[i, j]
      }
    }
  }
  if (!is.null(names(groups)))
  {
    rownames(A) <- colnames(A) <- names(groups)
  }
  return(A)
}
