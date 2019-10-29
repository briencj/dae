#devtools::test("dae")
context("Ameasures")

cat("#### Test for mt.Vpred and designAmeasures using Smith's small example\n")
test_that("SmithSmall", {
  skip_on_cran()
  library(dae)
  #'# Script to investigate the reduced design from Smith et al. (2015)

  #'## Set up designs
  mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
  field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
  field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                              levels = c("Y","W","G","M","D","E"))
  start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
  rownames(start.design) <- NULL
  start.design

  #'## Set gammas
  terms <- c("Variety", "Frep", "Frep:Fplot", "Mrep", "Mrep:Mday", "Mrep:Mday:Mord")
  gammas <- c(1, 0.1, 0.2, 0.3, 0.2, 1)
  names(gammas) <- terms
  print(gammas)
  
  #'## Examine A-optimality measures
  #'### Set up functions
  "invInf" <- function(design, gammas, ...)
  {
    #'## Set up matrices
    W <- model.matrix(~ -1 + Variety, design)
    ZMr <- model.matrix(~ -1 + Mrep, design)
    ZMrMd <- model.matrix(~ -1 + Mrep:Mday, design)
    ZFr <- model.matrix(~ -1 + Frep, design)
    ZFrFp <- model.matrix(~ -1 + Frep:Fplot, design)
    t <- length(levels(design$Variety))
    ng <- ncol(W)
    Vu <- gammas["Mrep"] * ZMr %*% t(ZMr) + gammas["Mrep:Mday"] * ZMrMd %*% t(ZMrMd) + 
      gammas["Frep"] * ZFr %*% t(ZFr) + gammas["Frep:Fplot"] * ZFrFp %*% t(ZFrFp)
    R <- diag(1, nrow(design))
    
    #'# Calculate variance matrix
    Vp <- mat.Vpred(W = W, Vu=Vu, R=R, ...)
    
    return(Vp)
  }
  
  
  #'## Calculate A measures
  Vp <- invInf(start.design, gammas)
  testthat::expect_true(all(abs(designAmeasures(Vp) - 1.521905) < 1e-6))
  n <- nrow(start.design)
  Vp.reduc <- invInf(start.design, gammas, 
                     eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
  testthat::expect_true(all(abs(rowSums(Vp.reduc)) < 1e-10))
  testthat::expect_true(abs(designAmeasures(Vp.reduc) - 
                              sum(diag(Vp.reduc)) * 2 / (nrow(Vp.reduc) - 1)) < 1.e-10)
  
  #'## Investigate the anatomy
  start2.canon <- designAnatomy(list(lab = ~ Mrep/Mday/Mord, 
                                     plot = ~ Frep/Fplot),
                                data = start.design, omit.projectors = "pcanon")
  testthat::expect_equal(degfree(start2.canon$Q[[2]]$`Mord[Mrep:Mday]&Fplot[Frep]`), 5)
})  

cat("#### Test for Ameasures using Cochran&Cox PBIBD2\n")
test_that("PBIBD2Ameasures", {
  skip_on_cran()
  library(dae)
  #'# PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 2nd edn Wiley, New York"
  
  #'## Input the design and randomize"
  Treatments <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
  PBIBD2.lay <- designRandomize(allocated = Treatments,
                                recipient = list(Blocks =6, Units = 4),
                                nested.recipients = list(Units = "Blocks"),
                                seed = 98177)
  
  #'## Compute the anatomy
  PBIBD2.canon <- designAnatomy(formulae = list(unit = ~Blocks/Units,
                                                trt = ~ Treatments),
                                which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                data = PBIBD2.lay, omit.projectors = "combined")
  
  
  testthat::expect_equal(PBIBD2.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2.canon$Q[[1]]$`Units[Blocks]`$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  
  #'## Get Ameasure
  #'## Set up matrices
  W <- model.matrix(~ -1 + Treatments, PBIBD2.lay)
  Vu <- with(PBIBD2.lay, fac.vcmat(Blocks, 5))
  R <- diag(1, nrow(PBIBD2.lay))
  
  #'## Calculate the variance matrix of the predictions
  Vp <- mat.Vpred(W = W, Vu=Vu, R=R)
  testthat::expect_true(abs(designAmeasures(Vp) - 0.5625) < 1e-6)

  #Test when no random terms
  Vpr <- mat.Vpredicts(target = ~ -1 + Treatments, fixed = ~ Blocks -1, design = PBIBD2.lay)
  testthat::expect_true(abs(designAmeasures(Vpr) - 0.5666667) < 1e-6)
  
    
  #'## Get reduced variance matix of predictions
  unit.str <- pstructure(~Blocks/Units, grandMean = TRUE, data = PBIBD2.lay)
  Q.B <- projector(unit.str$Q$Mean + unit.str$Q$Blocks)
  testthat::expect_equal(degfree(Q.B), 6)
  Vp.reduc <- mat.Vpred(W = W, Vu=Vu, R=R, eliminate = Q.B)
  testthat::expect_true(all(abs(rowSums(Vp.reduc)) < 1e-10))
  testthat::expect_true(abs(designAmeasures(Vp.reduc) - 
                              2/(4*PBIBD2.canon$Q[[1]]$`Units[Blocks]`$Treatments$adjusted$aefficiency)) 
                        < 1e-10)
})

cat("#### Test for mat.Vpredicts using using Smith's small example\n")
test_that("SmithSmallVpredicts", {
  skip_on_cran()
  library(dae)
  #'# Script to investigate the reduced design from Smith et al. (2015)
  
  #'## Set up designs
  mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
  field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
  field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                              levels = c("Y","W","G","M","D","E"))
  start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
  rownames(start.design) <- NULL
  start.design
  
  #'## Set gammas
  terms <- c("Variety", "Frep", "Frep:Fplot", "Mrep", "Mrep:Mday", "Mrep:Mday:Mord")
  gammas <- c(1, 0.1, 0.2, 0.3, 0.2, 1)
  names(gammas) <- terms
  print(gammas)
  
  #'## Set up matrices for fixed Variety effects
  W <- model.matrix(~ -1 + Variety, start.design)
  ZMr <- model.matrix(~ -1 + Mrep, start.design)
  ZMrMd <- model.matrix(~ -1 + Mrep:Mday, start.design)
  ZFr <- model.matrix(~ -1 + Frep, start.design)
  ZFrFp <- model.matrix(~ -1 + Frep:Fplot, start.design)
  t <- length(levels(start.design$Variety))
  ng <- ncol(W)
  Vu <- gammas["Mrep"] * ZMr %*% t(ZMr) + gammas["Mrep:Mday"] * ZMrMd %*% t(ZMrMd) + 
    gammas["Frep"] * ZFr %*% t(ZFr) + gammas["Frep:Fplot"] * ZFrFp %*% t(ZFrFp)
  R <- diag(1, nrow(start.design))
  
  #'# Calculate variance matrix
  Vp <- mat.Vpred(W = W, Vu=Vu, R=R)
 
  ## Specify matrices to calculate the variance matrix of the predicted fixed Variety effects 
  W <- model.matrix(~ -1 + Variety, start.design)
  Vu <- with(start.design, fac.vcmat(Mrep, gammas["Mrep"]) + 
               fac.vcmat(fac.combine(list(Mrep,Mday)), gammas["Mrep:Mday"]) + 
               fac.vcmat(Frep, gammas["Frep"]) + 
               fac.vcmat(fac.combine(list(Frep,Fplot)), gammas["Frep:Fplot"]))
  R <- diag(1, nrow(start.design))
  
  ## Calculate variance matrix
  Vp <- mat.Vpredicts(target = W, random = Vu, R = R)
  testthat::expect_equal(attr(Vp, which = "rank"), 5)
             
  ## Calculate information matrix
  C <- mat.Vpredicts(target = W, random = Vu, R = R, result = "inform")
  testthat::expect_equal(attr(C, which = "rank"), 5)
  
  #'## Calculate A measures
  AVpred <- designAmeasures(Vp)
  testthat::expect_true(all(abs(AVpred - 1.521905) < 1e-6))
  
  Vp <- mat.Vpredicts(target = W, random = Vu, R = R)
  AVpredict <- designAmeasures(Vp)
  testthat::expect_true(all(abs(AVpredict - AVpred) < 1e-6))
  
  Vp <- mat.Vpredicts(target = ~ -1 + Variety, 
                      fixed = ~ 1, 
                      random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                      G = as.list(gammas[c(4,5,2,3)]), 
                      R = R, design = start.design)
  Aform <- designAmeasures(Vp)
  testthat::expect_true(all(abs(Aform - AVpred) < 1e-6))
  
  ## Calculate the variance matrix of the predicted random Variety effects using formulae
  Vp <- mat.Vpredicts(target = ~ -1 + Variety, Gt = 1, 
                      fixed = ~ 1, 
                      random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                      G = as.list(gammas[c(4,5,2,3)]), 
                      R = R, design = start.design)
  Aran <- designAmeasures(Vp)
  testthat::expect_true(all(abs(Aran - 0.8564816) < 1e-6))
  
  
  ## Calculate the variance matrix of the predicted fixed Variety effects, 
  ## elminating the grand mean
  n <- nrow(start.design)
  Vp.reduc <- mat.Vpredicts(target = ~ -1 + Variety, 
                            random = ~ -1 + Mrep/Mday + Frep/Fplot, 
                            G = as.list(gammas[c(4,5,2,3)]), 
                            eliminate = projector(matrix(1, nrow = n, ncol = n)/n), 
                            design = start.design)
  testthat::expect_true(all(abs(rowSums(Vp.reduc)) < 1e-10))
  testthat::expect_true(abs(designAmeasures(Vp.reduc) - 
                              sum(diag(Vp.reduc)) * 2 / (nrow(Vp.reduc) - 1)) < 1.e-10)
})  

