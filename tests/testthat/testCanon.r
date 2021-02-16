#devtools::test("dae")
context("canonical")

cat("#### Test for  using Cochran&Cox PBIBD2\n")
test_that("PBIBD2", {
  skip_on_cran()
  library(dae)
  #'# PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 2nd edn Wiley, New York"
  
  #'## Input the design and randomize"
  Treatments <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
  PBIBD2.lay <- designRandomize(allocated = Treatments,
                                recipient = list(Blocks =6, Units = 4),
                                nested.recipients = list(Units = "Blocks"),
                                seed = 98177)
  
  #'## Show that projs.canon is deprecated
  testthat::expect_warning(PBIBD2.canon <- 
                             projs.canon(formulae = list(unit = ~Blocks/Units,
                                                         trt = ~ Treatments),
                                         which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                         labels = "terms", data = PBIBD2.lay))
})

cat("#### Test for designAnatomy using Cochran&Cox PBIBD2\n")
test_that("PBIBD2", {
  skip_on_cran()
  library(dae)
  #'# PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 2nd edn Wiley, New York"

  #'## Input the design and randomize"
  Treatments <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
  PBIBD2.lay <- designRandomize(allocated = Treatments,
                                recipient = list(Blocks =6, Units = 4),
                                nested.recipients = list(Units = "Blocks"),
                                seed = 98177)
  
  unit.struct <- pstructure(formula = ~Blocks/Units, labels = "terms", data = PBIBD2.lay)
  testthat::expect_equal(degfree(unit.struct$Q$Blocks), 5)
  testthat::expect_equal(degfree(unit.struct$Q$'Blocks:Units'), 18)
  trt.struct <- pstructure(formula = ~ Treatments, labels = "terms", data = PBIBD2.lay)
  testthat::expect_equal(degfree(trt.struct$Q$Treatments), 5)
  
  #Test for single structure
  unit.canon <- designAnatomy(formulae = list(units =  ~ Blocks/Units), data = PBIBD2.lay)
  testthat::expect_equal(unit.canon$Q[[1]]$`Units[Blocks]`, 18)
  summ <- summary(unit.canon)
  testthat::expect_equal(summ$decomp$` df`[1], 5)
  unit.canon <- designAnatomy(formulae = list(units =  ~ Blocks/Units), 
                              omit.projectors = "none", data = PBIBD2.lay)
  testthat::expect_equal(degfree(unit.canon$Q[[1]]$`Units[Blocks]`), 18)
  summ <- summary(unit.canon)
  testthat::expect_equal(summ$decomp$` df`[2], 18)
  
  #'## Compute the anatomy
  PBIBD2.canon <- designAnatomy(formulae = list(unit = ~Blocks/Units,
                                                trt = ~ Treatments),
                                which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                labels = "terms", data = PBIBD2.lay)
  testthat::expect_equal(PBIBD2.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2.canon$Q[[1]]$'Blocks:Units'$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  testthat::expect_equal(PBIBD2.canon$Q[[2]]$'Blocks&Treatments', 2)
  testthat::expect_equal(PBIBD2.canon$Q[[2]]$'Blocks&Residual', 3)
  testthat::expect_equal(PBIBD2.canon$Q[[2]]$'Blocks:Units&Treatments', 5)
  testthat::expect_equal(PBIBD2.canon$Q[[2]]$'Blocks:Units&Residual', 13)
  
  #'## unique plot levels
  PBIBD2.lay$AUnits <- with(PBIBD2.lay, fac.combine(list(Blocks,Units)))
  PBIBD2A.canon <- designAnatomy(list(unit = ~Blocks + AUnits,
                                      trt = ~ Treatments),
                                 labels = "terms", data = PBIBD2.lay,
                                 which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                 orthogonalize = "eigen")
  summary(PBIBD2A.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(PBIBD2A.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2A.canon$Q[[1]]$AUnits$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Treatments', 2)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Residual', 3)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'AUnits&Treatments', 5)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'AUnits&Residual', 13)
  
  
  unit.struct <- pstructure(formula = ~Blocks + AUnits, data = PBIBD2.lay,
                            labels = "terms", orthogonalize = "hybrid")
  testthat::expect_equal(degfree(unit.struct$Q$Blocks), 5)
  testthat::expect_equal(degfree(unit.struct$Q$AUnits), 18)

  PBIBD2A.canon <- designAnatomy(formulae = list(unit = ~Blocks + AUnits,
                                                 trt = ~ Treatments),
                                 labels = "terms", data = PBIBD2.lay, 
                                 which.criteria = c('aeff', 'xeff', 'eeff','order'))
  summary(PBIBD2A.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(PBIBD2A.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2A.canon$Q[[1]]$AUnits$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Treatments', 2)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Residual', 3)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'AUnits&Treatments', 5)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'AUnits&Residual', 13)
  
  PBIBD2A.canon <- designAnatomy(formulae = list(unit = ~Blocks/AUnits,
                                                 trt = ~ Treatments),
                                 labels = "terms", data = PBIBD2.lay,
                                 which.criteria = c('aeff', 'xeff', 'eeff','order'))
  summary(PBIBD2A.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(PBIBD2A.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_equal(PBIBD2A.canon$Q[[1]]$Blocks$Treatments$adjusted$xefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2A.canon$Q[[1]]$'Blocks:AUnits'$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Treatments', 2)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks&Residual', 3)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks:AUnits&Treatments', 5)
  testthat::expect_equal(PBIBD2A.canon$Q[[2]]$'Blocks:AUnits&Residual', 13)
  
  #Some tests of subsidiary functions
  ## manually obtain projectors for units
  Q.G <- projector(matrix(1, nrow=24, ncol=24)/24)                         
  Q.B <- projector(fac.meanop(PBIBD2.lay$Block) - Q.G)
  Q.BP <- projector(diag(1, nrow=24) - Q.B - Q.G)
  testthat::expect_equal(degfree(Q.B), 5)
  testthat::expect_equal(degfree(Q.BP), 18)
  
  ## manually obtain projector for trt
  Q.T <- projector(fac.meanop(PBIBD2.lay$Treatments) - Q.G)
  testthat::expect_equal(degfree(Q.T), 5)
  
  ##compute intrablock efficiency criteria
  effic <- proj2.efficiency(Q.BP, Q.T)
  effic
  effCrit <- efficiency.criteria(effic)
  testthat::expect_lt(abs(effCrit$aefficiency - 0.8823529), 1e-05)
  testthat::expect_lt(abs(effCrit$mefficiency - 0.9), 1e-05)
  testthat::expect_lt(abs(effCrit$sefficiency - 0.01875), 1e-05)
  testthat::expect_lt(abs(effCrit$xefficiency - 1), 1e-05)
  testthat::expect_equal(effCrit$order, 2)
  testthat::expect_equal(effCrit$dforthog, 3)
  
  ##obtain combined decomposition and summarize
  eff <- efficiencies.pcanon(PBIBD2A.canon)[[1]]
  testthat::expect_equal(length(eff$Blocks$Treatments), 2)
  testthat::expect_equal(length(eff$'Blocks:AUnits'$Treatments), 5)
  
})

cat("#### Test for BIBD example\n")
test_that("BIBD", {
  skip_on_cran()
  library(dae)
  "BIBD example"
  BIBD.lay <- fac.gen(list(run = 15, location = 4))
  BIBD.lay$trts <- factor(c(1,2,3,4, 1,2,5,6, 1,3,7,8, 1,4,9,10, 1,5,7,9, 
                            1,6,8,10, 2,3,6,9, 2,4,7,10, 2,5,8,10, 2,7,8,9, 
                            3,5,9,10, 3,6,7,10, 3,4,5,8, 4,5,6,7, 4,6,8,9))
  set.daeTolerance(1E-06)
  loc.trts.canon <- designAnatomy(formulae = list(loc=~ run/location, trts=~ trts), 
                                  data = BIBD.lay, labels = "terms", 
                                  omit.projectors="pc")
  summary(loc.trts.canon, which.criteria =c("aefficiency", "order"))
  testthat::expect_equal(length(loc.trts.canon$Q[[2]]), 4)
  set.daeTolerance(1E-08)
  
  LcombT <- loc.trts.canon$Q[[2]]
  testthat::expect_equal(degfree(LcombT[['run&Residual']]), 5)
  Q <- projector(LcombT[['run&Residual']])
  testthat::expect_equal(degfree(Q), 5)
  
  ## manually obtain projectors for units
  Q.G <- projector(matrix(1, nrow=60, ncol=60)/60)                         
  Q.r <- projector(fac.meanop(BIBD.lay$run) - Q.G)
  Q.rl <- projector(diag(1, nrow=60) - Q.r - Q.G)
  
  ## manually obtain projector for trt
  Q.t <- projector(fac.meanop(BIBD.lay$trts) - Q.G)
  
  ##compute intrablock efficiency criteria
  effic <- proj2.efficiency(Q.rl, Q.t)
  testthat::expect_equal(length(effic), 9)
  testthat::expect_lt(abs(unique(signif(effic, 6)) - 0.8333333), 1e-05)
  
  Q.r.t <- proj2.combine(Q.r, Q.t)$Qres
  testthat::expect_equal(degfree(Q.r.t), 5)
  Q.r.t <- Q.r %*% Q.t %*% Q.r/(1/6)
  Q.r.t <- Q.r - Q.r.t
  P <- projector(Q.r.t)
  testthat::expect_equal(degfree(P), 5)
  
})

cat("#### Test for Jarrett & Ruggiero example\n")
test_that("JarrettRuggiero", {
  skip_on_cran()
  library(dae)
  #'## Jarrett & Ruggiero example
  jr.lay <- fac.gen(list(Set=7, Dye=2,Array=3))
  jr.lay <- within(jr.lay, 
                   { Block <- factor(rep(1:7, each=6))
                   Plant <- factor(rep(c(1,2,3,2,3,1), times=7))
                   Sample <- factor(c(rep(c(2,1,2,2,1,1, 1,2,1,1,2,2), times=3), 2,1,2,2,1,1))
                   S1 <- Dye
                   Treat <- factor(c(1,2,4,2,4,1, 2,3,5,3,5,2, 3,4,6,4,6,3, 4,5,7,5,7,4, 
                                     5,6,1,6,1,5, 6,7,2,7,2,6, 7,1,3,1,3,7),
                                   labels=c("A","B","C","D","E","F","G"))
                   })
  
  array.plot.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye, 
                                                    plot = ~ Block/Plant),
                                    labels = "terms", data = jr.lay, 
                                    omit.projectors="none")
  testthat::expect_equal(length(array.plot.canon$Q[[2]]), 5)
  testthat::expect_true(class(array.plot.canon$Q[[2]]$Dye) == "projector")
  
  summ.default <- summary(array.plot.canon)
  testthat::expect_equal(nrow(summ.default$decomp), 5)
  testthat::expect_equal(ncol(summ.default$decomp), 7)
  summ.none <- summary(array.plot.canon, which = "none")
  testthat::expect_equal(ncol(summ.none$decomp), 4)
  summ.3 <- summary(array.plot.canon, which = c("aeff", "eeff", "xeff"))
  testthat::expect_equal(ncol(summ.3$decomp), 7)
  testthat::expect_true("xefficiency" %in% names(summ.3$decomp))
  
  
  array.plot.trt.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye,
                                                        plot = ~ Block/Plant/Sample,
                                                        trt = ~ Treat),
                                        labels = "terms", data = jr.lay, 
                                        omit.projectors="comb")
  
  testthat::expect_equal(length(array.plot.trt.canon$Q[[3]]), 7)
  testthat::expect_equal(array.plot.trt.canon$Q[[3]]$`Set:Array:Dye&Block:Plant:Sample`, 6)
  testthat::expect_lt(abs(array.plot.trt.canon$Q[[2]]$`Set:Array:Dye&Block:Plant`$Treat$adjusted$aefficiency - 0.58333333), 1e-05)
  testthat::expect_lt(abs(array.plot.trt.canon$Q[[1]]$`Set:Array:Dye`$`Block:Plant`$adjusted$aefficiency - 0.75), 1e-05)
  testthat::expect_true(array.plot.trt.canon$aliasing$Source == "Block:Plant:Sample")
  testthat::expect_true(array.plot.trt.canon$aliasing$Alias == "Block:Plant")
  testthat::expect_true(array.plot.trt.canon$aliasing$df == 14)
  testthat::expect_true(array.plot.trt.canon$aliasing$aefficiency == 1)
  
  summ.default <-summary(array.plot.trt.canon)
  testthat::expect_equal(nrow(summ.default$decomp), 7)
  testthat::expect_equal(ncol(summ.default$decomp), 9)
  summ.none <-summary(array.plot.trt.canon, which.criteria = "none")
  testthat::expect_equal(ncol(summ.none$decomp), 6)
  summ.all <-summary(array.plot.trt.canon, which.criteria = "all")
  testthat::expect_equal(ncol(summ.all$decomp), 13)
  summ.2 <-summary(array.plot.trt.canon, which.criteria =c("aefficiency", "order"))
  testthat::expect_equal(ncol(summ.2$decomp), 8)
  
  
  eff.adj <- efficiencies.pcanon(array.plot.trt.canon)
  testthat::expect_equal(length(eff.adj[[1]]), 3)
  testthat::expect_equal(unique(signif(eff.adj[[1]]$`Set:Array:Dye`$`Block:Plant`, 2)), 0.75)
  testthat::expect_equal(length(eff.adj[[1]]$`Set:Array:Dye`$`Block:Plant:Sample`), 6)
  testthat::expect_equal(length(eff.adj[[2]]), 3)
  testthat::expect_equal(length(eff.adj[[2]]$`Set:Array:Dye&Block:Plant`$Treat), 6)
  eff.pair <- efficiencies.pcanon(array.plot.trt.canon, which = "pairwise")
  testthat::expect_equal(length(eff.pair[[1]]), 3)
  testthat::expect_equal(unique(signif(eff.pair[[1]]$`Set:Array:Dye`$`Block:Plant`, 2)), 0.75)
  testthat::expect_equal(length(eff.pair[[1]]$`Set:Array:Dye`$`Block:Plant:Sample`), 20)
  testthat::expect_equal(length(eff.pair[[2]]), 3)
  testthat::expect_equal(length(eff.pair[[2]]$`Set:Array:Dye&Block:Plant`$Treat), 6)
  
})

cat("#### Test for designTwophaseAnatomies\n")
test_that("designTwophaseAnatomies", {
  skip_on_cran()
  library(dae)
  #'## Jarrett & Ruggiero example
  jr.lay <- fac.gen(list(Set=7, Dye=2,Array=3))
  jr.lay <- within(jr.lay, 
                   { 
                     Block <- factor(rep(1:7, each=6))
                     Plant <- factor(rep(c(1,2,3,2,3,1), times=7))
                     Sample <- factor(c(rep(c(2,1,2,2,1,1, 1,2,1,1,2,2), times=3), 2,1,2,2,1,1))
                     Treat <- factor(c(1,2,4,2,4,1, 2,3,5,3,5,2, 3,4,6,4,6,3, 4,5,7,5,7,4, 
                                       5,6,1,6,1,5, 6,7,2,7,2,6, 7,1,3,1,3,7),
                                     labels=c("A","B","C","D","E","F","G"))
                   })
  
  testthat::expect_error(jr.anat <- designTwophaseAnatomies(
    formulae = list(array = ~ (Set:Array)*Dye, 
                    plot = ~ Block/Plant),
    data = jr.lay))
  jr.anat <- designTwophaseAnatomies(formulae = list(array = ~ (Set:Array)*Dye,
                                                     plot = ~ Block/Plant/Sample,
                                                     trt = ~ Treat),
                                     data = jr.lay)  
  testthat::expect_equal(length(jr.anat), 4)
  testthat::expect_true(all(names(jr.anat) == c("twophase","first","cross","units")))
  testthat::expect_true(all(attr(jr.anat, which = "titles") == 
                              c("Anatomy for the full two-phase design",
                                "Anatomy for the first-phase design", 
                                "Anatomy for the cross-phase, treatments design", 
                                "Anatomy for the combined-units design" )))
  
  #Test for getting only two anatomies
  jrfc.anat <- designTwophaseAnatomies(formulae = list(array = ~ (Set:Array)*Dye,
                                                       plot = ~ Block/Plant/Sample,
                                                       trt = ~ Treat),
                                       which.designs = c("first", "cross"), data = jr.lay)
  testthat::expect_equal(length(jrfc.anat), 4)
  testthat::expect_true(all(names(jrfc.anat) == c("twophase","first","cross","units")))
  testthat::expect_true(is.null(jrfc.anat[["twophase"]]))
  testthat::expect_true(!is.null(jrfc.anat[["first"]]))
  testthat::expect_true(!is.null(jrfc.anat[["cross"]]))
  testthat::expect_true(is.null(jrfc.anat[["units"]]))
  testthat::expect_true(all(attr(jrfc.anat, which = "titles") == 
                              c("Anatomy for the full two-phase design",
                                "Anatomy for the first-phase design", 
                                "Anatomy for the cross-phase, treatments design", 
                                "Anatomy for the combined-units design" )))
  
  testthat::expect_silent(jrfc.anat <- designTwophaseAnatomies(formulae = list(array = ~ (Set:Array)*Dye,
                                                                               plot = ~ Block/Plant/Sample,
                                                                               trt = ~ Treat),
                                                               which.designs = c("first", "cross"), 
                                                               printAnatomies = FALSE, data = jr.lay))
})

cat("#### Test for designTwophaseAnatTitles\n")
test_that("designTwophaseAnatTitles", {
  skip_on_cran()
  library(dae)
  options(width = 95)
  
  # Generate first-phase sytematic design
  ph1.sys <- cbind(fac.gen(list(Expressive = c("Yes", "No"), Patients = 4, Occasions = 2)),
                   fac.gen(list(Motions = c("active", "passive")), times = 8))
  
  # Generate the two-phase systematic design
  ph2.sys <- cbind(fac.gen(list(Raters = 74, Viewings = 16)),
                   fac.gen(list(Trainings = 2, 16), times = 37),
                   rep(ph1.sys, times =74))
  
  #'## Randomize the two-phase design
  ph2.lay <- designRandomize(allocated = ph2.sys[c("Trainings", "Expressive", "Patients",
                                                   "Occasions", "Motions")],
                             recipient = ph2.sys[c("Raters", "Viewings")],
                             except = "Viewings",
                             seed = 15674)
  testthat::expect_equal(nrow(ph2.lay), 1184)
  testthat::expect_equal(ncol(ph2.lay), 7)
  
  # Produce the anatomies of the design and check that titles are correct
  jfh.anat <- designTwophaseAnatomies(formulae = list(rate  = ~ Raters * Viewings,
                                                     video = ~ (Expressive/Patients)*Occasions,
                                                     alloc = ~ Trainings * Motions),
                                     titles = c(NA, NA,
                                               "Anatomy for allocated factors with second-phase units",
                                                NA),
                                     data = ph2.lay)  
  testthat::expect_equal(length(jfh.anat), 4)
  testthat::expect_true(all(names(jfh.anat) == c("twophase","first","cross","units")))
  testthat::expect_true(all(attr(jfh.anat, which = "titles") == 
                              c("Anatomy for the full two-phase design",
                                "Anatomy for the first-phase design", 
                                "Anatomy for allocated factors with second-phase units", 
                                "Anatomy for the combined-units design" )))
})

cat("#### Test for Baby pseudoterm example\n")
test_that("Baby", {
  skip_on_cran()
  library(dae)
  #'## Baby pseudoterm example
  pseudo.lay <- data.frame(pl = factor(1:12),
                           ab = factor(rep(1:4, times=3)),
                           a = factor(rep(1:2, times=6)))

  trt.unit <- pstructure(formula = ~ ab+a, labels = "terms", data = pseudo.lay)
  testthat::expect_equal(length(trt.unit$Q), 1)
  testthat::expect_equal(degfree(trt.unit$Q$ab), 3)
  
  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~a+ab), pseudo.lay, 
                                labels = "terms")
  summary(pseudo.canon)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&a', 1)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab', 2)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&Residual', 8)
  
  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~ab+a), pseudo.lay,
                                labels = "terms")
  summary(pseudo.canon)
  testthat::expect_true(all(pseudo.canon$aliasing$Source == c("a", "a")))
  testthat::expect_true(all(pseudo.canon$aliasing$Alias == c("ab", "## Aliased")))
  testthat::expect_true(all(pseudo.canon$aliasing$aefficiency == c(1, 0)))
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab', 3)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&Residual', 8)

  
  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~ab+a), pseudo.lay,
                                labels = "terms", orhtogonalize = "eigen")
  summary(pseudo.canon)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab', 3)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&Residual', 8)
  
})

cat("#### Test for pseudoreplicated N experiment\n")
test_that("pseudoN", {
  skip_on_cran()
  library(dae)

  #'## A design with the same RCBD in each Area
  sameRCBD.lay <- designRandomize(allocated = fac.gen(list(Nitrogen = 3, 3, Varieties = 24)),
                                  recipient = list(Areas = 3, Blocks = 3, Plots = 24),
                                  nested.recipients = list(Plots = "Blocks"),
                                  seed = 6411)
  sameRCBD.lay <- within(sameRCBD.lay, 
                         {
                           Plot <- fac.combine(list(Blocks,Plots))
                         })
  testthat::expect_equal(nrow(sameRCBD.lay), 216)
  testthat::expect_equal(ncol(sameRCBD.lay), 6)
  
  #'### nested numbering of levels
  sameRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas*(Blocks/Plots),
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "terms", data = sameRCBD.lay, 
                                  grandMean = TRUE)
  summary(sameRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(sameRCBD.canon$Q[[2]]) == c(1,2,2,23,46,4,46,92)))
  #'### unique numbering of levels
  sameRCBD.canon <- designAnatomy(formulae = list(unit = ~  Areas*(Blocks + Plot),
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "terms", data = sameRCBD.lay) 
  summary(sameRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(sameRCBD.canon$Q[[2]]) == c(2,2,23,46,4,46,92)))
  
  #'## A design with different RCBDs in each Area
  diffRCBD.lay <- designRandomize(allocated = fac.gen(list(Nitrogen = 3, 3, Varieties = 24)),
                                  recipient = list(Areas = 3, Blocks = 3, Plots = 24),
                                  nested.recipients = list(Plots = "Blocks"),
                                  seed = 6467)
  diffRCBD.lay <- within(diffRCBD.lay, 
                         {
                           Block <- fac.combine(list(Areas,Blocks))
                           Plot <- fac.combine(list(Areas,Blocks,Plots))
                         })
  testthat::expect_equal(nrow(diffRCBD.lay), 216)
  testthat::expect_equal(ncol(diffRCBD.lay), 7)
  
  #'### nested numbering of levels
  diffRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas/Blocks/Plots,
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "terms", data = diffRCBD.lay, 
                                  grandMean = TRUE)
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(1,2,6,23,46,138)))
  #'### unique numbering of levels
  diffRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas + Block + Plot,
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "terms", data = diffRCBD.lay)
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(2,6,23,46,138)))
  #'### Use eigen and diff instead of hybrid
  diffRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas + Block + Plot,
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "terms", data = diffRCBD.lay, 
                                  grandMean = TRUE, 
                                  orthogonalize = c("eigen", "diff"))
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(1,2,6,23,46,138)))
  
})


cat("#### Test for Repeated LSD for Housewives example\n")
test_that("Housewives", {
  skip_on_cran()
  library(dae)
  data("LSRepeatHwife.dat")
  
  Hwife.canon <- designAnatomy(formulae = list(units = ~ Month:Week + Month:Hwife + 
                                                          Month:Week:Hwife),
                               labels = "terms", data = LSRepeatHwife.dat)
  summary(Hwife.canon)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Month:Week`, 7)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Month:Hwife`, 6)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Month:Week:Hwife`, 18)
})

cat("#### Test for Preece examples\n")
test_that("Preece", {
  skip_on_cran()
  library(dae)
  #'## Preece examples with two sets of treatments in a BIBD
  #'### two tiers
  preece1.lay <- fac.gen(list(block=10, plot=3))
  preece1.lay$T1 <- factor(c(1,3,4, 2,4,5, 3,5,1, 4,1,2, 5,2,3, 
                             1,2,5, 2,3,1, 3,4,2, 4,5,3, 5,1,4))
  preece1.lay$T2 <- factor(c(1,2,5, 2,3,1, 3,4,2, 4,5,3, 5,1,4, 
                             1,4,3, 2,5,4, 3,1,5, 4,2,1, 5,3,2))
  preece1.lay$T3 <- factor(c(1,4,3, 2,5,4, 3,1,5, 4,2,1, 5,3,2, 
                             1,5,2, 2,1,3, 3,2,4, 4,3,5, 5,4,1))

  trt.struct <- pstructure(formula = ~ T1+T2, labels = "terms", data = preece1.lay)
  testthat::expect_equal(length(trt.struct$Q), 2)
  testthat::expect_equal(degfree(trt.struct$Q$T1), 4)
  testthat::expect_equal(degfree(trt.struct$Q$T2), 4)
  
  preece1.canon <- designAnatomy(formulae = list(plot= ~ block/plot, trt= ~ T1+T2), 
                                 labels = "terms", data = preece1.lay)
  summary(preece1.canon)
  testthat::expect_equal(length(preece1.canon$Q[[1]]), 2)
  testthat::expect_equal(length(preece1.canon$Q[[1]]$block), 3)
  testthat::expect_equal(length(preece1.canon$Q[[1]]$`block:plot`), 3)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$block$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$block$T2$adjusted$aefficiency - 0.0952), 1e-04)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$`block:plot`$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$`block:plot`$T2$adjusted$aefficiency - 0.7619), 1e-04)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&T1', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&T2', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&Residual', 1)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block:plot&T1', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block:plot&T2', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block:plot&Residual', 12)
  
  preece2.canon <- designAnatomy(formulae = list(plot= ~ block/plot, trt= ~ T1+T3),
                                 labels = "terms", data = preece1.lay)
  summary(preece2.canon)
  testthat::expect_equal(length(preece2.canon$Q[[1]]), 2)
  testthat::expect_equal(length(preece2.canon$Q[[1]]$block), 2)
  testthat::expect_equal(length(preece2.canon$Q[[1]]$`block:plot`), 3)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$block$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$`block:plot`$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$`block:plot`$T3$adjusted$aefficiency - 0.8571), 1e-04)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block&T1', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block&Residual', 5)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block:plot&T1', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block:plot&T3', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block:plot&Residual', 12)
  
  #'## three tiers
  preece3.recip <- list(occasions=10, analysis=3)
  preece3.nest <- list(analysis = "occasions")
  preece3.tier.lay <- designRandomize(allocated = preece1.lay, 
                                      recipient = preece3.recip, 
                                      nested.recipients = preece3.nest, 
                                      seed = 80851)
  preece3.canon <- designAnatomy(formulae = list(lab= ~ occasions/analysis, 
                                                 plot= ~ block/plot, trt= ~ T1+T2),
                                 labels = "terms", data = preece3.tier.lay)
  summary(preece3.canon)
  testthat::expect_equal(length(preece3.canon$Q[[2]]), 2)
  testthat::expect_equal(length(preece3.canon$Q[[2]]$'occasions&block'), 3)
  testthat::expect_equal(length(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'), 3)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions&block'$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions&block'$T2$adjusted$aefficiency - 0.0952), 1e-04)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'$T2$adjusted$aefficiency - 0.7619), 1e-04)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions&block&T1', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions&block&T2', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions&block&Residual', 1)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&T1', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&T2', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&Residual', 12)
  
  preece3.canon <- designAnatomy(formulae = list(lab= ~ occasions/analysis, 
                                                 plot= ~ block/plot, trt= ~ T1+T3),
                                 labels = "terms", data = preece3.tier.lay)
  summary(preece3.canon)
  testthat::expect_equal(length(preece3.canon$Q[[2]]), 2)
  testthat::expect_equal(length(preece3.canon$Q[[2]]$'occasions&block'), 2)
  testthat::expect_equal(length(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'), 3)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions&block'$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece3.canon$Q[[2]]$'occasions:analysis&block:plot'$T3$adjusted$aefficiency - 0.8571), 1e-04)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions&block&T1', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions&block&Residual', 5)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&T1', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&T3', 4)
  testthat::expect_equal(preece3.canon$Q[[3]]$'occasions:analysis&block:plot&Residual', 12)
  testthat::expect_true(all(preece3.canon$aliasing$Source == "T3"))
  testthat::expect_true(all(preece3.canon$aliasing$df == c(4,4,0,4)))
  testthat::expect_true(all(preece3.canon$aliasing$Alias == c("T1", "## Information remaining",
                                                              "## Aliased", "T1")))
  testthat::expect_true(all(preece3.canon$aliasing$In == c("trt", "trt", "occasions&block",              
                                                           "occasions:analysis&block:plot")))
  testthat::expect_true(all(abs(preece3.canon$aliasing$aefficiency - c(0.02777778, 0.97222222, 
                                                                       0.00000000, 0.02702703)) < 1e-05))
  
  
})

cat("#### Test for Mostafa's green wall experiment in 2014\n")
test_that("Mostafa", {
  skip_on_cran()
  library(dae)
  #Mostafa's green wall experiment in 2014
  data(gwall.lay)
  options(width = 100, nwarnings = 150)
  set.daeTolerance(1e-06,1e-06)
  pot.treat.canon <- designAnatomy(formulae = list(pot = ~ Rows*Cols,
                                                   trt = ~ Species*Irrigation*Media + 
                                                     First/(SpeCarry*IrrCarry*MedCarry)), 
                                   labels = "terms", data = gwall.lay, 
                                   keep.order=TRUE)
  summary(pot.treat.canon)
  testthat::expect_equal(length(pot.treat.canon$Q[[2]]), 17)
  testthat::expect_equal(pot.treat.canon$Q[[2]]$`Rows:Cols&Residual`, 69)
  testthat::expect_lt(abs(pot.treat.canon$Q[[1]]$`Rows:Cols`$`Species:Irrigation:Media`$adjusted$aefficiency - 0.8240828), 1e-05)
  testthat::expect_lt(abs(pot.treat.canon$Q[[1]]$`Rows:Cols`$`First:SpeCarry:IrrCarry`$adjusted$aefficiency - 0.8320062), 1e-05)
  
  trt.struct <- pstructure(formula = ~ Species*Irrigation*Media + 
                           First/(SpeCarry*IrrCarry*MedCarry), 
                         labels = "terms", data = gwall.lay)
  testthat::expect_equal(length(trt.struct$Q), 11)
  testthat::expect_equal(degfree(trt.struct$Q$`Species:Irrigation:Media`), 8)
  testthat::expect_equal(degfree(trt.struct$Q$`First:SpeCarry:IrrCarry`), 8)
  
})

cat("#### Test for four-tiered corn example\n")
test_that("corn", {
  skip_on_cran()
  library(dae)
  #'## Four-tiered corn experiment with 3 structures including pseudoterms
  #'### Randomize field factors to lab factors
  corn1.recip <- list(Intervals=18, ConPlate=36)
  corn1.nest <- list(ConPlate = "Intervals")
  corn1.alloc <- fac.gen(list(Sites=3, Blocks=2, Plots=3, Lots=36))
  #corn1.alloc$Harvesters <- factor(1:3, each=36, times=6)
  corn1.lay <- designRandomize(allocated = corn1.alloc, 
                               recipient = corn1.recip, 
                               nested.recipients = corn1.nest, 
                               seed = 81505)
  
  #'### Randomize treatments to lab factors
  corn2.recip <- list(Intervals=18, Containers=9, Plates=4)
  corn2.nest <- list(Containers = "Intervals", Plates = c("Containers", "Intervals"))
  Treats <- factor(rep(1:9, each=4, times=18))
  corn2.lay <- designRandomize(allocated = Treats, 
                               recipient = corn2.recip, 
                               nested.recipients = corn2.nest, 
                               seed = 543205)
  corn2.lay <- corn2.lay[-1]
  
  #'### Randomize Harvesters to field factors
  corn3.recip <- list(Sites=3, Blocks=2, Plots=3, Lots=36)
  corn3.nest <- list(Blocks = "Sites", Plots = c("Blocks", "Sites"), Lots = c("Plots", "Blocks", "Sites"))
  Harvesters <- factor(rep(1:3, each=36, times=6))
  corn3.lay <- designRandomize(allocated = Harvesters, 
                               recipient = corn3.recip, 
                               nested.recipients = corn3.nest, 
                               seed = 135205)
  
  #'### Combine randomizations
  corn.lay <- cbind(corn1.lay, corn2.lay)
  corn.lay <- merge(corn.lay, corn3.lay)
  corn.lay <- corn.lay[c("Intervals", "Containers", "Plates", "ConPlate", 
                         "Sites", "Blocks", "Plots", "Lots", 
                         "Treats", "Harvesters")]
  corn.lay <- corn.lay[do.call(order, corn.lay), ]
  rownames(corn.lay) <- NULL
  
  #'### Make factors with unique levels
  corn.lay <- within(corn.lay, 
                     {
                       AContainers <- fac.combine(list(Intervals, Containers))
                       APlates <- fac.combine(list(AContainers, Plates))
                       ABlocks <- fac.combine(list(Sites, Blocks))
                       APlots <- fac.combine(list(ABlocks, Plots))
                       ALots <- fac.combine(list(APlots, Lots))
                     })
  
  #'## Check properties
  corn.canon <- designAnatomy(formulae = list(plate= ~ Intervals + AContainers + APlates, 
                                              field= ~ Sites + ABlocks + APlots + ALots,
                                              trts= ~ Sites*Harvesters*Treats),
                              labels = "terms", data = corn.lay)
  summary(corn.canon, which.criteria="aeff")
  testthat::expect_equal(corn.canon$Q[[3]]$'Intervals&APlots&Residual', 6)
  testthat::expect_equal(corn.canon$Q[[3]]$'AContainers&ALots&Residual', 72)
  testthat::expect_equal(corn.canon$Q[[3]]$'APlates&ALots', 486)
  
  #Create an example with double Residual
  corn2Res.canon <- designAnatomy(formulae = list(plate= ~ AContainers + APlates, 
                                                  field= ~ Sites + ABlocks + ALots,
                                                  fldtrts= ~ Harvesters, 
                                                  labtrts= ~ Treats),
                                  labels = "terms", data = corn.lay)
  summary(corn2Res.canon, which.criteria="none")
  testthat::expect_equal(corn2Res.canon$Q[[4]]$'AContainers&ALots&Residual&Treats', 8)
  testthat::expect_equal(corn2Res.canon$Q[[4]]$'AContainers&ALots&Residual&Residual', 146)
  
})


cat("#### Test for plaid example\n")
test_that("plaid", {
  skip_on_cran()
  library(dae)
  ph1.sys <- cbind(fac.gen(list(Expressive = c("Yes", "No"), Patients = 4, Occasions = 2)),
                   fac.gen(list(Motions = c("active", "passive")), times = 8))
  
  #'## Generate the two-phase systematic design
  ph2.sys <- cbind(fac.gen(list(Raters = 74, Viewings = 16)),
                   fac.gen(list(Trainings = 2, 16), times = 37),
                   rep(ph1.sys, times =74))
  
  #'## Randomize the two-phase design
  ph2.lay <- designRandomize(allocated = ph2.sys[c("Trainings", "Expressive", "Patients",
                                                   "Occasions", "Motions")],
                             recipient = ph2.sys[c("Raters", "Viewings")],
                             except = "Viewings",
                             seed = 15674)
  testthat::expect_true(all(dim(ph2.lay)==c(1184,7)))

  #'## Produce the anatomy of the design for the initial allocation model
  ph2.canon <- designAnatomy(formulae = list(rate  = ~ Raters * Viewings,
                                             video = ~ (Expressive/Patients)*Occasions,
                                             alloc = ~ Trainings * Motions),
                             data = ph2.lay)
  summ <- summary(ph2.canon, which.criteria = "aeff")
  testthat::expect_true(is.null(summ$aliasing))
  testthat::expect_true(all(dim(summ$decomp)==c(9,7)))
  
  #'## Convert the names of the factors to single capital letters
  ph2.L.lay <- ph2.lay
  names(ph2.L.lay)[match(c("Raters", "Viewings", "Trainings", "Expressive", "Patients", 
                           "Occasions", "Motions"), names(ph2.L.lay))] <- c("R", "V", "T", 
                                                                            "E", "P", "O", "M")
  # Produce an anatomy for the homogeneous allocation model
  # This anatomy is not the correct anatomy. 
  # However, it is an anatomy that did not display correctly 
  ph2.homog.canon <- designAnatomy(formulae = list(rate = ~ R * V, 
                                                   video = ~ E/P/O, #superfluous
                                                   inter = ~ R * O * (E/P), 
                                                   alloc = ~ T * M * (E/P)),
                                   data = ph2.L.lay)
  summ <- summary(ph2.homog.canon, which.criteria = "aeff")
  testthat::expect_true(is.null(summ$aliasing))
  testthat::expect_true(all(dim(summ$decomp)==c(17,9)))
  tab <- capture.output(print(summ))
  testthat::expect_true(sum(grepl("R#O", tab, fixed = TRUE))==3)
})


cat("#### Test for Piepho example with pseudoreplication\n")
test_that("Piepho", {
  skip_on_cran()
  library(dae)
  #'## Piepho example with pseudoreplication
  data('PiephoLSDRand')
  piepho.canon <- designAnatomy(formulae = list(lab= ~ Times/(Locations*Ovens),
                                                field= ~ Block/Plot/Sample,
                                                trt= ~ Harvest*Method), 
                                labels = "terms", data = Piepho_LSD_Rand)
  summary(piepho.canon, which.criteria="aeff")
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times&Block:Plot&Harvest", 3)                              
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Locations&Block", 2)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Locations&Block:Plot", 6)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Ovens&Block:Plot:Sample", 8)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Locations:Ovens&Block:Plot:Sample&Method",2)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Locations:Ovens&Block:Plot:Sample&Harvest:Method", 6)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times:Locations:Ovens&Block:Plot:Sample&Residual", 8)

  Piepho_LSD_Rand <- within(Piepho_LSD_Rand, 
                            {
                              ALocations <- fac.combine(list(Times, Locations))
                              AOvens <- fac.combine(list(Times, Ovens))
                              APositions <- fac.combine(list(Times, Locations, Ovens))
                              APlot <- fac.combine(list(Block, Plot))
                              ASample <- fac.combine(list(Block, Plot, Sample))
                            })
  piephoA.canon <- designAnatomy(formulae = list(lab= ~ Times + ALocations + AOvens + APositions,
                                                 field= ~ Block + APlot + ASample,
                                                 trt= ~ Harvest*Method),
                                 labels = "terms", data = Piepho_LSD_Rand)
  summary(piephoA.canon, which.criteria="aeff")
  testthat::expect_equal(piephoA.canon$Q[[3]]$"Times&APlot&Harvest", 3)                              
  testthat::expect_equal(piephoA.canon$Q[[3]]$"ALocations&Block", 2)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"ALocations&APlot", 6)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"AOvens&ASample", 8)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions&ASample&Method",2)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions&ASample&Harvest:Method", 6)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions&ASample&Residual", 8)
})


cat("#### FAME example with complicated pseudoterm structure\n")
test_that("FAME", {
  skip_on_cran()
  library(dae)
  "FAME example with complicated pseudoterm structure"
  data('FAME')
  FAME[c("Int1","Int2","Int3","Int4","Int5","Int6")] <- 
    lapply(FAME[c("Int1","Int2","Int3","Int4","Int5","Int6")], factor)
  
  fame.canon <- designAnatomy(formulae = list(lab= ~ Int1:Int2:Int3/Int4/Int5/Int6,
                                              field= ~ ((Block/Plot)*Depth)/Sample,
                                              trt= ~ Tillage*Depth*Method),
                              labels = "terms", data = FAME)
  summary(fame.canon, which.criteria="aeff")
  testthat::expect_equal(length(fame.canon$Q[[2]]), 14)
  testthat::expect_equal(length(fame.canon$Q[[3]]), 20)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3&Block:Plot:Depth:Sample&Residual', 1)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3:Int4:Int5&Block:Plot&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3:Int4:Int5&Block:Plot:Depth:Sample&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3:Int4:Int5:Int6&Block:Plot:Depth&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3:Int4:Int5:Int6&Block:Plot:Depth:Sample&Residual', 6)
  
  #Form factors with unique levels 
  FAME <- within(FAME, 
                 {
                   Int2 <- fac.combine(list(Int1,Int2))
                   Int3 <- fac.combine(list(Int2,Int3))
                   Int4 <- fac.combine(list(Int3,Int4))
                   Int5 <- fac.combine(list(Int4,Int5))
                   Int6 <- fac.combine(list(Int5,Int6))
                   Plot <- fac.combine(list(Block,Plot))
                   Sample <- fac.combine(list(Plot,Depth,Sample))
                 })
  fameA.canon <- designAnatomy(formulae = list(lab= ~ Int3 + Int4 + Int5 + Int6,
                                               field= ~ (Block + Plot)*Depth + Sample,
                                               trt= ~ Tillage*Depth*Method),
                               labels = "terms", data = FAME)
  summary(fameA.canon, which.criteria="aeff")
  testthat::expect_equal(length(fameA.canon$Q[[2]]), 14)
  testthat::expect_equal(length(fameA.canon$Q[[3]]), 20)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int3&Sample&Residual', 1)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int5&Plot&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int5&Sample&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int6&Plot:Depth&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int6&Sample&Residual', 6)
  
})


cat("#### Test for Piracicaba Euc pulp example\n")
test_that("EucPulp", {
  skip_on_cran()
  library(dae)
  "Piracicaba Euc pulp example- an example with no Residuals"
  data('EucPulp2x2')
  euc.canon <- designAnatomy(formulae = list(lab= ~ Runs/Positions,  
                                             cook= ~ Cookings/Samples,
                                             trt= ~ ((Kinds*Ages)/Lots/Batches)*Times),
                             labels = "terms", data = Euc_Pulp2x2)
  summary(euc.canon, which.criteria="aeff")
  testthat::expect_equal(length(euc.canon$Q[[3]]), 11)
  testthat::expect_equal(euc.canon$Q[[3]]$'Runs&Cookings&Kinds:Ages:Lots:Batches', 36)
  testthat::expect_equal(euc.canon$Q[[3]]$'Runs:Positions&Cookings:Samples&Kinds:Ages:Lots:Batches:Times', 180)
  
})


cat("#### Test for Split plot with rows and columns in main and split-plots\n")
test_that("SpliPlotRowsColumns", {
  skip_on_cran()
  library(dae)
  #'## Split plot with rows and columns in main and split-plots
  data('LS4.CC4x4.1a')
  split.recip <- LS4_CC4x4_1a[, 1:4]
  split.nest <- list(SubCols = "BigRows")
  split.trt <- LS4_CC4x4_1a[, 5:8]
  split.layout <- designRandomize(allocated = split.trt, 
                                  recipient = split.recip, 
                                  nested.recipients = split.nest, 
                                  seed = 7154)

  #'### Check design
  split.canon <- designAnatomy(formulae = list(plot= ~ ((BigRows/SubCols)*BigCols)*SubRows,
                                               trts= ~ Soils*Varieties*TreatB*MF),
                               labels = "terms", data = split.layout)
  summary(split.canon, which.criteria=c("aeff","ord"))
  testthat::expect_equal(length(split.canon$Q[[1]]), 11)
  testthat::expect_equal(split.canon$Q[[1]]$'BigRows:SubCols'$'Soils:TreatB:MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(split.canon$Q[[1]]$'BigRows:SubCols:BigCols'$'Soils:TreatB:MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(split.canon$Q[[1]]$'BigRows:SubCols:SubRows'$'Soils:Varieties:TreatB:MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(split.canon$Q[[1]]$'BigRows:SubCols:BigCols:SubRows'$'Soils:Varieties:TreatB:MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(split.canon$Q[[2]]$'BigRows:SubCols:BigCols:SubRows&Residual', 72)

  split.layout <- within(split.layout, 
                         {
                           ASubCols <- fac.combine(list(BigRows, SubCols))
                         })
  splitA.canon <- designAnatomy(formulae = list(plot= ~ ((BigRows + ASubCols)*BigCols)*SubRows,
                                                trts= ~ Soils*Varieties*TreatB*MF),
                                labels = "terms", data = split.layout)
  summary(splitA.canon, which.criteria=c("aeff","ord"))
  testthat::expect_equal(length(splitA.canon$Q[[1]]), 11)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols'$'Soils:TreatB:MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols:BigCols'$'Soils:TreatB:MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols:SubRows'$'Soils:Varieties:TreatB:MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols:BigCols:SubRows'$'Soils:Varieties:TreatB:MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(splitA.canon$Q[[2]]$'ASubCols:BigCols:SubRows&Residual', 72)
  
})


cat("#### Test for EXP249 - a two-phase, p-rep design\n")
test_that("Exp249", {
  skip_on_cran()
  library(dae)
  data(file="Exp249.mplot.sys")
  Exp249.mplot.sys$Blocks <- factor(rep(1:6, each = 44))
  
  #'## Expand design to rerandomize lines and to assign conditions to locations
  Exp249.recip <- list(Groups = 6, Columns = 4, Pairs = 11, Locations = 2)
  Exp249.nest <- list(Columns = c("Groups", "Pairs"),
                      Locations = c("Groups", "Columns", "Pairs"))
  Exp249.alloc <- data.frame(Lines = factor(rep(Exp249.mplot.sys$Lines, each=2), 
                                           levels=1:75),
                            Checks = fac.recode(rep(Exp249.mplot.sys$Lines, each=2), 
                                                newlevels=c(rep(3, 73), 1 , 2), 
                                                labels = c("NAM","Scout","Gladius")),
                            Conditions = factor(rep(1:2, times=264), 
                                                labels = c('0 NaCl','100 NaCl')))
  Exp249.lay <- designRandomize(allocated = Exp249.alloc, 
                                recipient = Exp249.recip, 
                                nested.recipients = Exp249.nest, 
                                seed = 51412)
  
  
  
  #'## Add second-phase factors 
  #'## (to which the first-phase factors have been systematically allocated)
  Exp249.lay <- cbind(fac.gen(list(Lanes = 24, Positions = 2:23)),
                      fac.gen(list(Zones = 6, Rows = 4, MainPosn = 11, Subplots = 2)), 
                      Exp249.lay)
  
  #'## Check design properties
  Exp249.canon <- designAnatomy(formulae = list(carts = ~(Zones*MainPosn)/Rows/Subplots, 
                                                tables = ~(Groups*Pairs)/Columns/Locations,
                                                treats = ~(Checks + Lines) * Conditions),
                                labels = "terms", data = Exp249.lay)
  summary(Exp249.canon)
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 5)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Zones&Groups'$'Lines'$adjusted$aefficiency - 0.1498311), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Zones:MainPosn:Rows&Groups:Pairs:Columns'$'Lines'$adjusted$aefficiency - 0.6639769), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Zones:MainPosn:Rows:Subplots&Groups:Pairs:Columns:Locations'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249.canon$Q[[3]]$'Zones:MainPosn:Rows&Groups:Pairs:Columns&Residual', 124)
  testthat::expect_equal(Exp249.canon$Q[[3]]$'Zones:MainPosn:Rows:Subplots&Groups:Pairs:Columns:Locations&Residual', 189)
  
  #'## Add factors and variates for new analysis
  Exp249.lay <- within(Exp249.lay, 
                       { xMainPosn <- as.numfac(MainPosn)
                       xMainPosn <- -(xMainPosn - mean(xMainPosn))
                       Mainplots <- fac.combine(list(Rows,MainPosn))
                       })

  #'## Check properties if only linear trend fitted
  Exp249.canon <- designAnatomy(formulae = list(cart = ~Zones/Mainplots/Subplots, 
                                                treat = ~xMainPosn + (Checks + Lines) * Conditions),
                                labels = "terms", data = Exp249.lay, 
                                orthogonalize = c("diff", "eigenmethods"))
  summary(Exp249.canon, which.criteria = c("aeff", "xeff", "eeff", "order"))
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 3)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Zones'$'Lines'$adjusted$aefficiency - 0.1499827), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Zones:Mainplots'$'Lines'$adjusted$aefficiency - 0.987877), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Zones:Mainplots:Subplots'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249.canon$Q[[2]]$'Zones:Mainplots&Residual', 183)
  testthat::expect_equal(Exp249.canon$Q[[2]]$'Zones:Mainplots:Subplots&Residual', 189)
  
  Exp249.lay <- within(Exp249.lay, 
                       {
                         AMainplots <- fac.combine(list(Zones, Mainplots))
                         ASubplots <- fac.combine(list(AMainplots,Subplots))
                       })
  set.daeTolerance(1e-06)
  Exp249A.canon <- designAnatomy(formulae = list(cart = ~Zones + AMainplots + ASubplots, 
                                                 treat = ~xMainPosn + (Checks + Lines) * Conditions),
                                 labels = "terms", data = Exp249.lay)
  summary(Exp249A.canon)
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 3)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'Zones'$'Lines'$adjusted$aefficiency - 0.1499827), 1e-05)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'AMainplots'$'Lines'$adjusted$aefficiency - 0.987877), 1e-05)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'ASubplots'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249A.canon$Q[[2]]$'AMainplots&Residual', 183)
  testthat::expect_equal(Exp249A.canon$Q[[2]]$'ASubplots&Residual', 189)
  
})

