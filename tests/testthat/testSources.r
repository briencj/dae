#devtools::test("dae")
context("canonical")

cat("#### Test for source formation\n")
test_that("sources", {
  skip_on_cran()
  library(dae)
  #((A*B)/C)*D
  ABCD.lay <- fac.gen(list(A = 3, B = 3, C = 3, D = 3))
  ABCD.struct <- pstructure(~ ((A*B)/C)*D, labels = "sources", data =ABCD.lay)
  ABCD.dat <- as.data.frame(ABCD.struct)
  testthat::expect_equal(names(ABCD.struct$Q), ABCD.dat$sources)
  testthat::expect_equal(ABCD.dat$sources,
                         c("A","B","A#B","C[A:B]","D","A#D","B#D","A#B#D","C#D[A:B]"))
  ABCD.struct <- pstructure(~ ((A*B)/C)*D, omit.projectors = TRUE, labels = "terms", data =ABCD.lay)
  ABCD.dat <- as.data.frame(ABCD.struct)
  testthat::expect_equal(ABCD.dat$df, c(2,2,4,18,2,4,4,8,36))
  
  
  #Generalized factor crossed with others
  genfac1.lay <- fac.gen(list(S = 3, A = 3, D = 3))
  gf1.struct <- pstructure(~ (S:A)*D, labels = "sources", data =genfac1.lay)
  testthat::expect_equal(unname(gf1.struct$sources), c("S:A","D","(S:A)#D"))
  
  #Generalized factor nested within another factor
  gf2.struct <- pstructure(~ D/(S:A), labels = "sources", data =genfac1.lay)
  testthat::expect_equal(unname(gf2.struct$sources), c("D","S:A[D]"))
  
  #Generalized factor nesting another factor
  gf3.struct <- pstructure(~ (S:A)/D, labels = "sources", data =genfac1.lay)
  testthat::expect_equal(unname(gf3.struct$sources), c("S:A","D[S:A]"))
  
  #Two generalized factor nested within another and each other
  genfac4.lay <- fac.gen(list(S = 3, A = 3, D = 3, B = 3, C = 3))
  gf4.struct <- pstructure(~ D/(S:A)/(B:C), labels = "sources", data =genfac4.lay)
  testthat::expect_equal(unname(gf4.struct$sources), c("D","S:A[D]","B:C[D:S:A]"))

  #Two generalized factor nested within another and each other - add grand mean
  gf4gm.struct <- pstructure(~ D/(S:A)/(B:C), labels = "sources", grandMean = TRUE, data =genfac4.lay)
  testthat::expect_equal(unname(gf4gm.struct$sources), c("Mean", "D","S:A[D]","B:C[D:S:A]"))
  testthat::expect_equal(nrow(as.data.frame(gf4gm.struct, omit.marginality = TRUE)), 4)

  #Single term
  one.struct <- pstructure(~ D, labels = "sources", data =genfac4.lay)

  #Single term with grand mean
  one.struct <- pstructure(~ D, labels = "sources", grandMean = TRUE, data =genfac4.lay)
  testthat::expect_equal(unname(one.struct$terms), c("Mean", "D"))
  testthat::expect_equal(unname(one.struct$sources), c("Mean", "D"))
  testthat::expect_equal(nrow(as.data.frame(one.struct, omit.marginality = TRUE)), 2)
})


cat("#### Test for designAnatomy for Thao designs\n")
test_that("Thao", {
  skip_on_cran()
  library(dae)

  Br <- 4
  Bc <- 4
  Sr<-4
  Sc<-4
  n <- Sr*Sc*Br*Bc
  v<- 4
  s<-4
  
  ##### Poset 2f design 2
  # Factor D and E are latinised using LCCD.
  LS1 <- factor(c(1:4, 4:1, 2,1,4,3,3,4,1,2 ))
  LS2 <- factor(c(1:4,2,1,4,3,3,4,1,2,4,3,2,1  ))
  LS3 <- factor(c(1:4, 3,4,1,2, 4:1, 2,1,4,3))
  LS4 <- factor(c(1:4, 2,1,4,3, 3,4,1,2, 4:1))
  LS4.CC4x4.2f.ran <- data.frame(
    A = factor(c(rep(c(1,2,1,2 ,2,1,2,1 ,1,2,1,2 ,2,1,2,1), each=16))), 
    B =  factor(c(rep(c(1,1,2,2 ,1,1,2,2 ,2,2,1,1 ,2,2,1,1), each=16))), 
    D = factor(c(rep(LS1[1:4],each=4), rep(LS1[5:8],each=4),rep(LS1[9:12],each=4),rep(LS1[13:16],each=4)
                 ,rep(LS2[1:4],each=4), rep(LS2[5:8],each=4),rep(LS2[9:12],each=4),rep(LS2[13:16],each=4)
                 ,rep(LS3[1:4],each=4), rep(LS3[5:8],each=4),rep(LS3[9:12],each=4),rep(LS3[13:16],each=4)
                 ,rep(LS4[1:4],each=4), rep(LS4[5:8],each=4),rep(LS4[9:12],each=4),rep(LS4[13:16],each=4))),
    
    E = factor(c(rep(LS1[1:4],times=4),rep(LS2[1:4],times=4),rep(LS3[1:4],times=4),rep(LS4[1:4],times=4)
                  ,rep(LS1[5:8],times=4), rep(LS2[5:8],times=4), rep(LS3[5:8],times=4),rep(LS4[5:8],times=4)
                  ,rep(LS1[9:12],times=4),rep(LS2[9:12],times=4),rep(LS3[9:12],times=4),rep(LS4[9:12],times=4)
                  ,rep(LS1[13:16],times=4),rep(LS2[13:16],times=4),rep(LS3[13:16],times=4) ,rep(LS4[13:16],times=4))))
  
  
  #generate layout
  LS4.CC4x4.2f.lay <- designRandomize(allocated=LS4.CC4x4.2f.ran,
                                      recipient = list(BigRows= Br, BigColumns=Bc, Rows=Sr, Columns=Sc),
                                      nested.recipients = list(Columns="BigColumns", Rows="BigRows"),
                                      seed = 550)
  
  #Compute anatomy
  LS4.CC4x4.2f.canon <- designAnatomy(list(unit = ~ (BigRows/Rows)*(BigColumns/Columns),
                                           trt = ~ A*B*D*E),
                                      data = LS4.CC4x4.2f.lay)
  testthat::expect_true(all(LS4.CC4x4.2f.canon$sources$unit == 
                              c("BigRows", "Rows[BigRows]", "BigColumns", "Columns[BigColumns]", 
                                "BigRows#BigColumns", "BigRows#Columns[BigColumns]",
                                "Rows#BigColumns[BigRows]", "Rows#Columns[BigRows:BigColumns]")))
  summary(LS4.CC4x4.2f.canon, which.criteria = c("aeff", "xeff", "eeff", "order"))

  #### Poset 3c, design 1 - an N-poset
  #Set up allocated factors
  LS4.CC4x4.3c.ran <- cbind(data.frame(A = factor(rep(c(1,1,2,2 ,1,1,2,2 ,2,2,1,1 ,2,2,1,1), each = 16)),
                                       B = factor(rep(c(1,2,1,2 ,2,1,2,1 ,1,2,1,2 ,2,1,2,1), each = 16))), 
                           fac.gen(list(D = 4, E = 4), times = 16))
  #generate layout and analyze
  LS4.CC4x4.3c.lay <- designRandomize(allocated=LS4.CC4x4.3c.ran, 
                                      recipient =  list(BigRows= Br, BigColumns=Bc, Rows=Sr, Columns=Sc),
                                      nested.recipients = list(Rows="BigRows", Columns=c("BigRows","BigColumns")),
                                      seed = 550)

  #Compute anatomy
  LS4.CC4x4.3c.canon <- designAnatomy(list(unit = ~ (BigRows/Rows)*BigColumns +((BigColumns*BigRows)/Columns)/Rows,
                                           trt = ~ A*B*D*E), 
                                      data =LS4.CC4x4.3c.lay)
  testthat::expect_true(all(LS4.CC4x4.3c.canon$sources$unit == 
                              c("BigRows", "Rows[BigRows]", "BigColumns", "BigRows#BigColumns", 
                                "Rows#BigColumns[BigRows]", "Columns[BigRows:BigColumns]", 
                                "Rows#Columns[BigRows:BigColumns]")))
  summary(LS4.CC4x4.3c.canon, which.criteria = c("aeff", "order"))
  
  ### Poset 4a
  LS4x4.CC5x5.4a.ran <- cbind(data.frame(A = factor(rep(c(1,2,1,2 ,2,1,2,1 ,1,2,1,2 ,2,1,2,1), each = 16)),
                                         B = factor(rep(c(1,1,2,2 ,1,1,2,2 ,2,2,1,1 ,2,2,1,1), each = 16))), 
                              fac.gen(list(D = 4, E = 4), times = 16))


  #generate layout and analyze
  LS4x4.CC5x5.4a.lay <- designRandomize(allocated=LS4x4.CC5x5.4a.ran,
                                        recipient = list(BigRows= Br, BigColumns=Bc, Rows=Sr, Columns=Sc), 
                                        nested.recipients = list(Columns=c("BigRows","BigColumns"), 
                                                              Rows=c("BigRows","BigColumns")),
                                        seed = 690)

  #Compute anatomy
  LS4x4.CC5x5.4a.lcanon <- designAnatomy(list(unit = ~ (BigRows*BigColumns)/(Rows*Columns),
                                              trt = ~ A*B*D*E), 
                                         data =LS4x4.CC5x5.4a.lay)
  testthat::expect_true(all(LS4x4.CC5x5.4a.lcanon$sources$unit == 
                              c("BigRows", "BigColumns", "BigRows#BigColumns", "Rows[BigRows:BigColumns]", 
                                "Columns[BigRows:BigColumns]", "Rows#Columns[BigRows:BigColumns]")))
  summary(LS4x4.CC5x5.4a.lcanon, which.criteria = c("aeff", "order"))
})

cat("#### Test for designAnatomy with sources and marginality using Cochran&Cox PBIBD2\n")
test_that("PBIBD2_sources", {
  skip_on_cran()
  library(dae)
  #'# PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 2nd edn Wiley, New York"
  
  #'## Input the design and randomize"
  Treatments <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
  PBIBD2.lay <- designRandomize(allocated = Treatments,
                                recipient = list(Blocks =6, Units = 4),
                                nested.recipients = list(Units = "Blocks"),
                                seed = 98177)
  
  #'## Test that projs.canon has been deprecated
  testthat::expect_warning(PBIBD2.canon <- 
                             projs.canon(formulae = list(unit = ~Blocks/Units,
                                                         trt = ~ Treatments),
                                         which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                         labels = "sources", data = PBIBD2.lay))

  #'##By differencing
  PBIBD2D.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                trt = ~ Treatments),
                which.criteria = c('aeff', 'xeff', 'eeff','order'),
                labels = "sources", orthogonalize = "diff", data = PBIBD2.lay)
  summary(PBIBD2D.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(PBIBD2D.canon$Q[[1]]$Blocks$Treatments$adjusted$aefficiency, 0.25)
  testthat::expect_lt(abs(PBIBD2D.canon$Q[[1]]$'Units[Blocks]'$Treatments$adjusted$aefficiency - 0.8824), 1e-04)
  testthat::expect_equal(PBIBD2D.canon$Q[[2]]$'Blocks&Treatments', 2)
  testthat::expect_equal(PBIBD2D.canon$Q[[2]]$'Blocks&Residual', 3)
  testthat::expect_equal(PBIBD2D.canon$Q[[2]]$'Units[Blocks]&Treatments', 5)
  testthat::expect_equal(PBIBD2D.canon$Q[[2]]$'Units[Blocks]&Residual', 13)
  testthat::expect_lt(abs(PBIBD2D.canon$Q[[1]]$'Units[Blocks]'$Treatments$adjusted$aefficiency - 
                            0.8824), 1e-04)
  testthat::expect_null(PBIBD2D.canon$aliasing$unit)
  testthat::expect_null(PBIBD2D.canon$aliasing$trt)
  PBIBD2.marg <- PBIBD2D.canon$marginality
  
  #'##By differencing with grand mean
  PBIBD2Dg.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  grandMean = TRUE, labels = "sources", orthogonalize = "diff", 
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  data = PBIBD2.lay)
  summary(PBIBD2Dg.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(PBIBD2Dg.canon$Q[[2]]$`Mean&Mean`, 1)
  
  #'##By eigenmethods
  PBIBD2E.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  labels = "sources", orthogonalize = "eigen", data = PBIBD2.lay)
  summary(PBIBD2E.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_true(attr(PBIBD2E.canon, which = "labels") == "terms")
  testthat::expect_lt(abs(PBIBD2E.canon$Q[[1]]$'Blocks:Units'$Treatments$adjusted$aefficiency - 
                            0.8824), 1e-04)

  #'##By eigenmethods with supplied marginality
  PBIBD2Em.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  marginality = PBIBD2D.canon$marginality, 
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  labels = "sources", orthogonalize = "eigen", data = PBIBD2.lay)
  summary(PBIBD2Em.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(length(PBIBD2Em.canon$marginality), 2)
  testthat::expect_true(attr(PBIBD2Em.canon, which = "labels") == "sources")
  
  #marginality list with trt NULL - same as supplying just unit marginality
  PBIBD2Em.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  marginality = PBIBD2.marg, 
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  labels = "sources", orthogonalize = "eigen", data = PBIBD2.lay)
  summary(PBIBD2Em.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(length(PBIBD2Em.canon$marginality), 2)
  testthat::expect_true(attr(PBIBD2Em.canon, which = "labels") == "sources")
  
  #full marginality list
  PBIBD2Em.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  marginality = PBIBD2.marg, 
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  labels = "sources", orthogonalize = "eigen", data = PBIBD2.lay)
  summary(PBIBD2Em.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(length(PBIBD2Em.canon$marginality), 2)
  testthat::expect_true(attr(PBIBD2Em.canon, which = "labels") == "sources")
  testthat::expect_true(all(names(PBIBD2Em.canon$marginality) == names(PBIBD2.marg)))
  testthat::expect_lt(abs(PBIBD2Em.canon$Q[[1]]$'Units[Blocks]'$Treatments$adjusted$aefficiency - 
                            0.8824), 1e-04)
  #Add grand mean term
  PBIBD2Egm.canon <- 
    designAnatomy(formulae = list(unit = ~Blocks/Units,
                                  trt = ~ Treatments),
                  grandMean = TRUE, marginality = PBIBD2.marg, 
                  which.criteria = c('aeff', 'xeff', 'eeff','order'),
                  labels = "sources", orthogonalize = "eigen", data = PBIBD2.lay)
  summary(PBIBD2Egm.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_equal(length(PBIBD2Em.canon$marginality), 2)
  testthat::expect_true(attr(PBIBD2Em.canon, which = "labels") == "sources")
  testthat::expect_true(all(PBIBD2Egm.canon$terms$unit == c("Mean", "Blocks", "Blocks:Units")))
  testthat::expect_true(all(PBIBD2Egm.canon$terms$trt == c("Mean", "Treatments")))
  testthat::expect_true(all(names(PBIBD2Em.canon$marginality) == names(PBIBD2.marg)))
  testthat::expect_lt(abs(PBIBD2Em.canon$Q[[1]]$'Units[Blocks]'$Treatments$adjusted$aefficiency - 
                            0.8824), 1e-04)
  #'## unique plot levels
  PBIBD2.lay$AUnits <- with(PBIBD2.lay, fac.combine(list(Blocks,Units)))
  PBIBD2A.canon <- designAnatomy(list(unit = ~Blocks + AUnits,
                                      trt = ~ Treatments),
                                 data = PBIBD2.lay, labels = "sources", 
                                 which.criteria = c('aeff', 'xeff', 'eeff','order'),
                                 orthogonalize = "hybrid")
  summary(PBIBD2A.canon, which.criteria = c('aeff', 'xeff', 'eeff','order'))
  testthat::expect_lt(abs(PBIBD2A.canon$Q[[1]]$'AUnits[Blocks]'$Treatments$adjusted$aefficiency - 
                            0.8824), 1e-04)
  
})

cat("#### Test for Jarrett & Ruggiero example\n")
test_that("JarrettRuggiero", {
  skip_on_cran()
  library(dae)
  #'## Jarrett & Ruggiero example
  jr.lay <- fac.gen(list(Set=7, Dye=2, Array=3))
  jr.lay <- within(jr.lay, 
                   { Block <- factor(rep(1:7, each=6))
                   Plant <- factor(rep(c(1,2,3,2,3,1), times=7))
                   Sample <- factor(c(rep(c(2,1,2,2,1,1, 1,2,1,1,2,2), times=3), 2,1,2,2,1,1))
                   S1 <- Dye
                   Treat <- factor(c(1,2,4,2,4,1, 2,3,5,3,5,2, 3,4,6,4,6,3, 4,5,7,5,7,4, 
                                     5,6,1,6,1,5, 6,7,2,7,2,6, 7,1,3,1,3,7),
                                   labels=c("A","B","C","D","E","F","G"))
                   })
  

  array.plot.trt.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye,
                                                        plot = ~ Block/Plant/Sample,
                                                        trt = ~ Treat),
                                        labels = "sources", data = jr.lay)
  
  testthat::expect_equal(length(array.plot.trt.canon$Q[[3]]), 7)
  testthat::expect_equal(array.plot.trt.canon$Q[[3]]$`(Set:Array)#Dye&Sample[Block:Plant]`, 6)
  testthat::expect_lt(abs(array.plot.trt.canon$Q[[2]]$`(Set:Array)#Dye&Plant[Block]`$Treat$adjusted$aefficiency - 0.58333333), 1e-05)
  testthat::expect_lt(abs(array.plot.trt.canon$Q[[1]]$`(Set:Array)#Dye`$`Plant[Block]`$adjusted$aefficiency - 0.75), 1e-05)
  
  summ.default <-summary(array.plot.trt.canon)
  testthat::expect_equal(nrow(summ.default$decomp), 7)
  testthat::expect_equal(ncol(summ.default$decomp), 9)
  testthat::expect_false(is.null(attr(summ.default$aliasing, which = "title")))
  summ.none <-summary(array.plot.trt.canon, which.criteria = "none")
  testthat::expect_equal(ncol(summ.none$decomp), 6)
  summ.all <-summary(array.plot.trt.canon, which.criteria = "all")
  testthat::expect_equal(ncol(summ.all$decomp), 13)
  summ.2 <-summary(array.plot.trt.canon, which.criteria =c("aefficiency", "order"))
  testthat::expect_equal(ncol(summ.2$decomp), 8)
 
  #with grand mean 
  apt.gm.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye,
                                                plot = ~ Block/Plant/Sample,
                                                trt = ~ Treat),
                                grandMean = TRUE, labels = "sources", data = jr.lay)
  summary(apt.gm.canon)
  testthat::expect_equal(apt.gm.canon$Q[[3]]$`Mean&Mean&Mean`, 1)
  
  #test eigen with supplied marginality for tier 2 only
  plot.marg <- array.plot.trt.canon$marginality$plot
  apt.tier2.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye,
                                                   plot = ~ Block/Plant/Sample,
                                                   trt = ~ Treat),
                                   marginality = list(plot = plot.marg), orthog = "eigen", 
                                   labels = "sources", data = jr.lay)
  summary(apt.tier2.canon)
  testthat::expect_true(all(apt.tier2.canon$sources$array == 
                              c("Set:Array", "Dye", "Set:Array:Dye")))
  testthat::expect_true(all(apt.tier2.canon$sources$plot == 
                              c("Block", "Plant[Block]", "Sample[Block:Plant]")))
  testthat::expect_equal(length(apt.tier2.canon$marginality), 3)
  testthat::expect_equal(nrow(apt.tier2.canon$marginality$plot), 3)

  apt.tier2.gm.canon <- designAnatomy(formulae = list(array = ~ (Set:Array)*Dye,
                                                      plot = ~ Block/Plant/Sample,
                                                      trt = ~ Treat),
                                      marginality = list(plot = plot.marg), orthog = "eigen", 
                                      grandMean = TRUE, labels = "sources", data = jr.lay)
  summary(apt.tier2.gm.canon)
  testthat::expect_true(all(apt.tier2.gm.canon$sources$array == 
                              c("Mean", "Set:Array", "Dye", "Set:Array:Dye")))
  testthat::expect_true(all(apt.tier2.gm.canon$sources$plot == 
                              c("Mean", "Block", "Plant[Block]", "Sample[Block:Plant]")))
  testthat::expect_equal(length(apt.tier2.gm.canon$marginality), 3)
  testthat::expect_equal(nrow(apt.tier2.gm.canon$marginality$plot), 3)
  testthat::expect_equal(apt.tier2.gm.canon$Q[[3]]$`Mean&Mean&Mean`, 1)
  
})


cat("#### Test for Baby pseudoterm example\n")
test_that("Baby", {
  skip_on_cran()
  library(dae)
  #'## Baby pseudoterm example
  pseudo.lay <- data.frame(pl = factor(1:12),
                           ab = factor(rep(1:4, times=3)),
                           a = factor(rep(1:2, times=6)))

  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~a+ab), 
                                labels = "sources", data = pseudo.lay)
  summary(pseudo.canon)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab[a]', 2)

  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~ab+a), 
                                labels = "sources", data = pseudo.lay)
  summ.hybrid <- summary(pseudo.canon)
  testthat::expect_equal(nrow(summ.hybrid$aliasing), 2)
  testthat::expect_equal(ncol(summ.hybrid$aliasing), 7)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab', 3)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&Residual', 8)
  
  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~ab+a), 
                                orthogonalize = "diff",
                                labels = "sources", data = pseudo.lay)
  summ.diff <- summary(pseudo.canon)
  testthat::expect_equal(nrow(summ.diff$aliasing), 2)
  testthat::expect_equal(ncol(summ.diff$aliasing), 7)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&ab', 3)
  testthat::expect_equal(pseudo.canon$Q[[2]]$'pl&Residual', 8)

  pseudo.canon <- designAnatomy(formulae = list(unit=~pl, trt=~ab+a), 
                                orthogonalize = "eigen",
                                labels = "sources", data = pseudo.lay)
  summ.eigen <- summary(pseudo.canon)
  testthat::expect_equal(nrow(summ.eigen$aliasing), 1)
  testthat::expect_equal(ncol(summ.eigen$aliasing), 7)
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
                                  data = sameRCBD.lay, 
                                  labels = "sources", grandMean = TRUE)
  testthat::expect_true(all(unlist(sameRCBD.canon$Q[[2]]) == c(1,2,2,23,46,4,46,92)))
  testthat::expect_equal(names(sameRCBD.canon$Q[[1]][4]), "Plots[Blocks]")
  #Test for use of label.swap
  summ.s <- summary(sameRCBD.canon, which = c('aeff','order'))
  testthat::expect_equal(summ.s$decomp$Source.unit[4], "Plots[Blocks]")
  summ.t <- summary(sameRCBD.canon, which = c('aeff','order'), labels.swap = TRUE)
  testthat::expect_equal(summ.t$decomp$Term.unit[4], "Blocks:Plots")
  testthat::expect_equal(names(sameRCBD.canon$Q[[1]][4]), "Plots[Blocks]")

  #'### unique numbering of levels
  sameRCBD.canon <- designAnatomy(formulae = list(unit = ~  Areas*(Blocks + Plot),
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "sources", data = sameRCBD.lay) 
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
                                  data = diffRCBD.lay, 
                                  labels = "sources", grandMean = TRUE)
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(1,2,6,23,46,138)))
  #'### unique numbering of levels
  diffRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas + Block + Plot,
                                                  trt = ~ Nitrogen*Varieties),
                                  labels = "sources", data = diffRCBD.lay)
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(2,6,23,46,138)))
  #'### Use eigen and diff instead of hybrid
  diffRCBD.canon <- designAnatomy(formulae = list(unit = ~ Areas + Block + Plot,
                                                  trt = ~ Nitrogen*Varieties),
                                  data = diffRCBD.lay, 
                                  labels = "sources", grandMean = TRUE, 
                                  orthogonalize = c("eigen", "diff"))
  summary(diffRCBD.canon, which = c('aeff','order'))
  testthat::expect_true(all(unlist(diffRCBD.canon$Q[[2]]) == c(1,2,6,23,46,138)))
  
})


cat("#### Test for Repeated LSD for Housewives example\n")
test_that("Housewives", {
  skip_on_cran()
  library(dae)
  data("LSRepeatHwife.dat")
  
  Hwife.struct <- pstructure( ~ Month:Week + Month:Hwife + 
                                Month:Week:Hwife,
                              labels = "sources", data = LSRepeatHwife.dat)
  testthat::expect_equal(nrow(Hwife.struct$aliasing), 2)
  testthat::expect_equal(ncol(Hwife.struct$aliasing), 10)
  
  Hwife.canon <- designAnatomy(formulae = list(units = ~ Month:Week + Month:Hwife + 
                                                          Month:Week:Hwife),
                               labels = "sources", data = LSRepeatHwife.dat)
  summ <- summary(Hwife.canon)
  testthat::expect_equal(nrow(summ$aliasing), 2)
  testthat::expect_equal(ncol(summ$aliasing), 7)
  testthat::expect_equal(nrow(Hwife.canon$aliasing), 2)
  testthat::expect_equal(ncol(Hwife.canon$aliasing), 11)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Month:Week`, 7)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Month:Hwife`, 6)
  testthat::expect_equal(Hwife.canon$Q[[1]]$`Week#Hwife[Month]`, 18)

  testthat::expect_error(Hwife.diff.canon <- designAnatomy(formulae = list(units = ~ Month:Week + Month:Hwife + 
                                                      Month:Week:Hwife),
                                    labels = "sources", 
                                    orthogonalize = "diff", data = LSRepeatHwife.dat))
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

  preece1.canon <- designAnatomy(formulae = list(plot= ~ block/plot, trt= ~ T1+T2), 
                                 labels = "sources", data = preece1.lay)
  summary(preece1.canon, which.criteria = c("aeff", "order"))
  testthat::expect_equal(length(preece1.canon$Q[[1]]), 2)
  testthat::expect_equal(length(preece1.canon$Q[[1]]$block), 3)
  testthat::expect_equal(length(preece1.canon$Q[[1]]$`plot[block]`), 3)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$block$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$block$T2$adjusted$aefficiency - 0.0952), 1e-04)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$`plot[block]`$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece1.canon$Q[[1]]$`plot[block]`$T2$adjusted$aefficiency - 0.7619), 1e-04)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&T1', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&T2', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'block&Residual', 1)
  testthat::expect_equal(preece1.canon$Q[[2]]$'plot[block]&T1', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'plot[block]&T2', 4)
  testthat::expect_equal(preece1.canon$Q[[2]]$'plot[block]&Residual', 12)
  
  preece2.canon <- designAnatomy(formulae = list(plot= ~ block/plot, trt= ~ T1+T3),
                                 labels = "sources", data = preece1.lay)
  summary(preece2.canon)
  testthat::expect_equal(length(preece2.canon$Q[[1]]), 2)
  testthat::expect_equal(length(preece2.canon$Q[[1]]$block), 2)
  testthat::expect_equal(length(preece2.canon$Q[[1]]$`plot[block]`), 3)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$block$T1$adjusted$aefficiency - 0.1666667), 1e-05)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$`plot[block]`$T1$adjusted$aefficiency - 0.8333333), 1e-05)
  testthat::expect_lt(abs(preece2.canon$Q[[1]]$`plot[block]`$T3$adjusted$aefficiency - 0.8571), 1e-04)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block&T1', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'block&Residual', 5)
  testthat::expect_equal(preece2.canon$Q[[2]]$'plot[block]&T1', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'plot[block]&T3', 4)
  testthat::expect_equal(preece2.canon$Q[[2]]$'plot[block]&Residual', 12)
  
})

cat("#### Test for Mostafa's green wall experiment in 2014\n")
test_that("Mostafa", {
  skip_on_cran()
  library(dae)
  #Mostafa's green wall experiment in 2014
  data(gwall.lay)
  options(width = 100, nwarnings = 150)
  set.daeTolerance(1e-06, 1e-06)
  pot.treat.canon <- designAnatomy(formulae = list(pot = ~ Rows*Cols,
                                                   trt = ~ Species*Irrigation*Media + 
                                                     First/(SpeCarry*IrrCarry*MedCarry)), 
                                   labels = "sources", data = gwall.lay, 
                                   keep.order=TRUE)
  summary(pot.treat.canon)
  testthat::expect_equal(length(pot.treat.canon$Q[[2]]), 17)
  testthat::expect_equal(pot.treat.canon$Q[[2]]$`Rows#Cols&Residual`, 69)
  testthat::expect_lt(abs(pot.treat.canon$Q[[1]]$`Rows#Cols`$`Species#Irrigation#Media`$adjusted$aefficiency - 0.8240828), 1e-05)
  testthat::expect_lt(abs(pot.treat.canon$Q[[1]]$`Rows#Cols`$`SpeCarry#IrrCarry[First]`$adjusted$aefficiency - 0.8320062), 1e-05)
  
  trt.struct <- pstructure(formula = ~ Species*Irrigation*Media + 
                             First/(SpeCarry*IrrCarry*MedCarry), 
                           labels = "sources", data = gwall.lay)
  testthat::expect_true(all(names(trt.struct$Q) == c("Species","Irrigation","Species#Irrigation",
                                                     "Media","Species#Media","Irrigation#Media",
                                                     "Species#Irrigation#Media","First",
                                                     "SpeCarry[First]","IrrCarry[First]",
                                                     "SpeCarry#IrrCarry[First]")))
  
  #use designAnatomy to call pstructure
  trt.canon <- designAnatomy(list(trt = ~ Species*Irrigation*Media + 
                                    First/(SpeCarry*IrrCarry*MedCarry)), 
                             labels = "sources", data = gwall.lay)
  testthat::expect_true(all(names(trt.canon$Q) == c("Species","Irrigation","Species#Irrigation",
                                                    "Media","Species#Media","Irrigation#Media",
                                                    "Species#Irrigation#Media","First",
                                                    "SpeCarry[First]","IrrCarry[First]",
                                                    "SpeCarry#IrrCarry[First]")))
  
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
                              labels = "sources", data = corn.lay)
  summary(corn.canon, which.criteria="aeff")
  testthat::expect_equal(corn.canon$Q[[3]]$`Intervals&APlots[Sites:ABlocks]&Residual`, 6)
  testthat::expect_equal(corn.canon$Q[[3]]$`AContainers[Intervals]&ALots[Sites:ABlocks:APlots]&Residual`, 
                         72)
  testthat::expect_equal(corn.canon$Q[[3]]$`APlates[Intervals:AContainers]&ALots[Sites:ABlocks:APlots]`, 
                         486)
  
  #Create an example with double Residual
  corn2Res.canon <- designAnatomy(formulae = list(plate= ~ AContainers + APlates, 
                                                  field= ~ Sites + ABlocks + ALots,
                                                  fldtrts= ~ Harvesters, 
                                                  labtrts= ~ Treats),
                                  labels = "sources", data = corn.lay)
  summary(corn2Res.canon, which.criteria="none")
  testthat::expect_equal(corn2Res.canon$Q[[4]]$`AContainers&ALots[Sites:ABlocks]&Residual&Treats`, 
                         8)
  testthat::expect_equal(corn2Res.canon$Q[[4]]$`AContainers&ALots[Sites:ABlocks]&Residual&Residual`, 
                         146)
  
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
                                labels = "sources", data = Piepho_LSD_Rand)
  summary(piepho.canon, which.criteria="aeff")
  testthat::expect_equal(piepho.canon$Q[[3]]$"Times&Plot[Block]&Harvest", 3)                              
  testthat::expect_equal(piepho.canon$Q[[3]]$"Locations[Times]&Block", 2)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Locations[Times]&Plot[Block]", 6)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Ovens[Times]&Sample[Block:Plot]", 8)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Locations#Ovens[Times]&Sample[Block:Plot]&Method",2)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Locations#Ovens[Times]&Sample[Block:Plot]&Harvest#Method", 6)
  testthat::expect_equal(piepho.canon$Q[[3]]$"Locations#Ovens[Times]&Sample[Block:Plot]&Residual", 8)

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
                                 labels = "sources", data = Piepho_LSD_Rand)
  summary(piephoA.canon, which.criteria="aeff")
  testthat::expect_equal(piephoA.canon$Q[[3]]$"Times&APlot[Block]&Harvest", 3)                              
  testthat::expect_equal(piephoA.canon$Q[[3]]$"ALocations[Times]&Block", 2)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"ALocations[Times]&APlot[Block", 6)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"AOvens[Times]&ASample[Block:APlot]", 8)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions[Times]&ASample[Block:APlot]&Method",2)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions[Times]&ASample[Block:APlot]&Harvest#Method", 6)
  testthat::expect_equal(piephoA.canon$Q[[3]]$"APositions[Times]&ASample[Block:APlot]&Residual", 8)
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
                              labels = "sources", data = FAME)
  summary(fame.canon, which.criteria="aeff")
  testthat::expect_equal(length(fame.canon$Q[[2]]), 14)
  testthat::expect_equal(length(fame.canon$Q[[3]]), 20)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int1:Int2:Int3&Sample[Block:Plot:Depth]&Residual', 1)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int5[Int1:Int2:Int3:Int4]&Plot[Block]&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int5[Int1:Int2:Int3:Int4]&Sample[Block:Plot:Depth]&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int6[Int1:Int2:Int3:Int4:Int5]&Plot#Depth[Block]&Residual', 3)
  testthat::expect_equal(fame.canon$Q[[3]]$'Int6[Int1:Int2:Int3:Int4:Int5]&Sample[Block:Plot:Depth]&Residual', 6)
  

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
                               labels = "sources", data = FAME)
  summary(fameA.canon, which.criteria="aeff")
  testthat::expect_equal(length(fameA.canon$Q[[2]]), 14)
  testthat::expect_equal(length(fameA.canon$Q[[3]]), 20)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int3&Sample[Block:Plot:Depth]&Residual', 1)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int5[Int3:Int4]&Plot[Block]&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int5[Int3:Int4]&Sample[Block:Plot:Depth]&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int6[Int3:Int4:Int5]&Plot#Depth[Block]&Residual', 3)
  testthat::expect_equal(fameA.canon$Q[[3]]$'Int6[Int3:Int4:Int5]&Sample[Block:Plot:Depth]&Residual', 6)
  
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
                             labels = "sources", data = Euc_Pulp2x2)
  summary(euc.canon, which.criteria="aeff")
  testthat::expect_equal(length(euc.canon$Q[[3]]), 11)
  testthat::expect_equal(euc.canon$Q[[3]]$'Runs&Cookings&Batches[Kinds:Ages:Lots]', 36)
  testthat::expect_equal(euc.canon$Q[[3]]$'Positions[Runs]&Samples[Cookings]&Batches#Times[Kinds:Ages:Lots]', 180)
  
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
                               labels = "sources", data = split.layout)
  summary(split.canon, which.criteria=c("aeff","ord"))
  testthat::expect_equal(length(split.canon$Q[[1]]), 11)
  testthat::expect_equal(split.canon$Q[[1]]$'SubCols[BigRows]'$'Soils#TreatB#MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(split.canon$Q[[1]]$'SubCols#BigCols[BigRows]'$'Soils#TreatB#MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(split.canon$Q[[1]]$'SubCols#SubRows[BigRows]'$'Soils#Varieties#TreatB#MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(split.canon$Q[[1]]$'SubCols#BigCols#SubRows[BigRows]'$'Soils#Varieties#TreatB#MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(split.canon$Q[[2]]$'SubCols#BigCols#SubRows[BigRows]&Residual', 72)

  Q.plot <- pstructure( ~ ((BigRows/SubCols)*BigCols)*SubRows,
                             labels = "sources", data = split.layout)
  names (Q.plot)
  
  split.layout <- within(split.layout, 
                         {
                           ASubCols <- fac.combine(list(BigRows, SubCols))
                         })
  splitA.canon <- designAnatomy(formulae = list(plot= ~ ((BigRows + ASubCols)*BigCols)*SubRows,
                                                trts= ~ Soils*Varieties*TreatB*MF),
                                labels = "sources", data = split.layout)
  summary(splitA.canon, which.criteria=c("aeff","ord"))
  testthat::expect_equal(length(splitA.canon$Q[[1]]), 11)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols[BigRows]'$'Soils#TreatB#MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols#BigCols[BigRows]'$'Soils#TreatB#MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols#SubRows[BigRows]'$'Soils#Varieties#TreatB#MF'$adjusted$aefficiency, 0.25)
  testthat::expect_equal(splitA.canon$Q[[1]]$'ASubCols#BigCols#SubRows[BigRows]'$'Soils#Varieties#TreatB#MF'$adjusted$aefficiency, 0.75)
  testthat::expect_equal(splitA.canon$Q[[2]]$'ASubCols#BigCols#SubRows[BigRows]&Residual', 72)
  
})


cat("#### Test for EXP249 - a two-phae, p-rep design\n")
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
                                labels = "sources", data = Exp249.lay)
  summary(Exp249.canon)
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 5)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Zones&Groups'$'Lines'$adjusted$aefficiency - 0.1498311), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Rows[Zones:MainPosn]&Columns[Groups:Pairs]'$'Lines'$adjusted$aefficiency - 0.6639769), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[2]]$'Subplots[Zones:MainPosn:Rows]&Locations[Groups:Pairs:Columns]'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249.canon$Q[[3]]$'Rows[Zones:MainPosn]&Columns[Groups:Pairs]&Residual', 124)
  testthat::expect_equal(Exp249.canon$Q[[3]]$'Subplots[Zones:MainPosn:Rows]&Locations[Groups:Pairs:Columns]&Residual', 189)
  
  #'## Add factors and variates for new analysis
  Exp249.lay <- within(Exp249.lay, 
                       { xMainPosn <- as.numfac(MainPosn)
                       xMainPosn <- -(xMainPosn - mean(xMainPosn))
                       Mainplots <- fac.combine(list(Rows,MainPosn))
                       })

  #'## Check properties if only linear trend fitted
  Exp249.canon <- designAnatomy(formulae = list(cart = ~Zones/Mainplots/Subplots, 
                                                treat = ~xMainPosn + (Checks + Lines) * Conditions),
                                data = Exp249.lay, 
                                labels = "sources", orthogonalize = c("diff", "eigenmethods"))
  summ <- summary(Exp249.canon)
  testthat::expect_equal(nrow(summ$aliasing), 2)
  testthat::expect_equal(ncol(summ$aliasing), 7)
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 3)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Zones'$'Lines'$adjusted$aefficiency - 0.1499827), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Mainplots[Zones]'$'Lines'$adjusted$aefficiency - 0.987877), 1e-05)
  testthat::expect_lt(abs(Exp249.canon$Q[[1]]$'Subplots[Zones:Mainplots]'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249.canon$Q[[2]]$'Mainplots[Zones]&Residual', 183)
  testthat::expect_equal(Exp249.canon$Q[[2]]$'Subplots[Zones:Mainplots]&Residual', 189)
  
  Exp249.lay <- within(Exp249.lay, 
                       {
                         AMainplots <- fac.combine(list(Zones, Mainplots))
                         ASubplots <- fac.combine(list(AMainplots,Subplots))
                       })
  Exp249A.canon <- designAnatomy(formulae = list(cart = ~Zones + AMainplots + ASubplots, 
                                                 treat = ~xMainPosn + (Checks + Lines) * Conditions),
                                 labels = "sources", data = Exp249.lay)
  summary(Exp249A.canon)
  testthat::expect_equal(length(Exp249.canon$Q[[1]]), 3)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'Zones'$'Lines[Checks]'$adjusted$aefficiency - 0.1499827), 1e-05)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'AMainplots[Zones]'$'Lines[Checks]'$adjusted$aefficiency - 0.987877), 1e-05)
  testthat::expect_lt(abs(Exp249A.canon$Q[[1]]$'ASubplots[Zones:AMainplots]'$'Conditions'$adjusted$aefficiency - 1), 1e-05)
  testthat::expect_equal(Exp249A.canon$Q[[2]]$'AMainplots[Zones]&Residual', 183)
  testthat::expect_equal(Exp249A.canon$Q[[2]]$'ASubplots[Zones:AMainplots]&Residual', 189)
  
 })


cat("#### Test for Brien and Payne 3-tier sensory experiment\n")
test_that("Sensory3tier", {
  skip_on_cran()
  library(dae)
  #Three-tier sensory experiment
  data("Need3.dat")
  
  #Do short names version
  names(Need3.dat)[match(c("Occasions", "Intervals", "Sittings", "Positions", "Judges",
                           "Squares", "Rows", "Columns", "Halfplots"), names(Need3.dat))] <- 
    c("Occ", "Int", "Sit", "Pos", "Jud", "Sqr", "Row", "Col", "HPlot")
  
  #'## Complete decomposition
  Eval.Field.Treat.canon <- designAnatomy(list(eval=~ ((Occ/Int/Sit)*Jud)/Pos, 
                                               field=~ (Row*(Sqr/Col))/HPlot,
                                               treats=~ Trellis*Method), 
                                          labels = "sources", data=Need3.dat)
  summary(Eval.Field.Treat.canon, which.criteria =c("aefficiency", "order"))
  testthat::expect_equal(names(Eval.Field.Treat.canon$Q[[1]]), 
                         c("Occ","Int[Occ]","Sit[Occ:Int]","Jud","Occ#Jud","Int#Jud[Occ]",
                           "Sit#Jud[Occ:Int]","Pos[Occ:Int:Sit:Jud]"))
  testthat::expect_equal(names(Eval.Field.Treat.canon$Q[[2]]), 
                         c("Occ&Sqr","Int[Occ]","Sit[Occ:Int]&Col[Sqr]","Sit[Occ:Int]&Residual",
                           "Jud","Occ#Jud","Int#Jud[Occ]&Row","Int#Jud[Occ]&Row#Sqr",
                           "Int#Jud[Occ]&Residual","Sit#Jud[Occ:Int]&Col[Sqr]",
                           "Sit#Jud[Occ:Int]&Row#Col[Sqr]","Sit#Jud[Occ:Int]&Residual",
                           "Pos[Occ:Int:Sit:Jud]&HPlot[Row:Sqr:Col]",
                           "Pos[Occ:Int:Sit:Jud]&Residual"))
  testthat::expect_equal(names(Eval.Field.Treat.canon$Q[[3]]), 
                         c("Occ&Sqr","Int[Occ]","Sit[Occ:Int]&Col[Sqr]&Trellis",
                           "Sit[Occ:Int]&Col[Sqr]&Residual","Sit[Occ:Int]&Residual",
                           "Jud","Occ#Jud","Int#Jud[Occ]&Row","Int#Jud[Occ]&Row#Sqr",
                           "Int#Jud[Occ]&Residual","Sit#Jud[Occ:Int]&Col[Sqr]&Trellis",
                           "Sit#Jud[Occ:Int]&Col[Sqr]&Residual","Sit#Jud[Occ:Int]&Row#Col[Sqr]&Trellis",
                           "Sit#Jud[Occ:Int]&Row#Col[Sqr]&Residual","Sit#Jud[Occ:Int]&Residual",
                           "Pos[Occ:Int:Sit:Jud]&HPlot[Row:Sqr:Col]&Method",
                           "Pos[Occ:Int:Sit:Jud]&HPlot[Row:Sqr:Col]&Trellis#Method",
                           "Pos[Occ:Int:Sit:Jud]&HPlot[Row:Sqr:Col]&Residual",
                           "Pos[Occ:Int:Sit:Jud]&Residual"))
  testthat::expect_equal(Eval.Field.Treat.canon$Q[[3]]$'Sit#Jud[Occ:Int]&Col[Sqr]&Trellis', 3)
  testthat::expect_equal(Eval.Field.Treat.canon$Q[[3]]$'Sit#Jud[Occ:Int]&Col[Sqr]&Residual', 3)
  testthat::expect_equal(Eval.Field.Treat.canon$Q[[3]]$'Sit#Jud[Occ:Int]&Row#Col[Sqr]&Trellis', 3)
  testthat::expect_equal(Eval.Field.Treat.canon$Q[[3]]$'Sit#Jud[Occ:Int]&Row#Col[Sqr]&Residual', 9)
  summ <- summary(Eval.Field.Treat.canon, which.criteria =c("aefficiency", "order"))
  testthat::expect_true(all(abs(summ$decomp[11:14,7] - 
                               c(0.07407407,0.66666667,0.88888889,1.00000000)) < 1e-06))

  })
  
