#devtools::test("dae")
context("canonical")

cat("#### Test for marginality algoritm using LRCCD\n")
test_that("LRCCD", {
   skip_on_cran()
   library(dae)
   data("LRRCD.lay")

   # All factor anatomy
   LRRCD.canon <- designAnatomy(formulae = list(carts = ~ (Block*MainPosn)/BLane/Cart,
                                                trts = ~ Line*Watering),
                                grandMean = TRUE, data = LRRCD.lay)
   summary(LRRCD.canon)
   testthat::expect_true(is.null(LRRCD.canon$aliasing))
   testthat::expect_equal(LRRCD.canon$Q[[2]]$`BLane[Block:MainPosn]&Residual`, 71)
   
   # Effect of removing `Block#MainPosn`
   LRRCD.Posn.canon <- designAnatomy(formulae = list(carts = ~ (Block+MainPosn)/BLane/Cart,
                                                     trts = ~ Line*Watering),
                                     grandMean = TRUE, data = LRRCD.lay)
   summary(LRRCD.Posn.canon)
   testthat::expect_true(is.null(LRRCD.Posn.canon$aliasing))
   testthat::expect_equal(LRRCD.Posn.canon$Q[[2]]$`Block#MainPosn#BLane&Residual`, 89)
   
   # Effect of using only a linear trend for MainPosn (Block:MainPosn:BLane defines a Mainunit)
   LRRCD.xPosn.canon <- designAnatomy(formulae = list(carts = ~ Block/MainPosn:BLane/Cart,
                                                      trts = ~ xMainPosn + Line*Watering),
                                      grandMean = TRUE, omit.projectors = "combined", data = LRRCD.lay)
   summary(LRRCD.xPosn.canon, which.criteria = c("aeff", "eeff", "xeff", "ord", "dfor"))
   testthat::expect_true(!is.null(LRRCD.xPosn.canon$aliasing))
   testthat::expect_equal(nrow(LRRCD.xPosn.canon$aliasing), 2)
   testthat::expect_equal(LRRCD.xPosn.canon$aliasing$Alias, c("xMainPosn", "## Information remaining"))
   testthat::expect_true(all(abs(LRRCD.xPosn.canon$aliasing$aefficiency - c(0.2749455, 0.9923205)) < 1e-05))
   testthat::expect_equal(LRRCD.xPosn.canon$Q[[2]]$`MainPosn:BLane[Block]&Residual`, 97)
   
   # Covariate followed by factor interaction
   print(p <- pstructure(formula = ~ xMainPosn + Line*Watering, data = LRRCD.lay, 
                         aliasing.print = FALSE), 
         which = "aliasing", which.criteria = c("aeff", "eeff", "xeff", "ord"))
   testthat::expect_true(!is.null(p$aliasing))
   testthat::expect_equal(nrow(p$aliasing), 2)
   testthat::expect_equal(p$aliasing$Alias, c("xMainPosn", "## Information remaining"))
   testthat::expect_true(all(abs(p$aliasing$aefficiency - c(0.2749455, 0.9923205)) < 1e-05))

   # Covariate followed by factor interaction - eigen orthogonalize
   testthat::expect_warning(
      print(p <- pstructure(formula = ~ xMainPosn + Line*Watering, data = LRRCD.lay, 
                            aliasing.print = FALSE, orthogonalize = "eigen"), 
            which = "aliasing", which.criteria = c("aeff", "eeff", "xeff", "ord")))
   testthat::expect_true(!is.null(p$aliasing))
   testthat::expect_equal(nrow(p$aliasing), 1)
   testthat::expect_equal(p$aliasing$Alias, "unknown")
   testthat::expect_true(abs(p$aliasing$aefficiency - 0.2749455) < 1e-05)
   
   # Factor interaction followed by a Covariate
   print(p <- pstructure(formula = ~ Line*Watering + xMainPosn, data = LRRCD.lay, 
                         aliasing.print = FALSE), 
         which = "aliasing", which.criteria = c("aeff", "eeff", "xeff", "ord"))
   testthat::expect_true(!is.null(p$aliasing))
   testthat::expect_equal(nrow(p$aliasing), 2)
   testthat::expect_equal(p$aliasing$Alias, c("Line", "## Information remaining"))
   testthat::expect_true(all(abs(p$aliasing$aefficiency - c(0.2749455, 0.7250545)) < 1e-05))

   ### Investigate a factor-covariate interaction
   LRRCD.BlkxPosn.canon <- designAnatomy(formulae = list(carts = ~ Block/MainPosn:BLane/Cart,
                                                         trts = ~ Block*xMainPosn + Line*Watering),
                                         grandMean = TRUE, omit.projectors = "combined", data = LRRCD.lay)
   summary(LRRCD.BlkxPosn.canon, which.criteria = c("aeff", "eeff", "xeff", "ord", "dfor"))
   testthat::expect_true(!is.null(LRRCD.BlkxPosn.canon$aliasing))
   testthat::expect_equal(nrow(LRRCD.BlkxPosn.canon$aliasing), 3)
   testthat::expect_equal(LRRCD.BlkxPosn.canon$aliasing$Alias, 
                          c("xMainPosn", "Block[xMainPosn]", "## Information remaining"))
   testthat::expect_equal(LRRCD.BlkxPosn.canon$aliasing$Source, rep("Line", 3))
   testthat::expect_true(all(abs(LRRCD.BlkxPosn.canon$aliasing$aefficiency - 
                                    c(0.2749455, 0.3605941, 0.9698102)) < 1e-05))
   
   # Factor-covariate and factor-factor interaction
   print(p <- pstructure(~ Block*xMainPosn + Line*Watering, data = LRRCD.lay, 
                         aliasing.print = FALSE), 
         which = "aliasing", which.criteria = c("aeff", "eeff", "xeff", "ord"))
   testthat::expect_true(!is.null(p$aliasing))
   testthat::expect_equal(nrow(p$aliasing), 3)
   testthat::expect_equal(p$aliasing$Alias, c("xMainPosn", "Block[xMainPosn]", "## Information remaining"))
   testthat::expect_equal(p$aliasing$Source, rep("Line", 3))
   testthat::expect_true(all(abs(p$aliasing$aefficiency - c(0.2749455, 0.3605941, 0.9698102)) < 1e-05))
   
   print(p <- pstructure(formula = ~ Block*xMainPosn + Line*Watering, data = LRRCD.lay, 
                         aliasing.print = FALSE, grandMean = TRUE), 
         which = "aliasing", which.criteria = c("aeff", "eeff", "xeff", "ord"))
   testthat::expect_true(!is.null(p$aliasing))
   testthat::expect_equal(nrow(p$aliasing), 3)
   testthat::expect_equal(p$aliasing$Alias, c("xMainPosn", "Block[xMainPosn]", "## Information remaining"))
   testthat::expect_equal(p$aliasing$Source, rep("Line", 3))
   testthat::expect_true(all(abs(p$aliasing$aefficiency - c(0.2749455, 0.3605941, 0.9698102)) < 1e-05))
   
})

cat("#### Test for marginality when not all treatment combinations are observed\n")
test_that("TwoxTwo3cells", {
  skip_on_cran()
  library(dae)
  data("Exp720.des")
  
  print(tab <- with(Exp720.des, table(Treatment,Soil,Sterilized)), zero.print = ".")
  testthat::expect_equal(sum(tab != 0), 14)
   
  #Simplified script for a case when the interaction is aliased with the main effects because only 3 or 2x2 cells
  pstr <- pstructure( ~ Soil*(Sterilized/Microbe),
                      data = subset(Exp720.des, Soil != "YP" & Nematode == "no"))
  marg <- pstr$marginality
  testthat::expect_equal(nrow(marg), 4)
  testthat::expect_equal(ncol(marg), 4)
  testthat::expect_true(all.equal(marg[upper.tri(marg, diag = TRUE)], c(1,0,1,0,1,1,1,1,1,1)))
  
  #Two Soils and Nematode = "no" -  smaller problem to investigate
  testthat::expect_warning(
    Exp720C.canon.part <- designAnatomy(formulae = list(unit = ~ Block/MainUnit/Cart,
                                                        trt  = ~ Soil*(Sterilized/Microbe)),
                                        data = subset(Exp720.des, Soil != "YP" & Nematode == "no")),
    regexp = "Soil:Sterilized is aliased with previous terms in the formula and has been removed")
  testthat::expect_true(all.equal(Exp720C.canon.part$aliasing$Alias, 
                                  c("Soil", "## Information remaining", "Sterilized", "## Aliased")))
  summ <- summary(Exp720C.canon.part)
  testthat::expect_true(all.equal(summ$decomp$Source.trt[-1], 
                                  c("Soil","Residual","Sterilized","Microbe[Sterilized]",
                                    "Soil#Sterilized#Microbe","Residual")))
  testthat::expect_true(all.equal(summ$decomp$Source.unit, 
                                  c("Block", "MainUnit[Block]", "MainUnit[Block]", 
                                    rep("Cart[Block:MainUnit]", 4))))
})


cat("#### Test for marginality when nested treatments\n")
test_that("SprayerRates", {
  skip_on_cran()
  library(dae)
  b <- 3
  t <- 6
  #'## Construct a systematic layout
  RCBD.sys <- cbind(fac.gen(generate = list(Blocks=b, Plots=t)),
                    fac.gen(generate = list(Pressure = c("140", "330"), 
                                            Speed = c("3.6", "2.6", "1.8")), times = b))
  
  #'## Obtain the randomized layout
  RCBD.lay <- designRandomize(allocated         = RCBD.sys[c("Pressure", "Speed")], 
                              recipient         = RCBD.sys[c("Blocks", "Plots")], 
                              nested.recipients = list(Plots = "Blocks"),
                              seed              = 353441)
  #'## Add nested factors
  RCBD.lay <- within(RCBD.lay, 
                     {
                       Treatments <- fac.combine(list(Pressure, Speed), combine.levels = TRUE)
                       Rates <- fac.recast(Treatments, 
                                           newlevels = c("2090", "2930", "4120", 
                                                         "2930", "4120", "5770"))
                     })
  RCBD.lay <- with(RCBD.lay, cbind(RCBD.lay, 
                                   fac.multinested(nesting.fac = Rates, 
                                                   nested.fac  = Treatments,
                                                   fac.prefix  = "Rate")))
  RCBD.canon <- designAnatomy(formulae  = list(units = ~ Blocks/Plots, 
                                               trts  = ~ Rates/(Rate2090 + Rate2930 + Rate4120 +
                                                                  Rate5770)),
                              grandMean = TRUE, data = RCBD.lay)
  summ <- summary(RCBD.canon, which.criteria = "aeff")
  marg <- RCBD.canon$marginality$trts
  testthat::expect_equal(nrow(marg), 3)
  testthat::expect_equal(ncol(marg), 3)
  testthat::expect_true(all.equal(marg[upper.tri(marg, diag = TRUE)], c(1,1,1,1,0,1)))
  
  testthat::expect_true(all.equal(summ$decomp$Source.units, c("Mean", "Blocks", rep("Plots[Blocks]",4))))
  testthat::expect_true(all.equal(summ$decomp$Source.trts, 
                                  c("Mean",NA,"Rates","Rate2930[Rates]","Rate4120[Rates]","Residual")))
  testthat::expect_true(all.equal(summ$decomp$df2, c(1,NA,3,1,1,10)))
})
  