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

