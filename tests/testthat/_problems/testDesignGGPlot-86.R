# Extracted from testDesignGGPlot.r:86

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "dae", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("designGGPlot")
cat("#### Test for designGGPlot using FHPain\n")
cat("#### Test for designGGPlot using SPLGrass\n")

# test -------------------------------------------------------------------------
skip_on_cran()
library(dae)
data("SPLGrass.dat")
tmp <- within(SPLGrass.dat, 
                {
                  Season <- fac.combine(list(Spring, Summer), combine.levels = "TRUE",
                                        sep = ",")
                  Treatment <- fac.combine(list(Period, Season), combine.levels = "TRUE",
                                           sep = "\n")
                })
plt <- designGGPlot(tmp, labels = "Treatment", 
                      row.factors = c("Rows", "SubRows"), column.factors = c("Columns", "SubColumns"),
                      facetstrips.switch = "y", facetstrips.placement = "outside.title", 
                      cellfillcolour.column = "Period", cellalpha = 0.75, label.size = 6, 
                      blockdefinition = cbind(2,2))
vdiffr::expect_doppelganger("Rows and Columns indexed by 2 factors", plt)
