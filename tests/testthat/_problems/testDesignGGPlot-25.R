# Extracted from testDesignGGPlot.r:25

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "dae", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("designGGPlot")
cat("#### Test for designGGPlot using FHPain\n")

# test -------------------------------------------------------------------------
skip_on_cran()
library(dae)
Pain.lay <- cbind(fac.gen(list(Expressiveness = 2, Patients = 4, Occasions = 2)),
                    fac.gen(list(Motions = c("active", "passive")), times = 8))
cell.colours <- c("lightblue","lightcoral","lightgoldenrod","lightgreen","lightgrey",
                    "lightpink","lightsalmon","lightcyan","lightyellow","lightseagreen",
                    "lightskyblue","lightslateblue","lightslategrey","lightsteelblue")
plt1 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                      row.factors = c("Expressiveness", "Patients"), 
                      column.factors = "Occasions",
                      facetstrips.switch = "y", 
                      title = NULL, title.size = 20, axis.text.size = 20, 
                      blockdefinition = cbind(1,5),
                      ggplotFuncs = list(theme(strip.placement = "outside")))
vdiffr::expect_doppelganger("rowFacets - switch", plt1)
