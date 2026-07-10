# Extracted from testInteractionABCplot.r:12

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "dae", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("interactionPlot")
cat("#### Test interaction.ABC.plot\n")

# test -------------------------------------------------------------------------
skip_on_cran()
library(dae)
data(ABC.Interact.dat)
testthat::expect_silent(plt <- interaction.ABC.plot(MOE, A, B, C, data=ABC.Interact.dat))
vdiffr::expect_doppelganger("ABCInteract-Plain", plt)
