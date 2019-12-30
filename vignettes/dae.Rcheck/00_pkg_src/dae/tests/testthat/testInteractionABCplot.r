#devtools::test("dae")
context("interactionPlot")

cat("#### Test interaction.ABC.plot\n")
test_that("interaction.ABC.plot", {
  skip_on_cran()
  library(dae)
  #'# PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 2nd edn Wiley, New York"
 
  data(ABC.Interact.dat)
  testthat::expect_silent(interaction.ABC.plot(MOE, A, B, C, data=ABC.Interact.dat))
  
  testthat::expect_silent(interaction.ABC.plot(MOE, A, B, C, 
                                               xlab = "Factor A", 
                                               ylab = "M.O.E.", data=ABC.Interact.dat))
  
  
  ABC.Interact.dat$se <- rep(c(0.5,1), each=4)
  testthat::expect_silent(
    interaction.ABC.plot(MOE, A, B, C, data=ABC.Interact.dat,
                         ggplotFunc=list(geom_errorbar(data=ABC.Interact.dat, 
                                                       aes(ymax=MOE+se, ymin=MOE-se), 
                                                       width=0.2))))

})

