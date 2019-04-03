#devtools::test("dae")
context("analysis")

cat("#### Test for Tukey.1df using Kirk example\n")
test_that("Tukey1dfKirk", {
  skip_on_cran()
  library(dae)
  #'## Dataset from Kirk (Table 7.2-1)
  kirk721 <- data.frame(su = factor(rep(1:8,each=4)),
                        A = factor(rep(1:4,8)),
                        y = c(    3,4,4,3,    2,4,4,5,    2,3,3,6, 3,3,3,5,
                                  1,2,4,7,    3,3,6,6, 4,4,5,10,    6,5,5,8)
  )
  
  #'### Add a unit factor
  kirk721$unit <- kirk721$A
  
  #'## Test tukey.1df - su + A
  #'### No error term
  fit <- aov(y ~ su + A, kirk721)
  kirk1df <- tukey.1df(fit,kirk721,"Residuals")
  testthat::expect_lt(abs(kirk1df$Tukey.SS - 1.773891), 1e-5)
  testthat::expect_lt(abs(kirk1df$Tukey.F - 1.279581), 1e-5)
  testthat::expect_lt(abs(kirk1df$Tukey.p - 0.271392), 1e-4)
  testthat::expect_lt(abs(kirk1df$Devn.SS - 27.72611), 1e-5)
  
  fit <- aov(y ~ su + A + Error(su), kirk721)
  kirk1df <- tukey.1df(fit,kirk721)
  testthat::expect_lt(abs(kirk1df$Tukey.F - 1.279581), 1e-5)
  
  fit <- aov(y ~ su + A + Error(su/unit), kirk721)
  kirk1df <- tukey.1df(fit,kirk721, "su:unit")
  testthat::expect_lt(abs(kirk1df$Tukey.F - 1.279581), 1e-5)

  #'### su not outside Error 
  fit <- aov(y ~ A + Error(su), kirk721)
  kirk1df <- tukey.1df(fit,kirk721)
  testthat::expect_gt(abs(kirk1df$Tukey.F - 1.279581), 1e-01)

  fit <- aov(y ~ A + Error(su/unit), kirk721)
  kirk1df <- tukey.1df(fit,kirk721, "su:unit")
  testthat::expect_gt(abs(kirk1df$Tukey.F - 1.279581), 1e-01)
  
})
