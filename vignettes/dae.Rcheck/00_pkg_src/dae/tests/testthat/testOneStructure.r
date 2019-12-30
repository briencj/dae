#devtools::test("dae")
context("analysis")

cat("#### Test for designAnatomy with single structure\n")
test_that("OneStructure", {
  skip_on_cran()
  library(dae)
  #'### Make a Latin square
  ls.ran <- designRandomize(allocated = data.frame(Trt = factor(designLatinSqrSys(7))),
                            recipient = list(Row = 7, Column = 7), 
                            seed = 354131)
  
  lsadd.canon <- designAnatomy(list(plt = ~ Row+Column, trt = ~ Trt), data = ls.ran)
  summadd <- summary(lsadd.canon)
  testthat::expect_warning(print(summadd))
  testthat::expect_equal(length(summadd),2)
  testthat::expect_true(all(summadd$decomp$Source.plt == c("Row", "Column")))
  testthat::expect_true(all(summadd$decomp$df1 == 6))
  testthat::expect_true(all(is.na(summadd$decomp$Source.trt)))
  testthat::expect_true(all(is.na(summadd$decomp$df2)))
  
  ls.canon <- designAnatomy(list(plt = ~ Row*Column, trt = ~ Trt), data = ls.ran)
  summ <- summary(ls.canon)
  testthat::expect_equal(attr(summ$decomp, which = "n"), 49)
  testthat::expect_equal(length(summ),2)
  testthat::expect_true(all(summ$decomp$Source.plt == c("Row", "Column", "Row#Column", "Row#Column")))
  testthat::expect_true(all(summ$decomp$df1 == c(6,6,36,36)))
  testthat::expect_true(all(summ$decomp$Source.trt[3:4] == c("Trt", "Residual")))
  testthat::expect_true(all(summ$decomp$df2[3:4] == c(6,30)))
  
  ls1.canon <- designAnatomy(list(plt = ~ Row+Column), data = ls.ran)
  summ1 <- summary(ls1.canon)
  testthat::expect_equal(length(summ1),2)
  testthat::expect_true(all(summ1$decomp$Source.plt == c("Row", "Column")))
  testthat::expect_true(all(summ1$decomp$df == 6))
  
  struct <- pstructure(~ Row+Column, data = ls.ran)
  
})
