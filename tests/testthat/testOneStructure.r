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


cat("#### Test for pstructure with difficult marginalitysingle structure\n")
test_that("PlaidInteractions", {
  skip_on_cran()
  library(dae)
  # Generate first-phase sytematic design
  ph1.sys <- cbind(fac.gen(list(Expressive = c("Yes", "No"), Patients = 4, Occasions = 2)),
                   fac.gen(list(Motions = c("active", "passive")), times = 8))
  
  # Generate the two-phase systematic design
  ph2.sys <- cbind(fac.gen(list(Raters = 74, Viewings = 16)),
                   fac.gen(list(Trainings = 2, 16), times = 37),
                   rep.data.frame(ph1.sys, times =74))
  
  # Randomize the two-phase design
  ph2.lay <- designRandomize(allocated = ph2.sys[c("Trainings", "Expressive", "Patients",
                                                   "Occasions", "Motions")],
                             recipient = ph2.sys[c("Raters", "Viewings")],
                             except = "Viewings",
                             seed = 15674)
  
  # Convert names of the factors to single capital letters
  ph2.L.lay <- ph2.lay
  names(ph2.L.lay)[match(c("Raters", "Viewings", "Trainings", "Expressive", "Patients", 
                           "Occasions", "Motions"), names(ph2.L.lay))] <- c("R", "V", "T", 
                                                                            "E", "P", "O", "M")
  
  #Test the neat formula
  terms <- attr(terms(~ T * M * E + T:M:E:P + R:(M * (E / P)), data = ph2.L.lay), 
                which = "term.labels")
  testthat::expect_equal(length(terms), 13)

  alloc.canon <- designAnatomy(list(alloc = ~ T * M * E + T:M:E:P + R:(M * (E / P))), 
                               data = ph2.L.lay)
  testthat::expect_true(all(alloc.canon$terms$alloc %in% terms))
  testthat::expect_true(all(names(alloc.canon$sources$alloc) %in% terms))
  testthat::expect_true(all(alloc.canon$sources$alloc %in% c("T", "M", "T#M", "E", "T#E", 
                                                             "M#E", "T#M#E", "P[T:M:E]", 
                                                             "R[T:M]", "R[T:E]", "P[T:E:R]", 
                                                             "M#E#R[T]", "M#P#R[T:E]")))
  
  #Test the simple formula
  terms <- attr(terms(~ (T + R) * M * (E / P), keep.order = TRUE, data = ph2.L.lay), 
                which = "term.labels")
  testthat::expect_equal(length(terms), 17)
  alloc.canon <- designAnatomy(list(alloc = ~ (T + R) * M * (E / P)), 
                               data = ph2.L.lay)
  testthat::expect_true(all(alloc.canon$terms$alloc %in% terms))
  testthat::expect_true(all(names(alloc.canon$sources$alloc) %in% terms))
  testthat::expect_true(all(alloc.canon$sources$alloc %in% c("T", "R[T]", "M", "T#M", "R#M[T]", 
                                                             "E", "P[E]", "T#E", "T#P[E]", 
                                                             "R#E[T]", "R#P[T:E]",
                                                             "M#E", "M#P[E]", "T#M#E", 
                                                             "T#M#P[E]", "R#M#E[T]", 
                                                             "R#M#P[T:E]")))
}) 
  