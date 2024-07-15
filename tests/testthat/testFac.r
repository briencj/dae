#devtools::test("dae")
context("Factors")

cat("#### Test for factor manipulation\n")
test_that("fac.recast", {
  skip_on_cran()
  library(dae)
  #'# Script to investigate the reduced design from Smith et al. (2015)
  
  #'## Set up designs
  Trt <- factor(rep(letters[1:5], 4))
  
  A <- fac.recode(Trt, c(3,1,2,2,1))
  A <- fac.recode(Trt, letters[c(3,1,2,2,1)])  
  A <- fac.recode(Trt, letters[c(3,1,2,2,1)], labels = letters[c(3,1,2)])  
  #You cannot reorder the levels vector without affecting the levels stored with fac.recode
  
  Treats <- factor(rep(1:4, 4), labels=letters[1:4])
  
  #### The following change the values stored in the factor vector
  #reduce the levels from 4 to 3 and without re-ordering the levels vector
  A <- fac.recast(Treats, newlevels = letters[c(3,1,1,2)])  
  testthat::expect_true(all(A[1:4] == letters[c(3,1,1,2)]))
  testthat::expect_true(all(levels(A) == letters[1:3]))
  
  #reassign the values in the factor vector using fac.recast without re-ordering the levels attribute
  A <- fac.recast(Treats, newlevels = letters[4:1])  
  testthat::expect_true(all(A[1:4] == letters[4:1]))
  testthat::expect_true(all(levels(A) == letters[1:4]))
  
  #reduce the levels from 4 to 3, with re-ordering the levels vector
  A <- fac.recast(Treats, letters[c(3,1,1,2)], levels.order = letters[c(3:1)])  
  testthat::expect_true(all(A[1:4] == letters[c(3,1,1,2)]))
  testthat::expect_true(all(levels(A) == letters[3:1]))
  
  #reassign the values in the factor vector using fac.recast, but with re-ordering the levels attribute
  A <- fac.recast(Treats, newlabels = letters[4:1])
  testthat::expect_true(all(A[1:4] == letters[4:1]))
  testthat::expect_true(all(levels(A) == letters[4:1]))
  
  #reassigns the values in the factor vector, without changing the levels attribute 
  A <- fac.recast(Treats, newlevels = LETTERS[4:1], newlabels = letters[1:4])
  testthat::expect_true(all(A[1:4] == letters[4:1]))
  testthat::expect_true(all(levels(A) == letters[1:4]))
  A <- fac.recast(Treats, newlevels = 4:1, newlabels = levels(Treats))
  testthat::expect_true(all(A[1:4] == letters[4:1]))
  testthat::expect_true(all(levels(A) == letters[1:4]))
  
  #change the levels vector and the relationship between the values in the factor vector and the levels attribute
  A <- fac.recast(Treats, levels.order = letters[4:1], newlabels = LETTERS[1:4])
  testthat::expect_true(all(A[1:4] == LETTERS[4:1]))
  testthat::expect_true(all(levels(A) == LETTERS[1:4]))
  
  #### The following change the levels attribute without changing the values in the factor vector
  #Change the levels attribute without affecting the values in the factor vector
  A <- fac.recast(Treats, levels.order = letters[4:1])
  testthat::expect_true(all(A[1:4] == letters[1:4]))
  testthat::expect_true(all(levels(A) == letters[4:1]))
  
})  

cat("#### Test for conversion to numeric\n")
test_that("as.numfac", {
  skip_on_cran()
  library(dae)
  #'# Script to investigate the reduced design from Smith et al. (2015)
  
  #Test factor with unequal replication of levels
  Trt <- factor(rep(1:4, c(4,3,2,3)))
  x <- as.numfac(Trt, center = TRUE)
  testthat::expect_true(all(unique(x) == c(-1.5,-0.5,0.5,1.5)))
  testthat::expect_equal(attr(x, which = "center"), 2.5)
  x <- as.numfac(Trt, center = TRUE, scale = TRUE)
  testthat::expect_equal(attr(x, which = "center"), 2.5)
  testthat::expect_true(abs(attr(x, which = "scale") - 1.290994) < 1e-04)
  x <- as.numfac(Trt, scale = TRUE)
  testthat::expect_true(abs(attr(x, which = "scale") - 3.162278) < 1e-04)
  
  #Test the supply of a numeric center and scale
  x <- as.numfac(Trt, center = 2)
  testthat::expect_true(all(unique(x) == c(-1,0,1,2)))
  testthat::expect_equal(attr(x, which = "center"), 2)
  x <- as.numfac(Trt, center = 2, scale = 2)
  testthat::expect_true(all(unique(x) == c(-0.5,0,0.5,1)))
  testthat::expect_equal(attr(x, which = "center"), 2)
  testthat::expect_equal(attr(x, which = "scale"), 2)
  x <- as.numfac(Trt, scale = 2)
  testthat::expect_true(all(unique(x) == c(0.5,1,1.5,2)))
  testthat::expect_equal(attr(x, which = "scale"), 2)
  x <- as.numfac(Trt, center = TRUE, scale = 2)
  testthat::expect_true(all(unique(x) == c(-0.75,-0.25,0.25,0.75)))
  testthat::expect_equal(attr(x, which = "center"), 2.5)
  testthat::expect_equal(attr(x, which = "scale"), 2)
  x <- as.numfac(Trt, center = 2, scale = TRUE)
  testthat::expect_true(all(abs(unique(x) - c(-0.7071068,-0,0.7071068,1.4142136  )) < 1e-04))
  testthat::expect_equal(attr(x, which = "center"), 2)
  testthat::expect_true(abs(attr(x, which = "scale") - 1.414214) < 1e-04)
  
  #Test when a numeric is supplied
  z <- rep(1:5, 3)
  x <- as.numfac(z, center = TRUE)
  testthat::expect_equal(attr(x, which = "center"), 3)
  x <- as.numfac(z, center = TRUE, scale = TRUE)
  testthat::expect_equal(attr(x, which = "center"), 3)
  testthat::expect_equal(sum(x^2), 12)
  testthat::expect_true(abs(attr(x, which = "scale") - 1.581139) < 1e-04)
  x <- as.numfac(z, scale = TRUE)
  testthat::expect_true(abs(attr(x, which = "scale") - 3.708099) < 1e-04)
})  

