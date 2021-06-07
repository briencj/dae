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

