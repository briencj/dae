#devtools::test("dae")
context("canonical")

cat("#### Test for  projector\n")
test_that("projector", {
  skip_on_cran()
  library(dae)
  
  #Tests for valid projector
  s <- matrix(c(1,2,2,1), nrow=2)
  testthat::expect_error(proj <- new("projector", .Data=s)) #not idempotent
  testthat::expect_error(proj.s <- projector(s)) #not idempotent
  testthat::expect_error(new("projector", .Data = matrix(1:6,2,3))) #not square

  #Tests for degfree
  m <- matrix(rep(0.5,4), nrow=2)
  testthat::expect_equal(degfree(projector(m)), 1)
  proj.m <- new("projector", data=m)
  testthat::expect_true(is.na(degfree(proj.m)))
  testthat::expect_equal(degfree(projector(proj.m)), 1)
  
  degfree(proj.m) <- 6
  testthat::expect_equal(degfree(proj.m), 6)
  degfree(proj.m) <- as.matrix(proj.m)
  testthat::expect_equal(degfree(proj.m), 1)
  degfree(proj.m) <- m
  testthat::expect_true(correct.degfree(proj.m))
  
  proj.m@.Data <- s
  testthat::expect_error(validObject(proj.m))
  testthat::expect_error(is.projector(proj.m))
  
  #Tests for projectors
  a <- factor(rep(3, 4))
  M.a <- fac.meanop(a)
  testthat::expect_true(is.projector(M.a))
  
  A <- factor(rep(1:2, each=6))
  B <- factor(rep(1:3, each=2, times=2))
  AB <- fac.combine(list(A,B))
  M.AB <- fac.meanop(AB)
  
  testthat::expect_true(validObject(M.AB))
  testthat::expect_equal(degfree(M.AB), 6)
 
})