#devtools::test("dae")
context("factor")

cat("#### Test for fac.nested\n")
test_that("fac.nested", {
  skip_on_cran()
  library(dae)
  
  #Set up a data.frame with two factors A & B and use fac.nested to get B
  lay <- data.frame(A = factor(rep(c(1:3), c(3,6,4)), labels = letters[1:3]))
  lay$B <-fac.nested(lay$A)
  testthat::expect_equal(length(levels(lay$B)),6)
  
  #Test for wqhen NAs are present
  A <- factor(rep(c(1:3), c(3,6,4)), labels = letters[1:3])
  A[c(4,9)] <- NA
  A <- c(A, NA)
  B <-fac.nested(A)
  testthat::expect_equal(length(levels(B)),4)
  testthat::expect_true(all(is.na(B[c(4,9,14)])))
  
  #Test when the number of levels of the nesting.fac is large (33478)
  data("TagDay")
  Watering <- fac.nested(TagDay, nested.labs = c("","a"))
  testthat::expect_equal(length(levels(Watering)),2)
  testthat::expect_true(all(table(Watering) == c(33478, 523)))
  
})
cat("#### Test for fac.multinested\n")
test_that("Multiple nesting", {
  skip_on_cran()
  library(dae)

  #Set up a data.frame with two factors A & B and use fac.nested to get B
  lay <- data.frame(A = factor(rep(c(1:3), c(3,6,4)), labels = letters[1:3]))
  lay$B <-fac.nested(lay$A)
  testthat::expect_equal(length(levels(lay$B)),6)
  
  #Add factors for B within each level of A
  lay2 <- cbind(lay, fac.multinested(lay$A))
  testthat::expect_true(all(letters[1:3] %in% names(lay2)))
  testthat::expect_equal(length(levels(lay2$b)),7)
  testthat::expect_equal(levels(lay2$b)[1],"rest")
  testthat::expect_true(all(lay2$a[1:4] == c("1", "2", "3", "rest")))
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summ <- summary(canon2)
  testthat::expect_true(all(c("A", "a[A]", "b[A]", "c[A]") == summ$decomp$Source))
  testthat::expect_true(all(c(2,2,5,3) == summ$decomp$df))
  
  #Add factors for B within each level of A, but with levels and outlabel given
  lay2 <- cbind(lay, fac.multinested(lay$A, nested.levs = seq(10,60,10), outlabel = "other"))
  testthat::expect_true(all(letters[1:3] %in% names(lay2)))
  testthat::expect_equal(length(levels(lay2$b)),7)
  testthat::expect_true(all(lay2$a[1:4] == c("10", "20", "30", "other")))
  canon2 <- designAnatomy(list(~A/(a+b+c)), data = lay2)
  summ <- summary(canon2)
  testthat::expect_true(all(c("A", "a[A]", "b[A]", "c[A]") == summ$decomp$Source))
  testthat::expect_true(all(c(2,2,5,3) == summ$decomp$df))
  
  #Set a value of A to missing
  lay2 <- lay
  lay2$A[7] <- NA
  lay2 <- cbind(lay2, fac.multinested(lay2$A, outlabel = "0"))
  testthat::expect_true(all(letters[1:3] %in% names(lay2)))
  testthat::expect_equal(length(levels(lay2$b)),6)
  testthat::expect_equal(levels(lay2$b)[1],"0")
  
  #Replicate the combinations of A and B three times and index them with the factor sample
  lay3 <- rbind(lay,lay,lay)
  lay3$sample <- with(lay3, fac.nested(fac.combine(list(A,B))))
  
  #Add factors for B within each level of A
  lay4 <- cbind(lay3, fac.multinested(nesting.fac = lay$A, nested.fac = lay$B))
  testthat::expect_true(all(letters[1:3] %in% names(lay4)))
  testthat::expect_equal(length(levels(lay4$b)),7)
  testthat::expect_equal(levels(lay4$b)[1],"rest")
  canon4 <- designAnatomy(list(~(A/(a+b+c))/sample), data = lay4)
  summ <- summary(canon4)
  testthat::expect_true(all(c("A", "a[A]", "b[A]", "c[A]", "a#b#c#sample[A]") == summ$decomp$Source))
  testthat::expect_true(all(c(2,2,5,3,26) == summ$decomp$df))
  
  #Add factors for sample within each combination of A and B
  lay5 <- with(lay4, cbind(lay4, 
                           fac.multinested(nesting.fac = a, fac.prefix = "a"),
                           fac.multinested(nesting.fac = b, fac.prefix = "b"),
                           fac.multinested(nesting.fac = c, fac.prefix = "c")))
  testthat::expect_equal(ncol(lay5),19)
  testthat::expect_true(all(unlist(lapply(lay5[paste0("b", 1:6)], 
                                          function(fac) length(levels(fac)))) == rep(4, 6)))
  testthat::expect_equal(levels(lay5$b)[1],"rest")
  testthat::expect_equal(levels(lay5$b1)[1],"rest")
  canon5 <- designAnatomy(list(~A/(a/(a1+a2+a3)+b/(b1+b2+b3+b4+b5+b6)+c/(c1+c2+c3))), data = lay5)
  summ <- summary(canon5)
  testthat::expect_true(all(rep(c(2,5,2,3,2), c(5,1,6,1,3)) == summ$decomp$` df`))
  
  #Add factors for sample within each level of A
  lay6 <- cbind(lay4, 
                fac.multinested(nesting.fac = lay4$A, nested.fac = lay$sample, fac.prefix = "samp"))
  testthat::expect_true(all(c(letters[1:3], paste0("samp",letters[1:3])) %in% names(lay6)))
  testthat::expect_equal(length(levels(lay6$b)),7)
  testthat::expect_equal(length(levels(lay6$sampb)),19)
  testthat::expect_equal(levels(lay6$b)[1],"rest")
  testthat::expect_equal(levels(lay6$sampb)[1],"rest")
  canon6 <- designAnatomy(list(~A/(a/sampa+b/sampb+c/sampc)), data = lay6)
  summ <- summary(canon6)
  testthat::expect_true(all(c("A", "a[A]", "sampa[A:a]", "b[A]", "sampb[A:b]", "c[A]", "sampc[A:c]") == 
                              summ$decomp$Source))
  testthat::expect_true(all(c(2,2,6,5,12,3,8) == summ$decomp$df))
})

