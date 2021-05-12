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


cat("#### Test for pstructure with factor nesting\n")
test_that("pstucture_fac.multinested", {
  skip_on_cran()
  library(dae)
  
  #'## Set constants
  nblks <- 6
  treat.levs <- c("Control","Dr","Na","LN")
  (ntreats <- length(treat.levs))
  lines.lev <- c("O. aust", "Calrose", paste0("Transgenic", 1:7))
  (nlines <- length(lines.lev))
  
  #'### Systematic allocation
  sys.lay <- cbind(
    fac.gen(list(Block = nblks, MainUnit = ntreats, Cart = nlines)),
    fac.gen(list(Treatment = treat.levs, Line = lines.lev), times = nblks))
  
  #'### Randomization
  rand.lay <- designRandomize(recipient = sys.lay[,1:3],
                              allocated = sys.lay[,4:5],
                              nested.recipients = list(MainUnit = "Block",
                                                       Cart = c("MainUnit", "Block")), 
                              seed = 82604)
  
  #'## Add nested factors
  #'### Line nested within Treatments
  rand.lay <- cbind(rand.lay, 
                    with(rand.lay, fac.multinested(nesting.fac = Treatment, nested.fac = Line, 
                                                   fac.prefix = "Line")))
  #'### Factors that remove contrast involving O. aust
  rand.lay <- within(rand.lay,
                     {
                       OaVsRest <- fac.uselogical(Line == "O. aust", labels = c("O. aust", "Other"))
                       OaTreat <- fac.recode(fac.combine(list(Line, Treatment)),
                                             c(levels(Treatment), rep("Other", 32)))
                     })
  #'### Factors for Lines within Treatments, excluding O. aust
  rand.lay <- within(rand.lay,
                     {
                       OaDr <- fac.uselogical(LineDr == "O. aust", labels = c("O. aust", "Other"))
                       OaControl <- fac.uselogical(LineControl == "O. aust", labels = c("O. aust", "Other"))
                       OaLN <- fac.uselogical(LineLN == "O. aust", labels = c("O. aust", "Other"))
                       OaNa <- fac.uselogical(LineNa == "O. aust", labels = c("O. aust", "Other"))
                     })
  
  #'## Investigate Treatment terms
  #'### Removal of O. aust from the Treatments*Line
  print(trt.str <- pstructure(~ OaVsRest/OaTreat + Treatment*Line, data = rand.lay), 
        which = "proj")
  testthat::expect_true(all(names(trt.str$Q) == c("OaVsRest", "OaTreat[OaVsRest]", "Treatment", 
                                                  "Line[OaVsRest]", "Treatment#Line[OaVsRest]")))
  testthat::expect_true(all(trt.str$aliasing$Source == c("Treatment", "Treatment", 
                                                         "Treatment#Line[OaVsRest]","Treatment#Line[OaVsRest]")))
  testthat::expect_true(all(trt.str$aliasing$Alias == c("OaTreat[OaVsRest]", 
                                                        "## Information remaining", "Treatment",
                                                        "## Information remaining")))
  
  #'### Removal of O. aust from remaining Lines nested within Treats
  print(trt.str <- pstructure(~ OaVsRest/OaTreat + Treatment/(LineControl + LineDr + 
                                                                LineLN + LineNa), 
                              which.criteria = c("aeff", "xeff", "eeff", "ord"), 
                              data = rand.lay), which = c("proj", "alias"))
  testthat::expect_true(all(names(trt.str$Q) == c("OaVsRest", "OaTreat[OaVsRest]", "Treatment", 
                                                  "LineControl[Treatment]", "LineDr[Treatment]", 
                                                  "LineLN[Treatment]", "LineNa[Treatment]")))
  testthat::expect_true(all(trt.str$aliasing$df == c(3,3, rep(c(1,1,1,7), times = 4))))
  testthat::expect_true(all(trt.str$aliasing$Alias[c(2,6,10,14,18)] == "## Information remaining"))
  
  #'### Treaments pooled over ALL lines but then separation of O. aust from remaining Lines, both nested within Treats
  print(trt.str <- pstructure(~ Treatment/(OaControl + LineControl + OaDr + LineDr +
                                             OaLN + LineLN + OaNa + LineNa), data = rand.lay),
        which = "proj")
  testthat::expect_true(all(names(trt.str$Q) == c("Treatment", "OaControl[Treatment]", 
                                                  "LineControl[Treatment:OaControl]", 
                                                  "OaDr[Treatment]", "LineDr[Treatment:OaDr]", 
                                                  "OaLN[Treatment]", "LineLN[Treatment:OaLN]", 
                                                  "OaNa[Treatment]", "LineNa[Treatment:OaNa]")))
  testthat::expect_true(all(trt.str$aliasing$df == c(3, rep(c(1,7), times = 4))))
  testthat::expect_true(is.null(trt.str$aliasing))
})

cat("#### Test for pstructure with generalized factors\n")
test_that("pstucture_genfac", {
  skip_on_cran()
  library(dae)

  pepalt.sys <- fac.gen(list(Rep = 2, Plate = 3, Side = 2, Boxrow = 2, Shelf = 4))
  pepalt.str <- pstructure( ~ (Shelf:Boxrow)*(Rep/(Side:Plate)), data = pepalt.sys)
  (sources <- pepalt.str$sources)
  testthat::expect_true(all(sources == c("Shelf:Boxrow", "Rep", "Side:Plate[Rep]", 
                                         "(Shelf:Boxrow)#Rep", "(Shelf:Boxrow)#(Side:Plate)[Rep]")))
  
  pepalt.str <- pstructure( ~ (Rep/Plate)*(Boxrow/(Shelf:Side)), data = pepalt.sys)
  (sources <- pepalt.str$sources)
  testthat::expect_true(all(sources == c("Rep", "Plate[Rep]", "Boxrow", "Shelf:Side[Boxrow]", 
                                         "Rep#Boxrow", "Rep#(Shelf:Side)[Boxrow]", 
                                         "Plate#Boxrow[Rep]", "Plate#(Shelf:Side)[Rep:Boxrow]")))
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
                                                             "E#R[T:M]", "P#R[T:M:E]")))
  
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
  