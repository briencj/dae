#devtools::test("dae")
context("designGGPlot")

cat("#### Test for designGGPlot using FHPain\n")
test_that("FHPain_designGGPlot", {
  skip_on_cran()
  library(dae)
  
  Pain.lay <- cbind(fac.gen(list(Expressiveness = 2, Patients = 4, Occasions = 2)),
                    fac.gen(list(Motions = c("active", "passive")), times = 8))
  
  #'## Plot the layout
  cell.colours <- c("lightblue","lightcoral","lightgoldenrod","lightgreen","lightgrey",
                    "lightpink","lightsalmon","lightcyan","lightyellow","lightseagreen",
                    "lightskyblue","lightslateblue","lightslategrey","lightsteelblue")
  #Test switch only for row facets
  plt1 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                      row.factors = c("Expressiveness", "Patients"), 
                      column.factors = "Occasions",
                      facetstrips.switch = "y", 
                      title = NULL, title.size = 20, axis.text.size = 20, 
                      blockdefinition = cbind(1,5),
                      ggplotFuncs = list(theme(strip.placement = "outside")))
  vdiffr::expect_doppelganger("rowFacets - switch", plt1)

  #Test switch and placement for row facets
  plt2 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                       row.factors = c("Expressiveness", "Patients"), 
                       column.factors = "Occasions", 
                       facetstrips.switch = "y", facetstrips.placement = "outside.title",
                       title = NULL, title.size = 20, axis.text.size = 20, 
                       blockdefinition = cbind(1,5))
  vdiffr::expect_doppelganger("rowFacets - switch & placement", plt2)
  
  #Test placement only for column facets
  plt3 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                       column.factors = c("Expressiveness", "Patients"), 
                       row.factors = "Occasions", facetstrips.placement = "outside.title",
                       title = NULL, title.size = 20, axis.text.size = 20, 
                       blockdefinition = cbind(1,5))
  vdiffr::expect_doppelganger("colFacets - placement", plt3)
  
  #Test x.axis.position only for column facets
  plt4 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                      column.factors = c("Expressiveness", "Patients"), 
                      row.factors = "Occasions", 
                      title = NULL, title.size = 20, axis.text.size = 20, 
                      blockdefinition = cbind(1,5), x.axis.position = "bottom")
  vdiffr::expect_doppelganger("colFacets - bottom", plt4)
  
  #Test switch and placement with x.axis.position for column facets
  plt5 <- designGGPlot(Pain.lay, labels = "Motions", label.size = 7, 
                       colour.values = cell.colours[2:3], 
                      column.factors = c("Expressiveness", "Patients"), 
                      row.factors = "Occasions", facetstrips.switch = "x",
                      facetstrips.placement = "outside.title",
                      title = NULL, title.size = 20, axis.text.size = 20, 
                      blockdefinition = cbind(1,5), x.axis.position = "bottom")
  vdiffr::expect_doppelganger("colFacets - bottom, switch & placement", plt5)
})


cat("#### Test for designGGPlot using SPLGrass\n")
test_that("SPLGrass_designGGPlot", {
  skip_on_cran()
  library(dae)
  
  data("SPLGrass.dat")
  
  tmp <- within(SPLGrass.dat, 
                {
                  Season <- fac.combine(list(Spring, Summer), combine.levels = "TRUE",
                                        sep = ",")
                  Treatment <- fac.combine(list(Period, Season), combine.levels = "TRUE",
                                           sep = "\n")
                })
  plt <- designGGPlot(tmp, labels = "Treatment", 
                      row.factors = c("Rows", "SubRows"), column.factors = c("Columns", "SubColumns"),
                      facetstrips.switch = "y", facetstrips.placement = "outside.title", 
                      cellfillcolour.column = "Period", cellalpha = 0.75, label.size = 6, 
                      blockdefinition = cbind(2,2))
  vdiffr::expect_doppelganger("Rows and Columns indexed by 2 factors", plt)
  
  plt2 <- designGGPlot(tmp, labels = "Treatment", 
                      row.factors = c("Rows", "SubRows"), column.factors = c("Columns", "SubColumns"),
                      facetstrips.switch = "y", 
                      cellfillcolour.column = "Period", cellalpha = 0.75, label.size = 6, 
                      blockdefinition = cbind(5,5), printPlot = FALSE)
  plt2 <- designBlocksGGPlot(plt2, nrows = 2, ncolumns = 2, blockdefinition = rbind(c(2,2)),
                             facetstrips.placement = "outside.title")
  vdiffr::expect_doppelganger("Using facetstrips.placement from designBlocksGGPlot", plt2)
})