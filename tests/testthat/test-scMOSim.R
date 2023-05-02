library(testthat)

# sc_omicData function
test_that("Passing a wrong string in 'omics' returns NA", {
  expect_message(sc_omicData(c("scR-seq"), rna_orig_counts), NA)
})

test_that("Passing an object in data which is neither a matrix or Seurat obj returns NA", {
  vector <- c(1,2,3)
  expect_message(sc_omicData(c("scATAC-seq"), vector), NA)
})

test_that("Passing 'scATAC-seq' as omic returns the expected subarray", {
  res<-sc_omicData("scATAC-seq")
  expected<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  expect_equal(res$`scATAC-seq`[1:50], expected)
})

test_that("Passing 'scRNA-seq' as omic and a Seurat obj as data returns a list", {
  scRNA <- sc_omicData("scRNA-seq")
  count <- scRNA[["scRNA-seq"]]
  Seurat_obj <- CreateSeuratObject(counts = count, assay = 'RNA')
  expect_type(sc_omicData(c("scRNA-seq"), c(Seurat_obj)), "list")
})

test_that("Passing an array ('scRNA-seq','scATAC-seq') as omic returns a list of length 2", {
  res<-sc_omicData(c("scRNA-seq","scATAC-seq"))
  expect_equal(length(res), 2)
})

#param_estimation function
test_that("param_estimation returns a list", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_type(param_estimation(omic_list, conditions, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2)),"list")
})

test_that("Not passing all optional arguments at once returns NA", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  conditions <- list(cellA = c(1:20), cellB = c(161:191))
  expect_message(param_estimation(omic_list, conditions, numberCells = c(10,20),sd = c(10^3, 10^2)), NA)
})


#sc_MOSim function
test_that("sc_MOSim returns a list with S4 obj as values", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-sc_MOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  expect_type(sim[[1]], "S4")
  
})

#sc_omicSim function
test_that("sc_omicSim returns a list", {
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-sc_MOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list(cellA= c(1:10), cellB = c(11:30))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_type(integration, "list")
})

test_that("sc_omicSim returns a list", {
  expected <- c("activator","NE","NE","NE","activator","activator","activator","NE","NE","NE")
  omic_list <- sc_omicData(c("scRNA-seq","scATAC-seq"))
  cell_types <- list(cellA = c(1:20), cellB = c(161:191))
  sim <-sc_MOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
  cell_types <- list(cellA= c(1:10), cellB = c(11:30))
  integration <- sc_omicSim(sim, cell_types, totalFeatures = 500)
  expect_equal(integration[["markers_cellA_cellB"]]$activity[1:10], expected)
})
