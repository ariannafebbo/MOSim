library(testthat)

test_that("Passing a wrong string in 'omics' returns NA", {
  expect_message(sc_omicData("scR-seq", rna_orig_counts), NA)
})

test_that("Passing an object which is not a matrix in datain returns NA", {
  vector <- c(1,2,3)
  expect_message(sc_omicData("scATAC-seq", vector), NA)
})

test_that("Passing 'scATAC-seq' as omic returns list", {
  expect_type(sc_omicData("scATAC-seq"), "list")
})

# TESTS
#test1 <- sc_omicData("scR-seq", rna_orig_counts)
# [1] "omic must be a either 'scRNA-seq' or 'scATAC-seq'"

#vector <- c(1,2,3)
# test2 <- sc_omicData("scATAC-seq", vector)
# [1] "data must be a matrix"

#test3 <- sc_omicData("scATAC-seq")
#it generates a named list with scATAC-seq as name, and a count matrix with peaks as rows as values of the list




test_that("param_estimation returns a list", {
  scRNA <- sc_omicData("scRNA-seq")
  scATAC <- sc_omicData("scATAC-seq")
  omic_list <- c(scRNA, scATAC)
  conditions <- list(cellA = c(1:160), cellB = c(161:270))
  expect_type(param_estimation(omic_list, conditions),"list")
})
