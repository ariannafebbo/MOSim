#' @param omic A string which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as rows and cells as columns. Alternatively MOSim allows user to estimate the input parameters from an existing count table by typing 'example_matrix'
#' @return a named list with omic type as name and the count matrix as value

#Required packages
suppressPackageStartupMessages({
  library(SPARSim)
  library(dplyr)
  library(Seurat)
  library(Signac)
  })

sc_omicData <- function(omic, data = NULL){
  
  if (omic != "scRNA-seq" && omic != "scATAC-seq"){
    
    print("omic must be a either 'scRNA-seq' or 'scATAC-seq'")
    return(NA)
  }
  
  
  if (is.null(data)){ 
    
      if (omic == "scRNA-seq"){ 
      
        ##scRNA##
        rna_orig_counts <- readRDS("/home/arifebbo/Desktop/MOSim/data/rna_orig_counts.rds")
      
        omic_list <- list("scRNA-seq" = rna_orig_counts)
        return(omic_list)
      
      } else if (omic =="scATAC-seq"){
      
        ##scATAC##
        atac_orig_counts <- readRDS("/home/arifebbo/Desktop/MOSim/data/atac_orig_counts.rds")
      
        omic_list <- list("scATAC-seq" = atac_orig_counts)
        return(omic_list)
      
      }
  } 
  
  if (! is.matrix(data)){
    
    print("data must be a matrix")
    return(NA)
    
  } else if (is.matrix(data)){
    print(omic)
    omic_list <- list()
    omic_list[[omic]] <- data 
    return(omic_list)
  }
}




#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 



sc_MOSim <- function(omics, numberCellTypes, numberCells = FALSE, mean = FALSE, sd = FALSE, sim_parameter = FALSE ){
}
