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
        rna_orig_counts <- readRDS("../data/rna_orig_counts.rds")
      
        omic_list <- list("scRNA-seq" = rna_orig_counts)
        return(omic_list)
      
      } else if (omic =="scATAC-seq"){
      
        ##scATAC##
        atac_orig_counts <- readRDS("../data/atac_orig_counts.rds")
      
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

scRNA <- sc_omicData("scRNA-seq")
scATAC <- sc_omicData("scATAC-seq")
omic_list <- c(scRNA, scATAC)



#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.

param_estimation <- function(omic, cellTypes){
  N_omic <- length(omic)
  
  lapply(omic, scran_normalization)
  #for(i in 1:N_omic){
    #param_est <- SPARSim_estimate_parameter_from_data(raw_data = omic[i], 
    #                                                norm_data = omic_norm, 
    #                                                conditions = cellTypes)
    #return(param_est)
    #}
}

prova <- param_estimation(omic_list, conditions)

omic_list[1]
prova <- param_estimation()

conditions <- list(cellA = c(1:160), cellB = c(161:270))



sc_MOSim <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL, sim_parameter = NULL ){
  
  lapply(omics, param_estimation)

  sim_parameter_matrix <- NULL
  if(!is.null(sim_parameter)){
    cat("Generating simulation parameters matrices...", "\n")
    
    return(SPARSim_sim_param)
  }
}
