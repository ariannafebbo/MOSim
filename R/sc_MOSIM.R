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






#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @return a named list with simulation parameters for each omic as values

param_estimation <- function(omic, cellTypes){
  
  N_omic <- length(omic)
  norm_list <- lapply(omic, scran_normalization)
  param_est_list <- list()
  
  for(i in 1:N_omic){

    param_est <- SPARSim_estimate_parameter_from_data(raw_data = omic[[i]],
                                                    norm_data = norm_list[[i]],
                                                    conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omic)[i])]] <- param_est
    
  }
  
  return(param_est_list)
  
}



#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @param numberCells vector of numbers. The numbers correspond to the number of cells the user wants to simulate per each cell type. The length of the vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be specified just if \code{numberCells} is specified.
#' @param output_sim_parameter boolean flag. If TRUE, the function will output a list of simulation parameter per each cell type. (Default: FALSE)
#' @return Or a list of count matrices. 1 per each omic. Or a list of Seurat obj. ?  


sc_MOSim <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL, output_sim_parameter = FALSE ){
  
  lapply(omics, param_estimation)

  sim_parameter_matrix <- NULL
  if(!is.null(sim_parameter)){
    cat("Generating simulation parameters matrices...", "\n")
    
    return(SPARSim_sim_param)
  }
}
