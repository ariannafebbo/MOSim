#Required packages
suppressPackageStartupMessages({
  library(SPARSim)
  library(dplyr)
  library(Seurat)
  library(Signac)
  library(stringr)
  })

#' @param omic A string which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as rows and cells as columns. Alternatively MOSim allows user to estimate the input parameters from an existing count table by typing 'example_matrix'
#' @return a named list with omics type as name and the count matrix as value


sc_omicData <- function(omics, data = NULL){
  
  if (omics != "scRNA-seq" && omics != "scATAC-seq"){
    
    print("omics must be a either 'scRNA-seq' or 'scATAC-seq'")
    return(NA)
    
  }
  
  if (is.null(data)){ 
    
      if (omics == "scRNA-seq"){ 
      
        ##scRNA##
        rna_orig_counts <- readRDS("../data/rna_orig_counts.rds")
        omics_list <- list("scRNA-seq" = rna_orig_counts)
        return(omics_list)
      
      } else if (omics =="scATAC-seq"){
      
        ##scATAC##
        atac_orig_counts <- readRDS("../data/atac_orig_counts.rds")
        omics_list <- list("scATAC-seq" = atac_orig_counts)
        return(omics_list)
      
      }
  } 
  
  if (! is.matrix(data)){
    
    print("data must be a matrix")
    return(NA)
    
  } else if (is.matrix(data)){
    
    print(omics)
    omics_list <- list()
    omics_list[[omics]] <- data 
    return(omics_list)
    
  }
}






#' @param omics named list containing the omics to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @return a named list with simulation parameters for each omics as values

param_estimation <- function(omics, cellTypes){
  
  N_omics <- length(omics)
  norm_list <- lapply(omics, scran_normalization)
  param_est_list <- list()
  
  for(i in 1:N_omics){

    param_est <- SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                    norm_data = norm_list[[i]],
                                                    conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
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

  
  if(missing(numberCells) && missing(mean) && missing(sd)){
    
    param_list <- param_estimation(omics, cellTypes)
    
    N_param <- length(param_list)
    sim_list <- list()
    
    for(i in 1:N_param){
      
      sim <- SPARSim_simulation(dataset_parameter = param_list[[i]])
      sim <- sim[["count_matrix"]]
      sim_list[[paste0("sim_", names(omics)[i])]] <- sim
      
    }
    
    seu_obj <- list()
    N_sim <- length(sim_list)
    
    for(i in 1:N_sim)
    seu_obj <- CreateSeuratObject(counts = sim_list[[i]], assay = )
    
  }
  # lapply(omics, param_estimation)
  # 
  # sim_parameter_matrix <- NULL
  # if(!is.null(sim_parameter)){
  #   cat("Generating simulation parameters matrices...", "\n")
  #   
  #   return(SPARSim_sim_param)
  # }
}

param <- param_estimation(omic_list, conditions)
sim <- SPARSim_simulation(param[[1]])
sim[[1]]

param[[1]]

names(prova)[1]
sub('sim_sc','',names(prova)[1])

