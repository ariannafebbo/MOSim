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
  
  if (! is.matrix(data) && class(data) != "Seurat"){
    
    print("data must be a matrix")
    return(NA)
    
  } else if (is.matrix(data)){
   
    omics_list <- list()
    omics_list[[omics]] <- data 
    return(omics_list)
    
  } else if (class(data) == "Seurat" && omics == "scRNA-seq"){
    
    omics_list <- list()
    counts <- data@assays[["RNA"]]@counts 
    omics_list[[omics]] <- counts
    return(omics_list)
    
  } else if (class(data) == "Seurat" && omics == "scATAC-seq"){
    
    omics_list <- list()
    counts <- data@assays[["ATAC"]]@counts 
    omics_list[[omics]] <- counts
    return(omics_list)
    
  }
}



#' @param omics named list containing the omics to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @return a named list with simulation parameters for each omics as values

param_estimation <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
  
  all_missing <- missing(numberCells) && missing(mean) && missing(sd)
  all_specified <- !missing(numberCells) && !missing(mean) && !missing(sd)
  
  if( !(all_missing || all_specified )){
  
    print("the user must either not provide the optional arguments or provide them all")
    return(NA)
    
  }
  
  N_omics <- length(omics)
  norm_list <- lapply(omics, scran_normalization)
  param_est_list <- list()
  
  for(i in 1:N_omics){

    param_est <- SPARSim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                    norm_data = norm_list[[i]],
                                                    conditions = cellTypes)
    param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
    
  }
  
  if(all_missing){
    
    return(param_est_list)
    
  } else if (all_specified){
    
      N_param_est_list<- length(param_est_list)
      N_cellTypes <- length(cellTypes)
      param_est_list_mod <- list()
    
      for(i in 1:N_param_est_list){
        cell_type_list <- list()
        
        for(j in 1:N_cellTypes){
        
        cond_param <- SPARSim_create_simulation_parameter(intensity = param_est_list[[i]][[j]][["intensity"]],
                                                          variability = param_est_list[[i]][[j]][["variability"]],
                                                          library_size = round(rnorm(n = numberCells[j], mean = mean[j], sd = sd[j])),
                                                          condition_name = param_est_list[[i]][[j]][["name"]],
                                                          feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
        cell_type_list[[names(cellTypes)[j]]] <- cond_param 
        
        }
        param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list 
      }
      
      return(param_est_list_mod)
  } 
  
}


scRNA <- sc_omicData("scRNA-seq")
scTATC <- sc_omicData("scATAC-seq")
omic_list <- c(scRNA, scTATC)
cell_types <- list(cellA = c(1:30), cellB = c(161:191))

prova_param_est_3 <- param_estimation(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
prova_param_est_4 <- param_estimation(omic_list, cell_types)



#' @param omics named list containing the omic to simulate as names, which can be "scRNA-seq" or "scATAC-seq, and the input count matrix as 
#' @param cellTypes list where the i-th element of the list contains the column indices for i-th experimental conditions. List must be a named list.
#' @param numberCells vector of numbers. The numbers correspond to the number of cells the user wants to simulate per each cell type. The length of the vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean per each cell type. Must be specified just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be specified just if \code{numberCells} is specified.
#' @return a list of Seurat object, one per each omic. 




sc_MOSim <- function(omics, cellTypes, numberCells = NULL, mean = NULL, sd = NULL){
    
    param_list <- param_estimation(omics, cellTypes, numberCells, mean, sd)
    
    N_param <- length(param_list)
    sim_list <- list()
    
    for(i in 1:N_param){
      
      sim <- SPARSim_simulation(dataset_parameter = param_list[[i]])
      sim <- sim[["count_matrix"]]
      sim_list[[paste0("sim_", names(omics)[i])]] <- sim
      
    }
    
    seu_obj <- list()
    N_sim <- length(sim_list)

    for(i in 1:N_sim){
      
    assay_name <- str_split(names(sim_list)[i], "-")[[1]][1]
    assay_name <- sub("sim_sc","",assay_name)
    seu <- CreateSeuratObject(counts = sim_list[[i]], assay = assay_name)
    seu_obj[[names(omics)[i]]] <- seu
    
    }
    
    return(seu_obj)
    
}

scRNA <- sc_omicData("scRNA-seq")
scTATC <- sc_omicData("scATAC-seq")
omic_list <- c(scRNA, scTATC)
cell_types <- list(cellA = c(1:160), cellB = c(161:270))


prova_param_est <- param_estimation(omic_list, cell_types, numberCells = c(1,2))
prova_param_est_2 <- param_estimation(omic_list, cell_types, numberCells = c(1,2), mean= c(2,3))
prova_param_est_3 <- param_estimation(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))


proma_sc_MOSim <- sc_MOSim(omic_list, cell_types, numberCells = c(10,20), mean = c(2*10^6, 2*10^3), sd = c(10^3, 10^2))
