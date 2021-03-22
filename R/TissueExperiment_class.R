library(spatstat)
library(FactoMineR)
library(foreach)
library(doParallel)

setClass("TissueExperiment", representation(r_range = "numeric",
                                            Contingency_table = "table",
                                            CA_analysis = "list",
                                            List_global_K_function = "data.frame",
                                            List_global_pcf_function = "data.frame",
                                            List_cluster_K_function = "list",
                                            List_cluster_pcf_function = "list"))


Create_TissueExperiment = function(sce,r_range = 0:250) {
  
  N_Images = length(unique(sce$ImageNumber))
  #Computing zero/first order statistics and CA object
  
  cat("Computing Correspondance Analysis object...")
  Contingency_table = table(sce$ImageNumber,colLabels(sce))
  colnames(Contingency_table) = paste("Cluster_",colnames(Contingency_table),sep = "")
  rownames(Contingency_table) = paste("Image_",rownames(Contingency_table),sep = "")
  
  CA_analysis = CA(Contingency_table,graph = F)
  
  cat("done ! \n")
  
  
  #Computing global Ripley K functions 
  cat("Computing global Ripley's K functions..")
  registerDoParallel(cores = metadata(sce)$N_core)
  List_Kest = foreach(i=1:N_Images) %dopar% {
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==i])

    
    Temp_ppp_object = ppp(x = Temp_location_data$X,y=Temp_location_data$Y,
                          window = owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Kest_temp = Kest(Temp_ppp_object,correction = "isotropic",r = r_range)
    Kest_temp
  }
  
  Merged_Kest = List_Kest[[1]]
  
  for (k in 2:length(List_Kest)) {
    Merged_Kest = suppressWarnings(bind.fv(Merged_Kest,List_Kest[[k]]$iso))
    
  }
  names(Merged_Kest) = c("r","theo",paste("Kest_Image_",1:length(List_Kest),sep = ""))
  fvlabels(Merged_Kest) = c("r","Kest[pois](r)",paste("Kest[Image_",1:length(List_Kest),"](r)",sep = ""))
  cat(" done ! \n")
  
  #Computing global pcf function 
  cat("Computing global pair Correlation Functions (pcf)..")
  
  List_pcf = foreach(i=unique(sce$ImageNumber)) %dopar% {
    #Estimating the pcf from the K function using the following strategy (Spatstat doc)
    # apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, estimate the derivative of Y, and solve
    
    pcf_temp = pcf(List_Kest[[i]],method = "b")
    
  }
  
  Merged_pcf = List_pcf[[1]]
  
  for (k in 2:length(List_pcf)) {
    Merged_pcf = suppressWarnings(bind.fv(Merged_pcf,List_pcf[[k]]$pcf))
    
  }
  names(Merged_pcf) = c("r","theo",paste("Pcf_Image_",1:length(List_pcf),sep = ""))
  fvlabels(Merged_pcf) = c("r","pcf[pois](r)",paste("pcf[Image_",1:length(List_pcf),"](r)",sep = ""))
  cat(" done ! \n")
  
  #Computing the cluster specific K functions
  
  cat("Computing the cluster-specific K and pair correlation functions in each image... \n")
  List_cluster_K_functions = c()
  List_cluster_pcf = c()
  
  for (k in 1:N_Images) {
    
    #Computing the K function for each cluster in a parallelised way..
    
    cat(paste("Processing of Image"),as.character(k),"...")
    List_Cluster_K_function_temp = foreach(i=unique(colLabels(sce))) %dopar% {
      
      Selected_cells = sce$ImageNumber==k & colLabels(sce)==i

      if (sum(Selected_cells)>0) {
        Temp_ppp_object = ppp(x = sce$Location_Center_X[Selected_cells],
                              y = sce$Location_Center_Y[Selected_cells],
                              window = owin(xrange = range(sce$Location_Center_X[Selected_cells]),
                                            yrange = range(sce$Location_Center_Y[Selected_cells])))
        
        Kest_temp = Kest(Temp_ppp_object,correction = "isotropic",r = r_range)
      }
      
      if (sum(Selected_cells)==0) {
        Kest_temp = fv(data.frame(r = r_range,iso = NA),argu = "r",valu = "iso")
      }
      Kest_temp
    }
    cat(" done !\n")
    
    
    #Aggregating the K functions from the same image 
    
    Cluster_K_function_temp_merged = List_Cluster_K_function_temp[[1]]
    
    for (i in 2:length(List_Cluster_K_function_temp)) {
      Cluster_K_function_temp_merged = suppressWarnings(bind.fv(Cluster_K_function_temp_merged,List_Cluster_K_function_temp[[i]]$iso))
      
    }
    names(Cluster_K_function_temp_merged) = c("r","theo",paste("Kest_cluster",1:length(List_Cluster_K_function_temp),sep = ""))
    fvlabels(Cluster_K_function_temp_merged) = c("r","pcf[pois](r)",paste("pcf[Cluster_",1:length(List_Cluster_K_function_temp),"](r)",sep = ""))
    List_cluster_K_functions[[k]] = Cluster_K_function_temp_merged
    
    
    #Generating and aggregating the pcf from the same image 
    
    if (sum(is.na(List_Cluster_K_function_temp[[1]]$iso))>0) {
      Cluster_pcf_temp_merged = fv(data.frame(r=r_range,theo=1,pcf_1 = NA),argu = "r",valu = "theo")
    }
    
    if (sum(is.na(List_Cluster_K_function_temp[[1]]$iso))==0) {
      Cluster_pcf_temp_merged = pcf.fv(List_Cluster_K_function_temp[[1]],method = "b")
    }
    

    for (i in 2:length(List_Cluster_K_function_temp)) {
      
      if (sum(is.na(List_Cluster_K_function_temp[[i]]$iso))>0) {
        x = NA
      }
      
      if (sum(is.na(List_Cluster_K_function_temp[[i]]$iso))==0) {
        x = suppressWarnings(pcf.fv(List_Cluster_K_function_temp[[i]],method = "b")$pcf)
      }
      
      Cluster_pcf_temp_merged = bind.fv(Cluster_pcf_temp_merged,x)
      
    }
    names(Cluster_pcf_temp_merged) = c("r","theo",paste("pcf_cluster",1:length(List_Cluster_K_function_temp),sep = ""))
    fvlabels(Cluster_pcf_temp_merged) = c("r","pcf[pois](r)",paste("pcf[Cluster_",1:length(List_Cluster_K_function_temp),"](r)",sep = ""))
    List_cluster_pcf[[k]] = Cluster_pcf_temp_merged
  }
  
  #Computing the spatial simpson index alpha 
  
  Spatial_alpha_index = foreach(k=1:N_Images,.combine = cbind) %dopar% {
    
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==k])
    
    #Computing the lambda parameters
    Area_size = spatstat::area(owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Lambda_parameter_global = nrow(Temp_location_data)/Area_size
    Lambda_parameter_cluster = table(factor(colLabels(sce))[sce$ImageNumber==k])/Area_size
    
    #Extracting the K and pair correlation function
    Temp_global_K_function = as.data.frame(Merged_Kest)[,k+2]
    
    Temp_cluster_K_function = as.data.frame(List_cluster_K_functions[[k]])
    Temp_cluster_K_function = Temp_cluster_K_function[,c(-1,-2)]
    
    #na_columns = colSums(is.na(Temp_cluster_K_function))==nrow(Temp_cluster_K_function)
    #Temp_cluster_K_function[,na_columns] = 0
        
    Temp_cluster_K_function = apply(X = Temp_cluster_K_function,MARGIN = 1,FUN = function(x) {x*Lambda_parameter_cluster^2})
    Temp_cluster_K_function = t(Temp_cluster_K_function)
    Temp_cluster_K_function = apply(X = Temp_cluster_K_function, MARGIN = 2,FUN = function(x) {x/(Lambda_parameter_global^2*Temp_global_K_function)})
    
    Alpha_index = 1- rowSums(Temp_cluster_K_function,na.rm = T)
    Alpha_index      
  }
  
  
  #Computing the spatial simpson index beta 
  
  Spatial_beta_index = foreach(k=1:N_Images,.combine = cbind) %dopar% {
    
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==k])
    
    #Computing the lambda parameters
    Area_size = spatstat::area(owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Lambda_parameter_global = nrow(Temp_location_data)/Area_size
    Lambda_parameter_cluster = table(factor(colLabels(sce))[sce$ImageNumber==k])/Area_size
    
    #Extracting the K and pair correlation function
    Temp_global_pcf = as.data.frame(Merged_pcf)[,k+2]
    
    Temp_cluster_pcf = as.data.frame(List_cluster_pcf[[k]])
    Temp_cluster_pcf = Temp_cluster_pcf[,c(-1,-2)]
    

    Temp_cluster_pcf = apply(X = Temp_cluster_pcf,MARGIN = 1,FUN = function(x) {x*Lambda_parameter_cluster^2})
    Temp_cluster_pcf = t(Temp_cluster_pcf)
    Temp_cluster_pcf = apply(X = Temp_cluster_pcf, MARGIN = 2,FUN = function(x) {x/(Lambda_parameter_global^2*Temp_global_pcf)})
    
    Beta_index = 1- rowSums(Temp_cluster_pcf,na.rm = T)
    Beta_index      
  }
  
  
  
}





