library(SingleCellExperiment)
library(doParallel)
library(foreach)
library(uwot)
library(N2R)
library(geometry)
library(dbscan)
library(FactoMineR)
library(igraph)
library(RColorBrewer)
library(readr)


##Two basic functions for convertion of numerical and string vectors to color vectors

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

cluster_to_color = function(cluster_vector,Defined_list_cluster = NULL) {
  
  cluster_vector = as.character(cluster_vector)
  List_unique_cluster = unique(cluster_vector)
  List_unique_cluster = List_unique_cluster[order(List_unique_cluster)]
  
  if (!is.null(Defined_list_cluster)) {
    List_unique_cluster = List_unique_cluster
  }
  
  N_clusters = length(unique(cluster_vector))
  
  optimal_palette = suppressWarnings(colorRamp(brewer.pal(N_clusters, "Spectral")))
  optimal_palette = optimal_palette((1:N_clusters)/N_clusters)
  optimal_palette = optimal_palette / 255
  optimal_palette = rgb(optimal_palette)
  
  color_cluster = cluster_vector
  
  for (k in 1:N_clusters) {
    selected_cluster = List_unique_cluster[k]
    color_cluster[cluster_vector==k] = optimal_palette[k]
  }
  return(color_cluster)
}


#Loading the data 

Load_IMC_data = function(Path_to_cell_file,Path_to_panel_file,
                         Name_target = "Target",Name_metal="Metal.Tag",
                         Regex_expression_channel = "Intensity_MeanIntensity_FullStackFiltered_c",
                         DNA_channel = NULL,Use_for_clustering = NULL) {
  
  Panel_table = read.csv(Path_to_panel_file)
  Panel_table = Panel_table[Panel_table$full==1,]
  
  
  cat("Loading the cell data ... ")
  Raw_data = readr::read_csv(Path_to_cell_file,progress = T)
  Raw_data = as.data.frame(Raw_data)
  cat("done ! \n ")
  
  #Extracting the mean intensity values
  
  Expression_data = Raw_data[,base::grepl(colnames(Raw_data),pattern = Regex_expression_channel)]
  l = colnames(Expression_data)
  l = strsplit(l,split = Regex_expression_channel)
  l = unlist(lapply(l,FUN = function(x) {x[2]}))
  l = as.numeric(l)

  colnames(Expression_data) = make.names(Panel_table[l,Name_target],unique = T)
  
  
  #Extracting the cell annotation table
  
  Cell_annotation = data.frame(ImageNumber = Raw_data$ImageNumber,
                               ObjectNumber = Raw_data$ObjectNumber,
                               Location_Center_X = Raw_data$Location_Center_X,
                               Location_Center_Y = Raw_data$Location_Center_Y,
                               Cell_size = Raw_data$AreaShape_Area)
  
  #Extracting the gene annotation table
  
  Gene_annotation = data.frame(Used_for_clustering = colnames(Expression_data)%in%Use_for_clustering,
                                 DNA_channel = colnames(Expression_data)%in%DNA_channel,
                                 row.names = colnames(Expression_data))
  return(list(Expression_data = Expression_data,
              Cell_annotation = Cell_annotation,
              Gene_annotation = Gene_annotation))


}


##Creating the object

Create_SCE = function(List_data,dimension = "2D",Bit_mode=16,N_core = 6) {
  
  sce = SingleCellExperiment(assays = list(Raw_intensity = as.matrix(t(List_data$Expression_data))),
                             colData = List_data$Cell_annotation,
                             rowData = List_data$Gene_annotation,
                             metadata = list(dimension = dimension,Bit_mode=Bit_mode,N_core = N_core))
  return(sce)
}


#Normalisation of the data

Poisson_normalization = function(sce,Bit_mode = 32,perform_batch_correction=FALSE,
                                 batch_vector=NULL,residual_normalisation = "Anscombe") {
  
  #Transforming the data back to count data
  Cell_size = sce$Cell_size
  Transformed_data = t(sce@assays@data[[1]])
  Transformed_data = Transformed_data * 2^metadata(sce)$Bit_mode
  Transformed_data = apply(Transformed_data,MARGIN = 2,FUN = function(x) {x*Cell_size})
  
  Transformed_data = round(Transformed_data)
  
  #Creating parallel backend
  cat(paste("Creating parallel backend using"),as.character(metadata(sce)$N_core),"cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  #Performing the poisson regression
  if (!perform_batch_correction) {
    cat("Fitting Poisson regressions ...")
    
    List_regression_model =foreach(i=colnames(Transformed_data)) %dopar% {
      Poisson_model = glm(Transformed_data[,i]~log(Cell_size),family = "poisson")
    }
  }
  cat(" done ! \n")
  
  #Extracting and normalizing residuals
  
  
  Residual_matrix =foreach(i=1:length(List_regression_model),.combine = cbind ) %dopar% {
    
    Fitted_values = List_regression_model[[i]]$fitted.values
    Real_values = Transformed_data[,i]
    
    if (!residual_normalisation%in%c("Anscombe","Pearson","Working","VST")){
      cat("No proper method for residual normalization provided. Using the Anscombe normalization method")
    }
    
    if (residual_normalisation=="Anscombe") {
      Normalised_residuals = 1.5 * (Real_values^(2/3)-Fitted_values^(2/3)) / (Fitted_values^(1/6))
    }
    
    if (residual_normalisation=="Pearson") {
      Normalised_residuals = (Real_values-Fitted_values)/sqrt(Fitted_values)
    }
    
    if (residual_normalisation=="Working") {
      Normalised_residuals = (Real_values-Fitted_values)/Fitted_values
    }
    
    if (residual_normalisation=="VST") {
      Normalised_residuals = (sqrt(Real_values)-sqrt(Fitted_values))/2
    }
    
    Normalised_residuals
    
    
  }
  
  #Scaling to 0-1 values
  Residual_matrix = apply(Residual_matrix,MARGIN = 2,FUN = function(x) {
    x = x-min(x)
    x = x/max(x)
  })
  
  Residual_matrix = as.data.frame(Residual_matrix)
  colnames(Residual_matrix) = colnames(Transformed_data)
  
  #Adding the matrix a new assay slot 
  assay(sce, "Normalized intensity",withDimnames = F) <- t(Residual_matrix)
  
  return(sce)
}


##UMAP embedding 

UMAP_embedding = function(sce,metric="correlation",n_neighbors = 30) {
  
  if (!"Normalized intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  data_to_project =assays(sce)[["Normalized intensity"]]
  data_to_project = t(data_to_project)
  
  umap_embedding = umap(data_to_project,n_neighbors = n_neighbors,pca = NULL,metric = metric,verbose=T)
  
  reducedDim(sce, "UMAP") <- umap_embedding
  return(sce)
}

###Graph based clustering


KNN_clustering = function(sce,K=30,clustering_method = "Louvain") {
  if (!"Normalized intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  
  data_to_cluster =assays(sce)[["Normalized intensity"]]
  data_to_cluster = t(data_to_cluster)
  
  Channel_for_clustering = rowData(sce)$Used_for_clustering
  
  #If no selection of the channels : using all channels
  if (is.null(Channel_for_clustering)) {
    Channel_for_clustering = rep(TRUE,ncol(data_to_cluster))
  }
  
  data_to_cluster = data_to_cluster[,Channel_for_clustering]
  
  cat(paste("Clustering of the data using",as.character(sum(Channel_for_clustering)),"channels \n"))
  
  KNN_graph_matrix =  N2R::Knn(as.matrix(data_to_cluster), K, nThreads=metadata(sce)$N_core, verbose=T, indexType='angular')
  KNN_graph_matrix_symmetric = (KNN_graph_matrix+t(KNN_graph_matrix))/2
  Final_graph <- igraph::graph_from_adjacency_matrix(KNN_graph_matrix_symmetric,mode='undirected',weighted=TRUE)
  
  if (clustering_method == "Louvain") {
    Clustering <- igraph::multilevel.community(Final_graph)
  }
  
  if (clustering_method == "Greedy") {
    Clustering <- igraph::cluster_fast_greedy(Final_graph,modularity = T)
  }
  
  if (clustering_method == "infomap") {
    Clustering <- igraph::infomap.community(Final_graph,modularity = T)
  }
  
  Clustering_group = membership(Clustering)
  colLabels(sce) = Clustering_group
  cat(paste(as.character(length(unique(Clustering_group))),"clusters have been identified \n"))
  return(sce)
}


##Plotting UMAP and if possible overlay clustering on it 

plot_UMAP = function(sce,overlay_cluster) {
  if (!"UMAP"%in%reducedDimNames(sce)) {
    stop("UMAP embedding has not been computed. Please compute it first \n")
  }
  
  Embedding = reducedDim(sce)
  if (class(Embedding)=="list") {
    Embedding = Embedding[["UMAP"]]
  }
  par(las=1,bty="l")
  
  if (is.null(colLabels(sce))) {
    plot(Embedding,pch=21,bg="red3", xlab="UMAP 1",ylab="UMAP 2",cex=0.7)
  }
  
  if (!is.null(colLabels(sce))) {
    plot(Embedding,pch=21,bg=cluster_to_color(colLabels(sce)), xlab="UMAP 1",ylab="UMAP 2",cex=0.7)
  }
}

##Plotting the expression of a given gene in a giveb image

Plot_gene_expression_spatial = function(sce,Image_number = 1,Gene=NULL) {
  if (is.null(Gene) | !Gene%in%rownames(sce@assays@data@listData$Raw_intensity) ) {
    stop("Please select a correct gene to plot ! \n")
  }
  
  if (!"Normalized intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  
  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_expression_data = as.numeric(sce@assays@data@listData$`Normalized intensity`[Gene,])
  Temp_size_data = sqrt(sce$Cell_size[sce$ImageNumber==Image_number])
  par(las=1,bty="l")
  plot(Temp_location_data,cex=Temp_size_data/10,pch=21,bg=color_convertion(Temp_expression_data))
  
}

##Plotting the expression of a given gene in a special image

Plot_cluster_spatial = function(sce,Image_number = 1) {
  if (is.null(colLabels(sce))) {
    stop("Please compute clustering first ! \n")
  }
  
  
  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_cluster_data = colLabels(sce)[sce$ImageNumber==Image_number]
  
  Dimension = metadata(sce)$dimension
  
  
  if (Dimension == "2D") {
    Temp_size_data = sqrt(sce$Cell_size[sce$ImageNumber==Image_number])
  }
  
  if (Dimension == "3D") {
    Temp_size_data = (sce$Cell_size[sce$ImageNumber==Image_number])^1/3
  }
  
  par(las=1,bty="l")
  plot(Temp_location_data,cex=Temp_size_data/10,pch=21,bg=cluster_to_color(Temp_cluster_data))
  
}


##Compute the neighboorhood graph  : graph_parameter -> either the number of neighbors for kNN, the distance threshold for fixed radius NN
Compute_neighborhood_graph = function(sce,graph_type,graph_parameter = 5) {
  
  List_graph = list()
  
  if (!graph_type%in%c("KNN","FRNN","Delaunay")) {
    stop("The type of graph required is not available. Please choose among KNN,Epsilon or Delaunay  ! \n")
  }
  
  
  if (graph_type == "KNN") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
    
      N = nrow(Temp_location_data)
      

      KNN_graph_matrix =  N2R::Knn(as.matrix(Temp_location_data),k = graph_parameter, nThreads=1, indexType='L2')
      Temp_graph = igraph::graph_from_adjacency_matrix(KNN_graph_matrix,mode = "undirected")
      Temp_graph = igraph::simplify(Temp_graph,remove.multiple = T,remove.loops = T)
      
      List_graph[[i]] = Temp_graph
    }
  }
  
  
  if (graph_type == "FRNN") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
      
      FRNN_list_edges = frNN(x = Temp_location_data,eps = graph_parameter)
      FRNN_list_edges = FRNN_list_edges$id
      FRNN_graph = igraph::graph_from_adj_list(FRNN_list_edges,mode = "all",duplicate = T)
      List_graph[[i]] = FRNN_graph
    }
  }
  
  if (graph_type == "Delaunay") {
    
    for (i in unique(sce$ImageNumber)) {
      Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==i],
                                      Y = sce$Location_Center_Y[sce$ImageNumber==i])
      

      Delaunay_triangulation = delaunayn(Temp_location_data,output.options = "Fn")
      Delaunay_triangulation_graph = igraph::graph_from_adj_list(Delaunay_triangulation$neighbours,mode = "all",duplicate = T)
      Delaunay_triangulation_graph =  igraph::simplify(Delaunay_triangulation_graph)
      List_graph[[i]] = Delaunay_triangulation_graph
    }
  }
  
 sce@metadata$List_graphs = List_graph
 return(sce)
   
}

##Loading an 'Object relationship' object and create the corresponding graph

Load_object_relationship = function(sce , Path_to_object_relationship) {
  
  Object_relationship = readr::read_csv(Path_to_object_relationship,progress = T)
  Object_relationship = as.data.frame(Object_relationship)
  
  Object_relationship_reshape = data.frame(ImageNumber = Object_relationship$`First Image Number`,
                                           Cell_number_1 = Object_relationship$`First Object Number`,
                                           Cell_number_2 = Object_relationship$`Second Object Number`)
  
  List_graph = list()
  
  for (i in unique(sce$ImageNumber)) {
    
    Object_relationship_reshape_temp = as.matrix(Object_relationship_reshape[Object_relationship_reshape$ImageNumber==i,-1])
    Temp_graph = igraph::graph_from_edgelist(el = Object_relationship_reshape_temp,directed = F)
    List_graph[[i]] = Temp_graph
  }
  
  sce@metadata$List_graphs = List_graph
  return(sce)
  
}

##Study of neibourhood : creation of a new object 


setClass("NeighborhoodObject", representation(Image_number = "character",Graph = "numeric"))

Create_NeighborhoodObject = function(sce) {
  
  if (is.null(sce@metadata$List_graphs)) {
    stop("Compute or load a neighborhood graph first ! \n")
  }
  
  
  
}







