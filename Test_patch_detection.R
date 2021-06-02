library(balagan)
library(pheatmap)
library(igraph)
library(spatstat)


color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)

}
List_data = Load_IMC_data(Path_to_cell_file = "Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/cell.csv",
                          Path_to_panel_file ="Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/panel.csv",
                          Regex_expression_channel = "Intensity_MeanIntensity_FullStackFiltered_c")


sce = Create_SCE(List_data,dimension = "2D",Bit_mode = 16,N_core = 6)
sce = Arcsinh_normalization(sce)
sce = Count_normalization(sce)

sce = UMAP_embedding(sce,assay_type = "Count_normalised_intensity")
sce = KNN_clustering(sce,K = 30,clustering_method = "Louvain",assay_type = "Count_normalised_intensity")


plot_embedding(sce,selected_embedding = "UMAP",overlay_cluster = TRUE)

Plot_cluster_gene_expression(sce,Gene = "Myeloperoxidase.MPO",assay_type = "Count_normalised_intensity" )
Plot_cluster_gene_expression(sce,Gene = "CD45RA",assay_type = "Raw_intensity" )
Plot_cluster_spatial(sce,Cex_parameter = 6,Image_number = 6)
Plot_gene_expression_spatial(sce,Cex_parameter = 6,assay_type = "Raw_intensity",Gene = "CD45RA",Image_number = )


pheatmap(cor(t(assay(sce,i = "Raw_intensity"))),clustering_method = "ward.D2")
pheatmap(cor(t(assay(sce,i = "Count_normalised_intensity"))),clustering_method = "ward.D2")
pheatmap(cor(t(assay(sce,i = "Arcsinh_transformed_intensity"))),clustering_method = "ward.D2")


x = aggregate(t(assay(sce,"Raw_intensity")),by=list(as.character(colLabels(sce))),FUN = mean)
rownames(x) =x$Group.1
x =x[,-1]
pheatmap(t(x),scale="row",clustering_method = "ward")

### Using the nuc-cyt function


List_data = Load_IMC_data_nuc_cyt(Path_to_cell_file = "~/Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/cell.csv",
                                              Path_to_nuc_file = "~/Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/nuc.csv",
                                              Path_to_cyt_file =  "~/Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/cyt.csv",
                          Path_to_panel_file ="~/Desktop/Lymph_node_panel_building/IMC_data/20210506_Test_panel_LN_1/analysis/cpout/panel.csv",
                          Regex_expression_channel = "Intensity_MeanIntensity_FullStackFiltered_c")


sce = Create_SCE(List_data,dimension = "2D",Bit_mode = 16,N_core = 6)
rowData(sce)$Used_for_clustering = !rownames(sce)%in%c("CD40L..CD154.","CD274..B7.H1..PD.L1." ,"AICDA","MX1","Ki.67","CD278..ICOS","Histone.H3K9Ac")

sce = Arcsinh_normalization(sce)
sce = Count_normalization(sce,residual_normalisation = "Anscombe")
sce = UMAP_embedding(sce,assay_type = "Count_normalised_intensity",metric = "cosine",n_neighbors = 15)
sce = KNN_clustering(sce,K = 30,assay_type = "Count_normalised_intensity",metric = "angular")

plot_embedding(sce) 


Plot_cluster_spatial(sce,Cex_parameter = 6,Image_number = 3)
Plot_gene_expression_spatial(sce,Cex_parameter = ,Image_number = 2,Gene = "CD31",
                             assay_type = "Count_normalised_intensity")

Selected_cells = sce$ImageNumber==3
Location_points = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                             Location_Y = sce$Location_Center_Y[Selected_cells])
Point_label = colLabels(sce)[Selected_cells]
List_values = levels(factor(Point_label))


#Define the graph : kNN graph
K_parameter = 15
KNN_graph_matrix = N2R::Knn(as.matrix(Location_points), K_parameter, nThreads = metadata(sce)$N_core, 
                       verbose = F, indexType = "L2")
KNN_graph_matrix = KNN_graph_matrix + t(KNN_graph_matrix)
Final_graph <- graph_from_adjacency_matrix(KNN_graph_matrix, 
                                           mode = "undirected", weighted = TRUE)

#Define the graph : Gabriel graph

Gabriel_graph_adj_list = spatgraph(ppp(x = Location_points$Location_X,y = Location_points$Location_Y,
                  window = owin(range(Location_points$Location_X),yrange = range(Location_points$Location_Y))),type = "gabriel")
Final_graph = graph_from_adj_list(Gabriel_graph_adj_list$edges,"all")



List_neighbor = neighborhood(graph = Final_graph)
List_neighbor_labels = lapply(List_neighbor,FUN = function(x) {Point_label[x]})
List_neighbor_labels = lapply(List_neighbor_labels,FUN = function(x) {table(factor(x,levels = List_values))/length(x)})

Table_neighbor_labels = matrix(data = unlist(List_neighbor_labels),ncol = length(List_values),byrow = T)
Entropy_labels = rowSums(-Table_neighbor_labels * log(Table_neighbor_labels),na.rm = T)

Final_graph = simplify(Final_graph, )
List_interactions = as_edgelist(Final_graph,)
List_weights = c()

#Only comparing the non-shared neighboorhood

List_weights = c()
for (k in 1:nrow(List_interactions)) {
  
  node_i = List_interactions[k,1]
  node_j = List_interactions[k,2]
  
  Neighbor_i = as.numeric(neighborhood(graph = Final_graph,nodes = node_i,order = 2)[[1]])
  Neighbor_j = as.numeric(neighborhood(graph = Final_graph,nodes = node_j,order = 2)[[1]])
  
  Neighbor_i_specific = Neighbor_i[!Neighbor_i%in%Neighbor_j]
  Neighbor_j_specific = Neighbor_j[!Neighbor_j%in%Neighbor_i]
  
  composition_i = Point_label[Neighbor_i_specific]
  composition_j = Point_label[Neighbor_j_specific]
  
  composition_i_normalized = table(factor(composition_i,levels = List_values))/length(List_values)
  composition_j_normalized = table(factor(composition_j,levels = List_values))/length(List_values)

  List_weights = c(List_weights,sum((composition_i_normalized*composition_j_normalized))/(sqrt(sum(composition_i_normalized^2))*sqrt(sum(composition_j_normalized^2))))
}
List_weights[is.na(List_weights)] = 0

#Transforming the weights 
List_weights_transformed = 1-exp(-List_weights^2/(2*quantile(List_weights,0.66)))
E(Final_graph)$weight = List_weights_transformed

Distance_matrix_temp = as.dist(1-as_adjacency_matrix(Final_graph,attr = "weight"))
x = dbscan::kNNdistplot(x = Distance_matrix_temp,k = 1)
  

library(dbscan)
x = dbscan(Distance_matrix_temp,eps = 0.02,minPts = 5)

Clustering_louvain= cluster_fast_greedy(Final_graph,weights = List_weights_transformed,modularity = T,membership = T)

Plot_cluster_spatial(sce,Image_number = 3,Cex_parameter = 6)

plot.igraph(Final_graph,layout = as.matrix(Location_points),vertex.label=NA,vertex.size=1.5,
            edge.width = E(Final_graph)$weight,vertex.color=string.to.colors(x$cluster==7))


Test_patch_LN_degree_1 = Patch_detection(sce = sce,graph_type = "Radius",graph_parameter=30,neighbor_degree = 1,clustering_method = "greedy")
Test_patch_LN_degree_2 = Patch_detection(sce = sce,graph_type = "Radius",graph_parameter=20,neighbor_degree = 2,clustering_method = "greedy")

k = 3
Plot_cluster_spatial(sce,Image_number = k,Cex_parameter = 6)
Selected_image = sce$ImageNumber==k
plot(sce$Location_Center_X[Selected_image],sce$Location_Center_Y[Selected_image],pch=21,bg=string.to.colors(Test_patch_LN_degree_2[[k]]))


