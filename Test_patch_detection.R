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
sce = KNN_clustering(sce,K = 15,clustering_method = "Louvain",assay_type = "Count_normalised_intensity",)


plot_embedding(sce,selected_embedding = "UMAP",overlay_cluster = TRUE)

Plot_cluster_gene_expression(sce,Gene = "Myeloperoxidase.MPO",assay_type = "Count_normalised_intensity" )
Plot_cluster_gene_expression(sce,Gene = "CD9",assay_type = "Count_normalised_intensity" )
Plot_cluster_spatial(sce,Cex_parameter = 6,Image_number = 3)
Plot_gene_expression_spatial(sce,Cex_parameter = 6,assay_type = "Count_normalised_intensity",Gene = "CD9",Image_number = 3)


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
sce = KNN_clustering(sce,K = 30,assay_type = "Count_normalised_intensity",metric = "angular")


Plot_cluster_spatial(sce,Cex_parameter = 6,Image_number = 3)
Plot_gene_expression_spatial(sce,Cex_parameter = ,Image_number = 3,Gene = "CD9",
                             assay_type = "Count_normalised_intensity")

Selected_cells = sce$ImageNumber==3
Location_points = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                             Location_Y = sce$Location_Center_Y[Selected_cells])
Point_label = colLabels(sce)[Selected_cells]
List_values = levels(factor(Point_label))

x  = CE_interaction_tensor(sce, type_output = "Index", Perform_symmetrization = T)

Dixon_tensor = Dixon_interaction_tensor(sce)
CE_tensor = CE_interaction_tensor(sce)
CP_decomposition = Compute_CP_decomposition(CE_tensor,num_components = 3,Center_data = T,Scale_data = F)


Visualize_CP(CP_decomposition, k = 1)


###Patch analysis 

Patch_list = Patch_detection(sce,graph_type = "Radius",graph_parameter = 30,neighbor_degree = 1,clustering_method = "saturation",distance_composition = "Bray-Curtis",image_list = 3)
k=6
Selected_cells = sce$ImageNumber==k
Location_points = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                             Location_Y = sce$Location_Center_Y[Selected_cells])
Point_label = Patch_list[[k]]
Plot_cluster_spatial(sce,Image_number = 3,Cex_parameter = 6,Provided_cluster = Patch_list[[1]])
Plot_cluster_spatial(sce,Image_number = 3,Cex_parameter = 6)

#Getting the cluster composition of each patch

Table_patch = c()
for (k in unique(sce$ImageNumber)) {
  
  List_patch_temp = Patch_list[[k]]
  Composition_patch = table(List_patch_temp,factor(colLabels(sce))[sce$ImageNumber==k])
  Composition_patch = Composition_patch[rownames(Composition_patch)!="0",]
  rownames(Composition_patch) = paste("Image_",k,"_patch_",rownames(Composition_patch),sep = "")
  Composition_patch = as.data.frame(unclass(Composition_patch))
  Composition_patch$Image = k
  
  Table_patch = rbind(Table_patch,Composition_patch)
}
library(FactoMineR)

CA_patch = CA(Table_patch[,-ncol(Table_patch)])
plot.CA(CA_patch,axes = c(1,3))
barplot(CA_patch$eig[,2])
Patch_size = rowSums(Table_patch[,-ncol(Table_patch)])

plot(CA_patch$row$coord[,c(1,2)],pch=21,bg=string.to.colors(Table_patch$Image),cex=sqrt(Patch_size)/10)
Plot_cluster_spatial(sce,Image_number = 3,Cex_parameter = 6,Specific_cluster = 10)

####

sce_Jana = readRDS("Desktop/Data_Jana/sce.rds")
temp_cluster = as.numeric(factor(sce_Jana$metacluster))
temp_cluster[temp_cluster%in%c(14:27)] = 14
colLabels(sce_Jana) = temp_cluster
sce_Jana$ImageNumber = sce_Jana$ImageNb
metadata(sce_Jana) = list(N_core = 6,dimension = "2D")

Patch_list = Patch_detection(sce_Jana,graph_type = "KNN",graph_parameter = 40,neighbor_degree = 1,
                             clustering_method = "saturation",distance_composition = "Bray-Curtis",image_list = 1:12)


k=12
Selected_cells = sce$ImageNumber==k
Location_points = data.frame(Location_X = sce_Jana$Location_Center_X[Selected_cells],
                             Location_Y = sce_Jana$Location_Center_Y[Selected_cells])
Point_label = Patch_list[[k]]
Plot_cluster_spatial(sce_Jana,Image_number = k,Cex_parameter = 3)
Plot_cluster_spatial(sce_Jana,Image_number = k,Cex_parameter = 3,Provided_cluster = Point_label)

#Getting the cluster composition of each patch


### Applying PenGe

library(PenGE)

Selected_cells = sce$ImageNumber==3
Location_points = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                             Location_Y = sce$Location_Center_Y[Selected_cells])
Point_label = colLabels(sce)[Selected_cells]
List_values = levels(factor(Point_label))

PPP_temp = ppp(x = Location_points$Location_X,y=Location_points$Location_Y,
               window = owin(range(Location_points$Location_X),range(Location_points$Location_Y)),
               marks = factor(Point_label))
p = length(unique(Point_label))
rv <- c(50,100, 200) # use these for everything
r1 <- rep( list(rv) , p) # intra-type interaction ranges
r2 <- rep( list(rv), p * (p-1)/2) # inter-type
Q <- make_Q_stepper_multi(x=PPP_temp, ranges1 = r1, ranges2 = r2)
fit <- fitGlbin_CV(Q,verb = 1,nlambda = 30,log_linear_lambda = F,mc.cores = 6)
M <- imatrix(fit, signed = TRUE)
image.imatrix(M, col = c("red","black","green","blue"), cex.axis = 1)
Rest <- residuals_cv(fit, x = PPP_temp, mc.cores = 6)






x = potential(fit,i=6,j=7)

eval_estimated_potential(at=200,fit = fit,k = 20,j = 6,i = 6)


