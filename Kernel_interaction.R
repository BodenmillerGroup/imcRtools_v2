
Sigma_parameter = 10
Distribution_metric = "Bhattacharyya"

N_cluster = max(colLabels(sce))
N_image = max(sce$ImageNumber)
List_score_matrix = c()



List_score_matrix=foreach(k = 1:N_image) %dopar% {
  
  cat(paste("Computing interactions for image",as.character(k),"\n"))
  Density_score_matrix = matrix(data = 0,nrow = max(colLabels(sce)),ncol= N_cluster)
  X_range = range(sce$Location_Center_X[sce$ImageNumber==k])
  Y_range = range(sce$Location_Center_Y[sce$ImageNumber==k])
  
  
  for (i in 1:N_cluster) {
    
    for (j in 1:N_cluster) {
      
      X_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k & colLabels(sce)==i ],
                             Y = sce$Location_Center_Y[sce$ImageNumber==k & colLabels(sce)==i ])
      X_pattern = ppp(X_pattern$X,X_pattern$Y,window = owin(X_range,yrange = Y_range))
      
      Y_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k & colLabels(sce)==j ],
                             Y = sce$Location_Center_Y[sce$ImageNumber==k & colLabels(sce)==j ])
      Y_pattern = ppp(Y_pattern$X,Y_pattern$Y,window = owin(X_range,yrange = Y_range))
      
      if (i!=j & npoints(X_pattern)>10 & npoints(Y_pattern)>10) {
        
        Density_X = density.ppp(X_pattern,sigma = Sigma_parameter)
        Density_Y = density.ppp(Y_pattern,sigma = Sigma_parameter)
        
        if (Distribution_metric=="Bhattacharyya") {
          Distribution_distance = sqrt(Density_X$v)*sqrt(Density_Y$v)
          Distribution_distance[is.na(Distribution_distance)] = 0
          Distribution_distance = -log(sum(Distribution_distance))
        }
        
        if (Distribution_metric=="KL") {
          Distribution_distance = Density_X$v * log(Density_X$v/Density_Y$v)
          Distribution_distance[is.na(Distribution_distance)] =0
          Distribution_distance[is.infinite(Distribution_distance)] =0
          Distribution_distance = sum(Distribution_distance)
        }
        
        
      }
      
      if (i==j & npoints(X_pattern)>10) {
        Distribution_distance = 0
      }
      Density_score_matrix[i,j] = Distribution_distance
      
    }
  }
  Density_score_matrix
}


for (k in 1:N_image) {
  Tensor_object[,,k] = List_score_matrix[[k]]
}

CP_KL = Compute_CP_decomposition(Interaction_tensor = Tensor_object,max_iter = 1000,tol = 1e-8,num_components = 10)
pheatmap(x$List_typical_matrix[[3]],cluster_cols = F,cluster_rows = F)


P_value_table = matrix(data = 1,ncol = ncol(CP_KL$Score_matrix),nrow=ncol(Annotation_table))

for (i in 1:ncol(P_value_table)) {
  for (j in 1:nrow(P_value_table)) {
    test_temp = kruskal.test(CP_KL$Score_matrix[,i]~Annotation_table[,j])
    P_value_table[j,i] = test_temp$p.value
  }
}

rownames(P_value_table) = colnames(Annotation_table)
colnames(P_value_table) = 1:ncol(P_value_table)
P_value_table = -log10(P_value_table)
pheatmap(P_value_table,clustering_method = "ward")
boxplot(CP_KL$Score_matrix[,7]~Annotation_table$Relapse,outline=F)
pheatmap(CP_KL$List_typical_matrix[[7]],cluster_rows = F,cluster_cols = F,clustering_method = "ward")


