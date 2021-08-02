Pcf_interaction_tensor = function(sce,r_range=seq(0,300,length.out=50)) {
  

  N_cluster = max(colLabels(sce))
  N_image = max(sce$ImageNumber)

  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_pcf_tensors=foreach(k = 1:N_image) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    X_range = range(sce$Location_Center_X[sce$ImageNumber==k])
    Y_range = range(sce$Location_Center_Y[sce$ImageNumber==k])
    
    Tensor_object_image_temp = rand_tensor(modes = c(N_cluster, N_cluster, length(r_range)), drop = FALSE)
    Tensor_object_image_temp[,,]=0
    
    for (i in 1:N_cluster) {
      
      x = rep(1,length(r_range))
      for (j in 1:N_cluster) {
        
        Temp_point_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k & colLabels(sce)%in%c(i,j) ],
                                        Y = sce$Location_Center_Y[sce$ImageNumber==k & colLabels(sce)%in%c(i,j) ])
        Temp_point_pattern = ppp(Temp_point_pattern$X,Temp_point_pattern$Y,
                                 window = owin(X_range,yrange = Y_range),
                                 marks = factor(colLabels(sce)[sce$ImageNumber==k & colLabels(sce)%in%c(i,j)]))

        if (i!=j & sum(marks(Temp_point_pattern)==i)>15 & sum(marks(Temp_point_pattern)==j)>15 ) {
          Temp_Kest = Kcross(Temp_point_pattern,i = i,j=j,r = r_range,correction = "iso")
          Normalized_Lest = sqrt(Temp_Kest$iso/pi) - Temp_Kest$r
        }
        
        if (i==j & npoints(Temp_point_pattern)>10) {
          Temp_Kest = Kest(Temp_point_pattern ,r= r_range,correction = "iso")

          Normalized_Lest = sqrt(Temp_Kest$iso/pi) - Temp_Kest$r
          
        }
        
        Normalized_Lest[is.na(Normalized_Lest)] = 0
        Tensor_object_image_temp[i,j,] = Normalized_Lest
        
      }
    }
    #Cleaning and symmetrizing the CE index matrix
    Tensor_object_image_temp
    
    }
  
  cat("Storing the data into a three-order tensor \n")
  
  Tensor_object = rand_tensor(modes = c(N_cluster, N_cluster,length(r_range) ,N_image), drop = FALSE)
  for (k in 1:N_image) {
    Tensor_object[,,,k] = List_pcf_tensors[[k]]
  }
  

  x = cp(Tensor_object,num_components = 10)
  
  
  P_value_table_2 = matrix(data = 1,ncol = ncol(x$U[[1]]),nrow=ncol(Clinical_data_Jana))
  
  for (i in 1:ncol(P_value_table_2)) {
    for (j in 1:nrow(P_value_table_2)) {
      test_temp = kruskal.test(x$U[[4]][,i]~Clinical_data_Jana[,j])
      P_value_table_2[j,i] = test_temp$p.value
    }
  }
  rownames(P_value_table_2) = colnames(Clinical_data_Jana)
  colnames(P_value_table_2) = 1:ncol(P_value_table_2)
  P_value_table_2 = -log10(P_value_table_2)
  pheatmap(P_value_table_2,clustering_method = "ward.D2")
  
  boxplot(x$U[[4]][,5]~Clinical_data_Jana$HER2_Status,outline=F)
  
  return(Tensor_object)
  
}
