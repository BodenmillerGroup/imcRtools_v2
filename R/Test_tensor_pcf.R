Pcf_interaction_tensor = function(sce,r_range=seq(0,300,length.out=50)) {
  

  N_cluster = max(colLabels(sce))
  N_image = max(sce$ImageNumber)
  List_CE_matrix = c()
  
  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_pcf_tensors=foreach(k = 1:N_image) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    X_range = range(sce$Location_Center_X[sce$ImageNumber==k])
    Y_range = range(sce$Location_Center_Y[sce$ImageNumber==k])
    
    Tensor_object_image_temp = rand_tensor(modes = c(N_cluster, N_cluster, length(r_range)), drop = FALSE)
    Tensor_object_image_temp[,,]=1
    
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
          Temp_pcf = pcf(Temp_Kest,method = "b")
          x=Temp_pcf$pcf
        }
        
        if (i==j & npoints(Temp_point_pattern)>10) {
          Temp_Kest = Kest(Temp_point_pattern ,r= r_range,correction = "iso")
          Temp_pcf = pcf(Temp_Kest,method = "b")
          x=Temp_pcf$pcf
          
        }
        
        x[is.na(x)] = 1
        Tensor_object_image_temp[i,j,] = x
        
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
  

  x = cp(Tensor_object,num_components = 5,max_iter = 1000,tol = 1e-8)
  return(Tensor_object)
  
}
