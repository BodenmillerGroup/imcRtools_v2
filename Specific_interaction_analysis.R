Analyze_specific_interaction = function(sce_Jana,i=1,j=2,r_range=seq(0,300,length.out=100),type_function="pcf",N_min_points = 20,bandwidth_size = 5) {
  
  
  if (!type_function%in%c("Besag","Ripley","K","L","pcf")) {
    stop("Please select a correct name of function ! Please choose among the following : Besag,Ripley,K,L,pcf")
  }
  
  N_image = max(sce_Jana$ImageNumber)
  
  cat(paste("Creating parallel backend using"), as.character(metadata(sce_Jana)$N_core), "cores \n")
  registerDoParallel(metadata(sce_Jana)$N_core)
  
  Table_functions = foreach(k = 1:N_image,.combine = cbind) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    X_range = range(sce_Jana$Location_Center_X[sce_Jana$ImageNumber==k])
    Y_range = range(sce_Jana$Location_Center_Y[sce_Jana$ImageNumber==k])
    
    Temp_point_pattern = data.frame(X = sce_Jana$Location_Center_X[sce_Jana$ImageNumber==k & colLabels(sce_Jana)%in%c(i,j) ],
                                        Y = sce_Jana$Location_Center_Y[sce_Jana$ImageNumber==k & colLabels(sce_Jana)%in%c(i,j) ])
    Temp_point_pattern = ppp(Temp_point_pattern$X,Temp_point_pattern$Y,
                                 window = owin(X_range,yrange = Y_range),
                                 marks = factor(colLabels(sce_Jana)[sce_Jana$ImageNumber==k & colLabels(sce_Jana)%in%c(i,j)]))
        
    
    if (sum(marks(Temp_point_pattern)==i)<N_min_points | sum(marks(Temp_point_pattern)==j)<N_min_points){
      #If not enought points to compute the function : put the default values
      x = rep(NA,length(r_range))
    }
    
    if (sum(marks(Temp_point_pattern)==i)>=N_min_points | sum(marks(Temp_point_pattern)==j)>=N_min_points){
      
      if (i!=j) {
        
        if (type_function%in%c("Besag","L")) {
          
          L_function = Kcross(X = Temp_point_pattern,i = i,j=j,r = r_range,correction = "iso")
          L_function = sqrt(L_function$iso/L_function$r) - L_function$r
          x = L_function
        }
        
        if (type_function%in%c("Ripley","K")) {
          
          K_function = Kcross(X = Temp_point_pattern,i = i,j=j,r = r_range,correction = "iso")
          K_function = K_function$iso
          x = K_function$iso
        }
        
        
        if (type_function%in%c("pcf")) {
          
          Pcf = pcf(X = Temp_point_pattern,r = r_range,correction = "iso",bw = bandwidth_size)
          Pcf = Pcf$iso
          x = Pcf
        }
        
      }
      
      if (i==j) {
        
        
        if (type_function%in%c("Besag","L")) {
          
          L_function = Kest(X = Temp_point_pattern,r = r_range,correction = "iso")
          L_function = sqrt(L_function$iso/pi) - L_function$r
          x = L_function
        }
        
        if (type_function%in%c("Ripley","K")) {
          
          K_function = Kest(X = Temp_point_pattern,r = r_range,correction = "iso")
          K_function = K_function$iso
          x = K_function$iso
        }
        
        
        if (type_function%in%c("pcf")) {
          Pcf = pcf(X = Temp_point_pattern,r = r_range,correction = "iso",bw = bandwidth_size)
          Pcf = Pcf$iso
          x = Pcf
        }
        
      }
      
    }
    
    
    x
    
  }
  
  colnames(Table_functions) = 1:N_image
  Table_functions = Table_functions[-1,]
  #Visulisation of the functions
  y_function_range = round(quantile(as.numeric(Table_functions),c(0.001,0.999),na.rm=T))
  par(las=1,bty="l")
  plot(NULL,ylim=y_function_range,xlim=range(r_range),xaxs='i',yaxs='i',ylab="Function value",xlab="r")
  
  for (k in 1:N_image) {
    points(r_range[-1],Table_functions[,k],type="l",)
  }
  
}

  library(FactoMineR)
  Table_functions_filled = Table_functions
  Table_functions_filled[is.na(Table_functions_filled)]=1

  PCA_functions = PCA(t(Table_functions_filled),scale.unit = F,graph = F,ncp = 10)
  barplot(PCA_functions$eig[1:10,2])
  plot(PCA_functions$var$coord[,2],type="l")
  plot(PCA_functions$ind$coord[,c(1,3)])
  kruskal.test(PCA_functions$ind$coord[,1]~Clinical_data_Jana$Microinvasion)



  library("HDPenReg")
  fused_lasso_regression_cv = EMcvfusedlasso(y = as.numeric(Clinical_data_Jana$Microinvasion),X = t(Table_functions),model = "logistic",lambda1 = seq(0,15,length.out=30),lambda2 = seq(0,15,length.out=30))
  fused_lasso_regression = EMfusedlasso(y = as.numeric(Clinical_data_Jana$Microinvasion),X = t(Table_functions),model = "logistic",lambda1 = fused_lasso_regression_cv$lambda.optimal[1],lambda2 = fused_lasso_regression_cv$lambda.optimal[2])
  plot(fused_lasso_regression$coefficient,type="o",pch=21,bg="orange")
  
  
