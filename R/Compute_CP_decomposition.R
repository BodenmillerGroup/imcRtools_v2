#' @rdname Compute_CP_decomposition
#' @title Compute the Canonical Polyadic decomposition of a tensor object and reorder the data in a more understandable way
#'
#' @description Compute the Canonical Polyadic decomposition of a tensor object and reorder the data in a more understandable way
#'
#' @param Interaction_tensor a three-order tensor object,


#' @return Returns a list containing the score matrix (score for each sample in each CP dimension), the list of typical matrices (product of the left and right CP vectors) and a vector with the lambda values
#' @examples
#' CP_decomposition = Compute_CP_decomposition(sce)
#' @import rTensor

Compute_CP_decomposition = function(Interaction_tensor,num_components=5,max_iter = 1000,tol = 1e-6) {
  
  cat("Computing the CP decomposition : \n")
  Resulting_decomposition = cp(Interaction_tensor,num_components = num_components,max_iter = max_iter,tol = tol)
  
  par(las=1,bty="l",mfcol=c(1,2))
  barplot(Resulting_decomposition$lambdas,ylab=c("Lambda value"),xlab="CP dimension",cex.lab=1.3)
  plot(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data),pch=21,bg="orange",
       ylim=range(c(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))),
       xlim= range(c(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))),
       xaxs="i",yaxs='i',xlab="Estimated values",ylab="Observed values",cex.lab=1.3)
  abline(0,1,lwd=2,lty=2)
  R_coef = cor(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))
  legend("topleft",paste("R=",as.character(round(R_coef,2))),bty="n",cex=1.4)
  
  
  cat("Reorganizing the result of the CP decomposition")
  
  List_typical_matrix = c()
  for (k in 1:num_components) {
    List_typical_matrix[[k]] = Resulting_decomposition$U[[1]][,k]%*%t(Resulting_decomposition$U[[2]][,k])
    colnames(List_typical_matrix[[k]]) = 1:ncol(List_typical_matrix[[k]])
    rownames(List_typical_matrix[[k]]) = 1:nrow(List_typical_matrix[[k]])
   }
  
  Score_matrix = Resulting_decomposition$U[[3]]
  colnames(Score_matrix) = paste("CP_dimension_",as.character(1:ncol(Score_matrix)),sep = "")
  return(list(Score_matrix = Score_matrix,List_typical_matrix= List_typical_matrix,Lambda_vector = Resulting_decomposition$lambdas))
  
  
  
}
