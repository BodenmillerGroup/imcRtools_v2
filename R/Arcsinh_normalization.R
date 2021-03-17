#' @rdname Arcsinh_normalization
#' @title Transformation of raw intensity using arcsinh transform 
#'
#' @description This function transform the raw intensity 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cofactor_value numeric value describing the arcsinh cofactor (by default set to 50)


#' @return Returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a new assay slot called "Arcsinh_transformed_intensity"
#'
#' @examples
#' sce = Arcsinh_normalization(sce,cofactor_value = 50)
#' @import SingleCellExperiment
#' @export

Arcsinh_normalization = function(sce,cofactor_value = NULL) {
  
  if (is.null(cofactor_value)) {
    cofactor_value = 50
  }
  
  Transformed_data = asinh(assays(sce)[["Raw_intensity"]]/cofactor_value)
  #Adding the matrix as a new assay slot 
  assay(sce, "Arcsinh_transformed_intensity") <- (Transformed_data)
  
  return(sce)
}

